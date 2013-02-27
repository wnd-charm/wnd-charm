#include "SharedImageMatrix.h"
#include "cmatrix.h"

#include <unistd.h> // sysconf(), page_size
#include <errno.h>
#include <sys/types.h> // for dev_t, ino_t
#include <sys/stat.h>  // fstat, stat
#include <sys/mman.h>  // mmap, shm_open
#include "digest/md5.h"
#include "b64/encode.h"

/* global variable */
extern int verbosity;

static inline std::string string_format(const std::string &fmt, ...) {
    int size = 256;
    std::string str;
    va_list ap;
    while (1) {
        str.resize(size);
        va_start(ap, fmt);
        int n = vsnprintf((char *)str.c_str(), size, fmt.c_str(), ap);
        va_end(ap);
        if (n > -1 && n < size) {
            str.resize(n);
            return str;
        }
        if (n > -1)
            size = n + 1;
        else
            size *= 2;
    }
    return str;
}

struct shmem_data {
	uint32_t width, height;
	enum ColorModes ColorMode;
	uint8_t bits;
};

// storage and initialization for the object statics
size_t SharedImageMatrix::shmem_page_size = sysconf(_SC_PAGE_SIZE);
bool SharedImageMatrix::never_read = false;
bool SharedImageMatrix::disable_destructor_cache_cleanup = false;
pid_t SharedImageMatrix::PID = getpid();
std::string SharedImageMatrix::PID_string = string_format ("%lld",(long long)SharedImageMatrix::PID);



void SharedImageMatrix::SetShmemName ( ) {
	MD5::MD5 digest;
	uint8_t Message_Digest[MD5::HashSize];

	assert (!(cached_source.empty() && operation.empty()) && "Attempt to establish a shared memory name without a cached_source or operation field");
	if (!cached_source.empty()) digest.Update( (uint8_t *)(cached_source.data()), cached_source.length() );
	if (!operation.empty()) digest.Update( (uint8_t *)(operation.data()), operation.length() );
	digest.Result(Message_Digest);

	shmem_name = "/wndchrm";
	// base64-encode without newlines or padding (last two booleans)
	base64::encoder().encode ((char *)Message_Digest, MD5::HashSize, shmem_name, false, false);
	std::cout <<         "         SharedImageMatrix::SetShmemName: " << shmem_name << std::endl;
	
}

// This is a helper class method to calculate offsets into shared memory.
// Note that this is a class method declared as static in the header.
//   It is not declared static here because static in a .cpp means something else entirely.
const size_t SharedImageMatrix::calc_shmem_size (const unsigned int w, const unsigned int h, const enum ColorModes ColorMode, size_t &clr_plane_offset, size_t &shmem_data_offset) {
	size_t new_mat_size = w * h;
	size_t new_shmem_size = new_mat_size * sizeof (double);
	// Expand the size to be a multiple of the page size.
	new_shmem_size = ( (1 + (new_shmem_size / shmem_page_size)) * shmem_page_size );
	// The color plane starts at a page boundary.
	clr_plane_offset = new_shmem_size;
	if (ColorMode != cmGRAY) new_shmem_size += new_mat_size * sizeof (HSVcolor);
	// the shmem_data_t struct is stored at the end to preserve page-size memory alignment for Eigen.
	new_shmem_size += sizeof (shmem_data);
	// Expand the total size to be a multiple of the page size.
	new_shmem_size = ( (1 + (new_shmem_size / shmem_page_size)) * shmem_page_size );
	shmem_data_offset = new_shmem_size - sizeof (shmem_data);
	return (new_shmem_size);
}


// override of parent class allocate method to use shared memory.
// This method sets up the shared memory to accomodate the matrixes, and maps the Eigen maps to use it for their data.
// It is also able to resize an existing shared memory object.
void SharedImageMatrix::allocate (unsigned int w, unsigned int h) {
 	std::cout << "-------- called SharedImageMatrix::allocate (" << w << "," << h << ") on " << shmem_name << std::endl;

	assert (shmem_page_size && "Memory page size is undefined!");
	assert (shmem_fd > -1 && "Shared memory file descriptor is invalid!");
	// Additional checks:
	//   There cannot be any allocate calls on a cached matrix - the allocation is done by fromCache()
	//   A matrix was read from cache if shmem_name == cached_source (was_cached is true)
	assert (!was_cached && "Attempt to call allocate on a cached object!");
	

	// calculate the size of the required shared memory block
	size_t new_shmem_size, clr_plane_offset, shmem_data_offset;
	new_shmem_size = calc_shmem_size (w, h, ColorMode, clr_plane_offset, shmem_data_offset);
	std::cout << " shmem_size: " << new_shmem_size << " pages: " << new_shmem_size / shmem_page_size << std::endl;

	// Map shared memory object for writing
	// Unmap any pre-existing memory and re-map
	// FIXME: It may be good to only resize expanding memory. Shrinking memory should only need a remap of the Eigen Map.
	if (mmap_ptr != MAP_FAILED) {
		if ( munmap (mmap_ptr, shmem_size) ) {
			std::cout << "munmap error: " << strerror(errno) << std::endl;
			exit (-1);
		}
		// OS X forces us to delete the segment entirely because only one call to ftruncate can be made on each segment.
		// FIXME: It may be good to escape this on non-OS X systems.
		if (shmem_fd > -1) close (shmem_fd);
		shmem_fd = -1;
		shm_unlink (shmem_name.c_str());

		shmem_fd = shm_open(shmem_name.c_str(), O_RDWR | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
		if (shmem_fd < 0) {
			std::cout << "shm_open error: " << strerror(errno) << std::endl;
			exit (-1);
		}
		shmem_size = 0;
	}

	// size the file backing the memory object
	if (ftruncate(shmem_fd, new_shmem_size) == -1) {
		std::cout << "ftruncate error: " << strerror(errno) << std::endl;
		exit (-1);
	}

	mmap_ptr = (byte *)mmap(NULL, new_shmem_size, PROT_READ | PROT_WRITE, MAP_SHARED, shmem_fd, (off_t) 0);
	if (mmap_ptr == MAP_FAILED) {
		std::cout << "mmap error: " << strerror(errno) << std::endl;
		exit (-1);
	}
	shmem_size = new_shmem_size;
	// remap the data for the object to use the mmap_ptr.
	remap_pix_plane ( (double *)mmap_ptr, w, h);
	if (ColorMode != cmGRAY) remap_clr_plane ((HSVcolor *)(mmap_ptr + clr_plane_offset), w, h);
			
	
	// If this memory gets read from cache, we wont know the size of the matrix,
	// so we store the shmem_data at the end of the shmem to tell us how to read it later.
	// However, to ensure that read-validation passes only after a successful call to Cache(),
	// we set this region of memory to 0.
	memset((mmap_ptr + shmem_data_offset), 0, sizeof (shmem_data));
}

// This returns a CacheStatus enum explaining what happened with the cache
// The operation field and the cached_source field determine the object's cache name (shmem_name)
// In the returned enum, cache_unknown means there was an error
// cache_read means the SharedImageMatrix was read in from cache.
// cache_write means the SharedImageMatrix has exclusive write access.
//   The object is empty at this point.
// The chached_source_in and operation_in parameters are optional.
//   If not empty, they will replace their private couterparts before calculating the shmem_name

void SharedImageMatrix::fromCache (const std::string cached_source_in, const std::string operation_in) {	
	// Make sure the shared memory has a name
	if (!cached_source_in.empty()) cached_source = cached_source_in;
	if (!operation_in.empty()) operation = operation_in;
std::cout <<         "-------- called SharedImageMatrix::fromCache on [" << cached_source << "]->[" << operation << "]" <<std::endl;;
	SetShmemName();

	// Shared memory access with concurrent processes
	// Only one process should be able to write the shared memory.  All others must wait if necessary and then read only.
	// There are lots of issues revolving around concurrency using POSIX "real time" methods. TL; DNR: Don't bother to try it.
	// The POSIX standard says that fcntl locks on shared memory are "unspecified".
	// Although they do currently work on linux, they do not currently work on OS X. The status of BSD, Solaris, etc. is unknown.
	// POSIX named/shared semaphores are unusable in a real world containing uncaught and uncatchable signals for two reasons:
	//      1) A process that dies with a locked semaphore will cause all other processes waiting for that semaphore
	//       to wait forever.  Uncaught signals cause deadlocks.
	//      2) Uncaught signals prevent semaphores from being opened/created in a known state.  A process may assume
	//       that its creating a semaphore in a known state, but it may be re-attaching to a pre-existing semaphore
	//       locked by a dead process.  No way to know which.
	//    2.5) The standard is ambiguous with regards to unlinking semaphores.  As worded, a sem_unlink() will cause
	//       subsequent sem_open to attach to a new semaphore even though the unlinked semaphore won't actually be destroyed
	//       until a later call to sem_close (directly or indirectly b/c of process termination).  This suggests that we may
	//       have two *different* named semaphores on the same system with the exact same name!
	//       From <http://pubs.opengroup.org/onlinepubs/7908799/xsh/sem_unlink.html>:
	//         "Calls to sem_open() to re-create or re-connect to the semaphore refer to a new semaphore after sem_unlink() is called"
	//       This is clearly bogus, especially if its actually implemented this way.
	//    2.8) There is no mechanism by which to list and clean up "stale" POSIX semaphores. ipcs/ipcrm only works with SysV semaphores.
	// Anonymous/unnamed POSIX semaphores sound nice, but they can only be used for related processes (not unrelated),
	//    and in any case, they are unimplemented on OS X.
	// SysV semaphores are well supported can revert/undo their state when a process with a locked semaphore terminates - they remain a viable option.
	// However, the SysV interface is deprecated in favor of the unworkable POSIX named semaphore interface.  Progress!
	// Yet another option is to use the POSIX thread API to get a rwlock on the shared memory, storing the lock structure in the shared memory,
	// and ensuring that it is only initialized once using the once call.  Unfortunately, this is implementation-dependent because support for
	// the shared flag for an rwlock is optional, and is not supported on OS X.
	// OS X in other words leaves with a single workable option: Use a file lock (or SysV semaphores).

	// Use WORMfile to determine if we should read or write. This creates a lockfile with the same name as our shared memory in the /tmp/ directory.
	// Since this is a Write-Once-Read-Many file, asking for a write lock may result in a read lock. A stale (empty) file will result in a write lock.
	// Asking it to wait until some kind of lock can be made will ensure that we only need to deal with immediate reads or writes at this point.
	lock_file.path = std::string ("/tmp").append (shmem_name);
	// if the never_read flag is set, we unlink the lock_file before opening it, so we always get a write-lock.
	if (never_read) unlink (lock_file.path.c_str());
	// Use the reopen method since the constructor was never called. Send optional parameters for not-read-only, wait for lock.
	lock_file.reopen (false, true);
		
	// This section is responsible for doing an immediate read or write.
	// The lock is already taken care of if it had an error.
	// Otherwise, we have either a readlock or a write lock, but no open shmem.
	switch (lock_file.status) {
		case WORMfile::WORM_RD: {
std::cout << "cache_read" << std::endl;
			assert (!never_read && "SharedImageMatrix Class is set to never read, but the lockfile still exists after unlinking it!");
			// This is an immediate read, and we already have a read lock
			shmem_fd = shm_open(shmem_name.c_str(), O_RDONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
			if (shmem_fd < 0) {
				error_str = std::string ("shm_open error when reading: ") + strerror(errno);
				break;
			}

			// Get the shared memory size from the file size.
			struct stat the_stat;
			if ( fstat (shmem_fd, &the_stat) == -1 ) {
				error_str = std::string ("fstat error determining length of existing shmem: ") + strerror(errno);
				break;
			}
			shmem_size = (size_t)the_stat.st_size;
std::cout << "shmem_size: " << shmem_size << std::endl;
			if (! shmem_size > sizeof (shmem_data) ) {
				error_str = std::string ("size of existing shmem is 0");
				break;
			}

			// mmap the shared memory segment
			mmap_ptr = (byte *)mmap(NULL, shmem_size, PROT_READ, MAP_SHARED, shmem_fd, (off_t)0);
			if (mmap_ptr == MAP_FAILED) {
				error_str = std::string ("mmap error when mapping existing shmem: ") + strerror(errno);
				break;
			}
			// get the shmem_data object stored at the end of shared memory
			// This will tell us how to read the matrixes and reconstruct the SharedImageMatrix
			size_t shmem_data_offset = shmem_size - sizeof (shmem_data);
			shmem_data *stored_shmem_data = (shmem_data *)(mmap_ptr + shmem_data_offset);

			// ensure that the memory size is correct.
			size_t stored_mat_shmem_size, stored_clr_plane_offset, stored_shmem_data_offset;
			stored_mat_shmem_size = calc_shmem_size (stored_shmem_data->width, stored_shmem_data->height,
				stored_shmem_data->ColorMode, stored_clr_plane_offset, stored_shmem_data_offset
			);
			if (stored_mat_shmem_size < shmem_size) {
				error_str = string_format ("error when mapping existing shmem: stored data requires %lu bytes, but shmem size is %lu bytes",
					(unsigned long)stored_mat_shmem_size, (unsigned long)shmem_size);
				break;
			}
			
			// Looks like we have a valid matrix stored, so create the cached result.
			// remap the data for the object to use the mmap_ptr, keeping the rest of the object where it was.
			remap_pix_plane ( (double *)mmap_ptr, stored_shmem_data->width, stored_shmem_data->height);
			if (ColorMode != cmGRAY) remap_clr_plane ((HSVcolor *)(mmap_ptr + stored_clr_plane_offset), stored_shmem_data->width, stored_shmem_data->height);
			if (width != stored_shmem_data->width || height != stored_shmem_data->height) {
				error_str = string_format ("error when mapping existing shmem: recovered w,h (%u, %u) doesn't match that in shmem (%u, %u)",
					(unsigned int)width, (unsigned int)height,
					(unsigned int)stored_shmem_data->width, (unsigned int)stored_shmem_data->height);
				break;
			}
			ColorMode = stored_shmem_data->ColorMode;
			bits = stored_shmem_data->bits;
			
			// The pixels are read-only.
			WriteablePixelsFinish ();
			if (ColorMode != cmGRAY) WriteableColorsFinish ();
			
			// Indicate that the object came from cache, and currently has no operation
			cached_source = shmem_name;
			operation = "";
			was_cached = true;
			cache_status = csREAD;
			
			// close the lockfile
			// Since this is a WORM file, we don't need to maintain an active readlock
			lock_file.finish();


std::cout << "recovered size: " << width << "," << height << std::endl;
			error_str = "";
			return;
		}
		break;
		
		case WORMfile::WORM_WR:
std::cout << "cache_write" << std::endl;
			// Since we have a write-lock, we open the shmem for writing, creation, and truncation.
			shmem_fd = shm_open(shmem_name.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
			if (shmem_fd < 0) {
				error_str = string_format ("shm_open error when writing: %s (%d)", strerror(errno), errno);
				break;
			}
			was_cached = false;
			error_str = "";
			cache_status = csWRITE;
			// Do not release the write-lock or close the lockfile at this point - only after we're done writing by calling Cache().
			return;
		break;
		
		default:
			error_str = string_format ("error getting lock on lockfile '%s': %s", lock_file.path.c_str(), strerror(lock_file.status_errno));
			was_cached = false;
		break;
	}

	// Not having returned at this point means we had an error.
	cache_status = csERROR;
	// Leave a clean slate for a subsequent attempt: no lock file, no shmem.
	// unmap if there is an mmap_ptr
	if (mmap_ptr != MAP_FAILED) munmap (mmap_ptr, shmem_size);
	mmap_ptr = (byte *)MAP_FAILED;
	// close shmem_fd if its open
	if (shmem_fd > -1) close (shmem_fd);
	shmem_fd = -1;
	// unlink the shmem file
	shm_unlink (shmem_name.c_str());
	// close and unlink the lock
	lock_file.finish();
	unlink (lock_file.path.c_str());
	lock_file = WORMfile();

	return;
}

// This method must be called to finalize an object that will be stored in the cache.
// It must not be called on objects retrieved from the cache.
void SharedImageMatrix::Cache ( ) {

	assert (!was_cached && "Called Cache on an object that was read from cache.");
std::cout <<         "-------- called SharedImageMatrix::Cache on [" << cached_source << "]->[" << operation << "]" <<std::endl;;
	SetShmemName();
	
	// make the shmem_data region valid.
	shmem_data *stored_shmem_data = (shmem_data *)(mmap_ptr + shmem_size - sizeof(shmem_data));
	stored_shmem_data->width = width;
	stored_shmem_data->height = height;
	stored_shmem_data->ColorMode = ColorMode;
	stored_shmem_data->bits = bits;
	
	// make the lockfile valid (zero-length files count as "stale")
	write (lock_file.fd(),  PID_string.data(), PID_string.length());
	// this closes the lockfile, releasing all locks.
	lock_file.finish();
	
	// This object should now be indistinguishable from a cached object
	// This is not strictly true because the memory is still mapped with write permissions.
	// However, attempting a write by calling WriteablePixels or allocate will result in run-time assertions.
	WriteablePixelsFinish ();
	if (ColorMode != cmGRAY) WriteableColorsFinish ();
	cached_source = shmem_name;
	operation = "";
	was_cached = true;
}


// load an image of any supported format
// This attempts to read the file from cache, and calls the superclass to do the work if we need to read it in.
int SharedImageMatrix::OpenImage(char *image_file_name,
		int downsample, rect *bounding_rect,
		double mean, double stddev) {

	std::string open_operation;
	int fildes = open (image_file_name, O_RDONLY);
	if (fildes > -1) {
		struct stat the_stat;
		fstat (fildes, &the_stat);
		close (fildes);
		open_operation = string_format ("Open_%lld_%lld", (long long)the_stat.st_dev, (long long)the_stat.st_ino);
	} else {
		open_operation = string_format ("Open_%s", image_file_name);
	}

	if (bounding_rect && bounding_rect->x >= 0)
		open_operation.append ( string_format ("_Sub_%d_%d_%d_%d",
			bounding_rect->x, bounding_rect->y,
			bounding_rect->x+bounding_rect->w-1, bounding_rect->y+bounding_rect->h-1
		));
	if (downsample>0 && downsample<100)  /* downsample by a given factor */
		open_operation.append ( string_format ("_DS_%lf_%lf",((double)downsample)/100.0,((double)downsample)/100.0) );
	if (mean>0)  /* normalize to a given mean and standard deviation */
		open_operation.append ( string_format ("_Nstd_%lf_%lf",mean,stddev) );

	// We're setting an empty cached_source since we don't have a shared memory name to begin with.
	fromCache ("", open_operation);
	
	if (cache_status == csWRITE) {
		int ret = ImageMatrix::OpenImage (image_file_name,downsample,bounding_rect,mean,stddev);
		Cache();
		return (ret);
	} else if (cache_status == csREAD) {
		return (1);
	} else {
		std::cerr << "Errors while recovering cache:" << error_str << std::endl;
		exit (-1);
	}

}


// This is a general transform method that returns a new image matrix by applying the specified transform.
// Of course this returns a new SharedImageMatrix as an ImageMatrix pointer.
// We have to declare the virtual method this way to make the override work.
ImageMatrix &SharedImageMatrix::transform (const ImageTransform *transform) const {
	SharedImageMatrix *matrix_OUT = new SharedImageMatrix;
	matrix_OUT->fromCache (this->GetShmemName(), transform->name);

	if (matrix_OUT->Status() == csWRITE) {
		transform->execute (this, matrix_OUT);
		matrix_OUT->Cache();
		return (*matrix_OUT);
	} else if (matrix_OUT->Status() == csREAD) {
		return (*matrix_OUT);
	} else {
		std::cerr << "Errors while recovering cache:" << matrix_OUT->Error() << std::endl;
		exit (-1);
	}
}


// The destructor will unlink the shared memory unless DisableDestructorCacheCleanup(true) class method has been called.
SharedImageMatrix::~SharedImageMatrix () {
std::cout << "SharedImageMatrix DESTRUCTOR for " << shmem_name << std::endl;

	WriteablePixelsFinish();
	remap_pix_plane (NULL, 0, 0);

	if (ColorMode != cmGRAY) WriteableColorsFinish();
	remap_clr_plane (NULL, 0, 0);

	if (mmap_ptr != MAP_FAILED) munmap (mmap_ptr, shmem_size);
	mmap_ptr = (byte *)MAP_FAILED;
	// close shmem_fd if its open
	if (shmem_fd > -1) close (shmem_fd);
	shmem_fd = -1;
	// close the lockfile
	lock_file.finish();
	// unlink the shmem file and the lockfile
	if (!disable_destructor_cache_cleanup) {
std::cout << "    unlinking POSIX shared memory " << shmem_name << " and lockfile " << lock_file.path << std::endl;
		shm_unlink (shmem_name.c_str());
		unlink (lock_file.path.c_str());
	}
	cache_status = csUNKNOWN;

}
