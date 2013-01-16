#include <iostream> // used for debug output from instantiator methods
#include <iomanip>
#include <cmath>
#include <fcntl.h>
#include <unistd.h> // sysconf(), page_size
#include <errno.h>
#include <sys/types.h> // for dev_t, ino_t
#include <sys/stat.h>  // fstat, stat
#include <sys/mman.h>  // mmap, shm_open
#include <semaphore.h> // semaphores

// #include "digest/sha1.h"
#include "digest/md5.h"
#include "b64/encode.h"
#include "cmatrix.h"
#include "FeatureTransforms.h"
#include "colors/FuzzyCalc.h" // for definition of compiler constant COLORS_NUM
#include "transforms/fft/bcb_fftw3/fftw3.h"

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
	uint8_t ColorMode;
	uint8_t bits;
};

// storage and initialization for the object static storing the page size.
size_t SharedImageMatrix::shmem_page_size = sysconf(_SC_PAGE_SIZE);


void SharedImageMatrix::SetShmemName ( ) {
	MD5::MD5 digest;
	uint8_t Message_Digest[MD5::HashSize];

	assert ((cached_source.length() || operation.length()) && "Attempt to establish a shared memory name without a cached_source or operation field");
	if (cached_source.length()) digest.Update( (uint8_t *)(cached_source.data()), cached_source.length() );
	if (operation.length()) digest.Update( (uint8_t *)(operation.data()), operation.length() );
	digest.Result(Message_Digest);

	shmem_name = "/wndchrm";
	// base64-encode without newlines or padding (last two booleans)
	base64::encoder().encode ((char *)Message_Digest, MD5::HashSize, shmem_name, false, false);
	std::cout <<         "         SharedImageMatrix::SetShmemName: " << shmem_name << std::endl;
	
}

const size_t SharedImageMatrix::calc_shmem_size (const unsigned int w, const unsigned int h, size_t &clr_plane_offset, size_t &shmem_data_offset) {
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


void SharedImageMatrix::allocate (unsigned int w, unsigned int h) {
 	std::cout << "-------- called SharedImageMatrix::allocate (" << w << "," << h << ") on " << shmem_name << std::endl;

	assert (shmem_page_size && "Memory page size is undefined!");
	assert (shmem_fd > -1 && "Shared memory file descriptor is invalid!");
	// Additional checks:
	//   There cannot be any allocate calls on a cached matrix - the allocation is done by fromCache()
	//   A matrix was read from cache if shmem_name == cached_source (was_cached is true)
	assert (!was_cached && "Attempt to call allocate on a cached object!");
	

	// calculate the size of the required shared memory block
	size_t new_mat_size = w * h;
	size_t new_shmem_size = new_mat_size * sizeof (double);
	// Expand the size to be a multiple of the page size.
	new_shmem_size = ( (1 + (new_shmem_size / shmem_page_size)) * shmem_page_size );
	// The color plane starts at a page boundary.
	size_t clr_plane_offset = new_shmem_size;
	if (ColorMode != cmGRAY) new_shmem_size += new_mat_size * sizeof (HSVcolor);
	// the shmem_data_t struct is stored at the end to preserve page-size memory alignment for Eigen.
	new_shmem_size += sizeof (shmem_data);
	// Expand the total size to be a multiple of the page size.
	new_shmem_size = ( (1 + (new_shmem_size / shmem_page_size)) * shmem_page_size );
	size_t shmem_data_offset = new_shmem_size - sizeof (shmem_data);
	std::cout << " shmem_size: " << new_shmem_size << " pages: " << new_shmem_size / shmem_page_size << std::endl;

	// Map shared memory object for writing
	// Unmap any pre-existing memory and re-map
	if (mmap_ptr != MAP_FAILED) {
		if ( munmap (mmap_ptr, shmem_size) ) {
			std::cout << "munmap error: " << strerror(errno) << std::endl;
			exit (-1);
		}
		// OS X actually forces us to delete the segment entirely because only one call to ftruncate can be made on each segment.
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
	// remap the data for the object to use the mmap_ptr, keeping the rest of the object where it was.
	remap_pix_plane ( (double *)mmap_ptr, w, h);

	if (ColorMode != cmGRAY) remap_clr_plane ((HSVcolor *)(mmap_ptr + clr_plane_offset), w, h);
			
	
	// If this memory gets read from cache, we wont know the size of the matrix,
	// so we have to store the shmem_data at the end of the shmem
	shmem_data *stored_shmem_data = (shmem_data *)(mmap_ptr + shmem_data_offset);
	stored_shmem_data->width = width;
	stored_shmem_data->height = height;
	stored_shmem_data->ColorMode = ColorMode;
	stored_shmem_data->bits = bits;
}

// This returns a CacheStatus enum explaining what happened with the cache
// The operation field and the cached_source field determine the object's cache name (shmem_name)
// cache_unknown means there was an error
// cache_read means the SharedImageMatrix was read in from cache.
// cache_write means the SharedImageMatrix has exclusive write access.
//   The object is empty at this point.

const CacheStatus SharedImageMatrix::fromCache () {
	CacheStatus cache_status = cache_unknown;
	std::string error_str;
	
std::cout <<         "-------- called SharedImageMatrix::fromCache on [" << cached_source << "]->[" << operation << "]" <<std::endl;;
	// Make sure the shared memory has a name
	SetShmemName();

	// Shared memory access with concurrent processes
	// Only one process should be able to write the shared memory.  All others must wait and/or read only.
	// The POSIX standard says that fcntl locks on shared memory are "unspecified".
	// Although they do currently work on linux, they do not currently work on OS X. The status of BSD, Solaris, etc. is unknown.
	// The POSIX standard may or nay not support this behavior in the future.
	// Instead, we will use POSIX semaphores to ensure exclusive write access.
	// The semaphore has the same name as the shared memory segment, and it has an associated value.
	// For our purposes, the semaphore can be in two states: value > 0 (shmem is readable) or value == 0 (shmem is being updated by another process).
	// We have to be careful not to leak semaphores, so we have to create->unlink->trywait.
	// The semaphore will be created in an unlocked state (value = 1), immediately unlinked, and then call try_wait on it.
	// if the try_wait call succeeds, it also locks the semaphore, indicating we have exclusive read or write access.
	// if it does not succeed and gives an EAGAIN error, it means we could not get a lock because of write operations, so we must call sem_wait, then read.
	// if the process quits at any time, the postponed sem_unlink call will take effect, so a subsequent call to sem_open will crate an unlocked semaphore.
	// The shmem lifetime events are then:
	// locked semaphore -> Create shmem -> mature/write shmem -> post semaphore -> readable shmem -> close semaphore -> unlink shmem
	// After a successful sem_open (O_CREAT)/try_wait(), an EEXIST on shm_open (O_RDWR|O_CREAT|O_EXCL) tells us the shmem is ready for reading, so we:
	//   shm_open (O_RDONLY) -> mmap (PROT_READ) ...
	// A successful shm_open (O_RDWR|O_CREAT|O_EXCL) tells us that we just created the shmem, and it is ready for exclusive writing.
	// Any other error indicates that "bad things happened"
	// Regardless of the nature of the "bad things" that did happen, the response is to unlink the sem and the shm and start all over.

	// This first section results in a locked semaphore and possibly a shmem_fd ready for writing.
	// If there were error with the semaphore or opening shmem, the semaphore is cleaned up.
	// open or create an unlocked semaphore
	shmem_sem = sem_open(shmem_name.c_str(), O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH, 1);
	int shmem_sem_locked = EAGAIN;
	if (shmem_sem == SEM_FAILED) {
		// This is an error - something bad happened
		error_str = std::string ("sem_open error: ") + strerror(errno);
	} else {
		// The sem_unlink call is postponed until we call sem_close or we terminate.
		// Calling it here ensures that it will be called if we terminate.
		sem_unlink (shmem_name.c_str());
		// Try to lock the semaphore.
		shmem_sem_locked = sem_trywait (shmem_sem);
		if (shmem_sem_locked == EAGAIN) {
			// another process has an exclusive semaphore lock, so we wait for it to finish
			if (sem_wait (shmem_sem) == 0) {
				// after a successful wait, the shmem should be ready for reading.
				cache_status = cache_read;				
			} else {
				// If sem_wait returned an error, then something bad happened.
				error_str = std::string ("sem_wait error: ") + strerror(errno);
				sem_close (shmem_sem);
				shmem_sem = SEM_FAILED;
			}

		} else if (shmem_sem_locked == 0) {
			// We have an exclusive lock on the semaphore, either because shmem is ready for reading or because it doesn't exist yet.
			// First, try to create/open shmem exclusively.
			// To force a write every time, unlink the shmem at this point:
			//shm_unlink (shmem_name.c_str());
			shmem_fd = shm_open(shmem_name.c_str(), O_RDWR | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
			if (shmem_fd > -1) {
				// We have exclusive write access because we just created the shmem.
				// The allocate() method takes care of using the shared memory at shmem_fd.
				errno = 0;
				cache_status = cache_write;
			} else if (shmem_fd < 0 && errno == EEXIST) {
				// not an error: shmem_fd is invalid because the shmem exists
				// This means we should be ready to read.
				errno = 0;
				cache_status = cache_read;
			} else {
				// This is an error - something bad happened
				error_str = std::string ("shm_open error when opening/creating: ") + strerror(errno);
				// Clean up the semaphore.
				sem_close (shmem_sem);
				shmem_sem = SEM_FAILED;
			}

		} else {
			// Some bad things happened trying to lock the semaphore.
			error_str = std::string ("sem_trywait error: ") + strerror(errno);
			// Clean up the semaphore.
			sem_close (shmem_sem);
			shmem_sem = SEM_FAILED;
		}
	}
	
	// This section is responsible for doing an immediate read or write.
	// The semaphore is already taken care of if it had an error.
	// Otherwise, the semaphore is locked.
	switch (cache_status) {
		case cache_unknown:
			// don't deal with errors yet.
			break;
		break;
		
		case cache_read: {
std::cout << "cache_read" << std::endl;
			// This is an immediate read.
			// We have a locked semaphore, but no open shmem file
			shmem_fd = shm_open(shmem_name.c_str(), O_RDONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
			if (!shmem_fd < 0) {
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
			// get the matrix sizes from the shmem_data object stored at the end of shared memory
			size_t shmem_data_offset = shmem_size - sizeof (shmem_data);
			shmem_data *stored_shmem_data = (shmem_data *)(mmap_ptr + shmem_data_offset);

			// ensure that the memory size is correct.
			size_t stored_mat_shmem_size, stored_clr_plane_offset, stored_shmem_data_offset;
			stored_mat_shmem_size = calc_shmem_size (stored_shmem_data->width, stored_shmem_data->height,
				stored_clr_plane_offset, stored_shmem_data_offset);
			if (stored_mat_shmem_size != shmem_size) {
				error_str = string_format ("error when mapping existing shmem: stored data requires %lu bytes, but shmem size id %lu bytes",
					(unsigned long)stored_mat_shmem_size, (unsigned long)shmem_size);
				break;
			}
			
			// Looks like we have a valid matrix stored, so create the cached result.
			ColorMode = static_cast<ColorModes> (stored_shmem_data->ColorMode);
			// remap the data for the object to use the mmap_ptr, keeping the rest of the object where it was.
			remap_pix_plane ( (double *)mmap_ptr, stored_shmem_data->width, stored_shmem_data->height);
			if (ColorMode != cmGRAY) remap_clr_plane ((HSVcolor *)(mmap_ptr + stored_clr_plane_offset), stored_shmem_data->width, stored_shmem_data->height);
			bits = stored_shmem_data->bits;
			if (width != stored_shmem_data->width || height != stored_shmem_data->height) {
				error_str = string_format ("error when mapping existing shmem: recovered w,h (%u, %u) doesn't match that in shmem (%u, %u)",
					(unsigned int)width, (unsigned int)height,
					(unsigned int)stored_shmem_data->width, (unsigned int)stored_shmem_data->height);
				break;
			}
			
			// The recovered object is as verified as we can manage, so release and unlink the semaphore.
			sem_post (shmem_sem);
			sem_close (shmem_sem);
			shmem_sem = SEM_FAILED;
			
			// The pixels are read-only.
			WriteablePixelsFinish ();
			if (ColorMode != cmGRAY) WriteableColorsFinish ();
			
			// Indicate that the object came from cache, and currently has no operation
			cached_source = shmem_name;
			operation = "";
			was_cached = true;

std::cout << "recovered size: " << width << "," << height << std::endl;
			return (cache_status);
		}
		break;
		
		case cache_write:
std::cout << "cache_write" << std::endl;
			return (cache_status);
		break;
		
		case cache_wait:
		// this can't happen here.
		break;
	}

	// Not having returned at this point means we had an error while reading.
	// We have a locked semaphore, possibly an open shmem_fd, possibly an mmap_ptr to clean up.
	// unmap if there is an mmap_ptr
	if (mmap_ptr != MAP_FAILED) munmap (mmap_ptr, shmem_size);
	mmap_ptr = (byte *)MAP_FAILED;
	// close shmem_fd if its open
	if (shmem_fd > -1) close (shmem_fd);
	shmem_fd = -1;
	// unlink the shmem file
	shm_unlink (shmem_name.c_str());
	// close and unlink the semaphore
	if (shmem_sem != SEM_FAILED) sem_close (shmem_sem);
	shmem_sem = SEM_FAILED;

	// report the error and exit.  Since we have a clean slate, may want to try again from the beginning.
	std::cerr << "Errors while recovering cache (cout):" << std::endl;
	std::cerr << error_str << std::endl;
	exit (-1);
	return (cache_status);
}

void SharedImageMatrix::Cache ( ) {

std::cout <<         "-------- called SharedImageMatrix::Cache on [" << cached_source << "]->[" << operation << "]" <<std::endl;;
	SetShmemName();

	// Since shared memory was set up in the fromCache call,
	// all that's left to do here is make sure its finalized and ready to use by others

}



int SharedImageMatrix::OpenImage(char *image_file_name,            // load an image of any supported format
		int downsample, rect *bounding_rect,
		double mean, double stddev) {

	int fildes = open (image_file_name, O_RDONLY);
	if (fildes > -1) {
		struct stat the_stat;
		fstat (fildes, &the_stat);
		close (fildes);
		operation = string_format ("Open_%lld_%lld", (long long)the_stat.st_dev, (long long)the_stat.st_ino);
	} else {
		operation = string_format ("Open_%s", image_file_name);
	}

	if (bounding_rect && bounding_rect->x >= 0)
		operation.append ( string_format ("_Sub_%d_%d_%d_%d",
			bounding_rect->x, bounding_rect->y,
			bounding_rect->x+bounding_rect->w-1, bounding_rect->y+bounding_rect->h-1
		));
	if (downsample>0 && downsample<100)  /* downsample by a given factor */
		operation.append ( string_format ("_DS_%lf_%lf",((double)downsample)/100.0,((double)downsample)/100.0) );
	if (mean>0)  /* normalize to a given mean and standard deviation */
		operation.append ( string_format ("_Nstd_%lf_%lf",mean,stddev) );

	CacheStatus cache_status = fromCache ( );
	
	if (cache_status == cache_write) {
		int ret = ImageMatrix::OpenImage (image_file_name,downsample,bounding_rect,mean,stddev);
		Cache();
		return (ret);
	} else if (cache_status == cache_read) {
		return (1);
	} else {
		std::cerr << "Errors while recovering cache" << std::endl;
		exit (-1);
	}

}

SharedImageMatrix::~SharedImageMatrix () {
std::cout << "SharedImageMatrix DESTRUCTOR for " << shmem_name << std::endl;
	if (mmap_ptr != MAP_FAILED) munmap (mmap_ptr, shmem_size);
	mmap_ptr = (byte *)MAP_FAILED;
	// close shmem_fd if its open
	if (shmem_fd > -1) close (shmem_fd);
	shmem_fd = -1;
	// unlink the shmem file
	shm_unlink (shmem_name.c_str());
	// close the semaphore
	sem_close (shmem_sem);
	shmem_sem = SEM_FAILED;

	WriteablePixelsFinish();
	remap_pix_plane (NULL, 0, 0);

	WriteableColorsFinish();
	remap_clr_plane (NULL, 0, 0);

}


void Transform::print_info() {

}

SharedImageMatrix* Transform::transform( SharedImageMatrix * matrix_IN ) {
	SharedImageMatrix *matrix_OUT = new SharedImageMatrix;
	matrix_OUT->cached_source = matrix_IN->shmem_name;
	matrix_OUT->operation = name;
	CacheStatus cache_status = matrix_OUT->fromCache ();
	
	if (cache_status == cache_write) {
		execute (matrix_IN, matrix_OUT);
		matrix_OUT->Cache();
		return (matrix_OUT);
	} else if (cache_status == cache_read) {
		return (matrix_OUT);
	} else {
		std::cerr << "Errors while recovering cache (cout):" << std::endl;
		exit (-1);
	}
}


EmptyTransform::EmptyTransform () {
	 Transform::name = "Empty";
};


void EmptyTransform::execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT ) {
	if( !( matrix_IN && matrix_OUT) )
		return;
	matrix_OUT->copy (*matrix_IN);
	std::cout << "Empty transform." << std::endl;
}

//===========================================================================

FourierTransform::FourierTransform () {
	 Transform::name = "Fourier";
}

/* fft 2 dimensional transform */
// http://www.fftw.org/doc/
void FourierTransform::execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT ) {
	if( !( matrix_IN && matrix_OUT) )
		return;
	
	matrix_OUT->copy (*matrix_IN);
	std::cout << "Performing transform " << name << std::endl;
	matrix_OUT->fft2();
}



//WNDCHARM_REGISTER_TRANSFORM(FourierTransform)

//===========================================================================

ChebyshevTransform::ChebyshevTransform () {
	 Transform::name = "Chebyshev";
}


void ChebyshevTransform::execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT ) {
	if( !( matrix_IN && matrix_OUT) )
		return;	
	matrix_OUT->copy (*matrix_IN);
	std::cout << "Performing transform " << name << std::endl;

	matrix_OUT->ChebyshevTransform(0);
}

//WNDCHARM_REGISTER_TRANSFORM(ChebyshevTransform)

//===========================================================================

WaveletTransform::WaveletTransform () {
	 Transform::name = "Wavelet";
};


void WaveletTransform::execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT ) {
	if( !( matrix_IN && matrix_OUT) )
		return;
	
	matrix_OUT->copy (*matrix_IN);
	std::cout << "Performing transform " << name << std::endl;
	matrix_OUT->Symlet5Transform();
}

//WNDCHARM_REGISTER_TRANSFORM(WaveletTransform)

//===========================================================================

EdgeTransform::EdgeTransform () {
	 Transform::name = "Edge";
}


void EdgeTransform::execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT ) {
	if( !( matrix_IN && matrix_OUT) )
		return;
	
	matrix_OUT->copy (*matrix_IN);
	std::cout << "Performing transform " << name << std::endl;
	matrix_OUT->EdgeTransform();
}

//WNDCHARM_REGISTER_TRANSFORM(EdgeTransform)

//===========================================================================

ColorTransform::ColorTransform () {
	 Transform::name = "Color";
}


void ColorTransform::execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT ) {
	if( !( matrix_IN && matrix_OUT) )
		return;
	
	matrix_OUT->copy (*matrix_IN);
	std::cout << "Performing transform " << name << std::endl;
	double temp_vec [COLORS_NUM+1];

	matrix_OUT->ColorTransform(temp_vec, 0);
	histogram_vals.assign(temp_vec, temp_vec+COLORS_NUM+1);
}

//WNDCHARM_REGISTER_TRANSFORM(ColorTransform)

//===========================================================================

HueTransform::HueTransform () {
	 Transform::name = "Hue";
}


void HueTransform::execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT ) {
	if( !( matrix_IN && matrix_OUT) )
		return;
	
	matrix_OUT->copy (*matrix_IN);
	std::cout << "Performing transform " << name << std::endl;
	matrix_OUT->ColorTransform(NULL,1);
}

//WNDCHARM_REGISTER_TRANSFORM(HueTransform)

