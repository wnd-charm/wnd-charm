#ifndef __SHARED_IMAGE_MATRIX_H_
#define __SHARED_IMAGE_MATRIX_H_

#include <string>
#include <sys/mman.h>  // mmap, shm_open, MAP_FAILED
#include "cmatrix.h"
#include "WORMfile.h"
/*! SharedImageMatrix
* inherits from ImageMatrix to store data in a named mmap for memory sharing b/w processes
*/
enum CacheStatus {csUNKNOWN, csREAD, csWRITE, csWAIT, csERROR};
class SharedImageMatrix: public ImageMatrix {
	public:
		SharedImageMatrix () : ImageMatrix () {
			cached_source = "";
			operation = "";
			shmem_name = "";
			was_cached = false;
			shmem_size = 0;
			shmem_fd = -1;
			mmap_ptr = (byte *)MAP_FAILED;
			error_str = "";
			cache_status = csUNKNOWN;
		}
		// N.B.: The source and the operation can't both be empty.
		void fromCache (const std::string cached_source_in, const std::string operation_in = "");
		void Cache ();

	// Accessors for read-only fields
		const std::string &Error () const {return (error_str);};
		const CacheStatus &Status() const {return (cache_status);};
		const std::string &GetShmemName() const {return shmem_name; };

	// Class methods
		static void DisableCacheRead (const bool status) {never_read = status;};
		static void DisableDestructorCacheCleanup (const bool status) {disable_destructor_cache_cleanup = status;};

	// Overrides of parent class methods
		virtual void allocate (unsigned int w, unsigned int h) ;
		virtual int OpenImage(char *image_file_name,            // load an image of any supported format
			int downsample, rect *bounding_rect,
			double mean, double stddev);
		virtual ImageMatrix &transform (const ImageTransform *transform) const;

		virtual ~SharedImageMatrix();                                 // destructor

	private:
		// class statics
		static size_t shmem_page_size;
		static bool never_read;
		static bool disable_destructor_cache_cleanup;
		static pid_t PID; // statics are per-process, even when linking against a shared library with statics
		static std::string PID_string; // statics are per-process, even when linking against a shared library with statics

		// private class methods
		static const size_t calc_shmem_size (const unsigned int w, const unsigned int h, const enum ColorModes ColorMode, size_t &clr_plane_offset, size_t &shmem_data_offset);
		// private instance methods
		void SetShmemName();

		// private object fields
		WORMfile lock_file;
		std::string cached_source;        // the shmem_name of the source.
		std::string operation;            // the operation on the cached_source
		std::string shmem_name;           // concatenated cached_source + operation, MD5 digest, Base 64-encoded
		bool was_cached;
		size_t shmem_size;
		int shmem_fd;
		byte *mmap_ptr;
		std::string error_str;           // String reporting errors with the cache.
		CacheStatus cache_status;

		// shared-memory layout:
		// last sizeof(shmem_data) bytes are shmem_data.
		// First pages are the array of doubles for pix_plane
		// Second set of pages are array of HSVColor for clr_plane if shmem_data->ColorMode != cmGRAY
		// The pix_plane matrix storage ends on a page boundary so that clr_plane can begin at a page boundary
		// The shmem_data storage does not necessarily begin at a page boundary, but it is at the end of the last page.
};

#endif // __SHARED_IMAGE_MATRIX_H_
