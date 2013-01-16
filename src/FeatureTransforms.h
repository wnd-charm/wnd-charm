#ifndef __TRANSFORMS_H_
#define __TRANSFORMS_H_

#include <string>
#include <sys/mman.h>  // mmap, shm_open
#include <semaphore.h> // semaphores
#include "wndchrm_error.h"
#include "cmatrix.h"


/*! SharedImageMatrix
* inherits from ImageMatrix to store data in a named mmap for memory sharing b/w processes
*/
enum CacheStatus { cache_unknown, cache_wait, cache_read, cache_write };

class SharedImageMatrix: public ImageMatrix {
	public:
		SharedImageMatrix () : ImageMatrix () {
			shmem_name = "";
			shmem_size = 0;
			shmem_fd = -1;
			mmap_ptr = (byte *)MAP_FAILED;
			shmem_sem = SEM_FAILED;
		}
		virtual ~SharedImageMatrix();                                 // destructor
		void SetShmemName( const std::string &final_op );
		const size_t calc_shmem_size (const unsigned int w, const unsigned int h, size_t &clr_plane_offset, size_t &shmem_data_offset);
		void allocate (unsigned int w, unsigned int h) ;
		const CacheStatus fromCache ( const std::string &final_op );
		void Cache ();
		virtual int OpenImage(char *image_file_name,            // load an image of any supported format
			int downsample, rect *bounding_rect,
			double mean, double stddev);

		static size_t shmem_page_size;

		std::string shmem_name;
		size_t shmem_size;
		int shmem_fd;
		byte *mmap_ptr;
		sem_t *shmem_sem;

		// shared-memory layout:
		// last sizeof(shmem_data) bytes are shmem_data.
		// First pages are the array of doubles for pix_data
		// Second set of pages are array of HSVColor for clr_data if shmem_data->ColorMode != cmGRAY
};

/*! Transform
 *  defines the interface for all inheriting transform classes
 *  Turns any class that inherits this interface into a singleton
 */
class Transform {
	public:
		SharedImageMatrix* transform( SharedImageMatrix * matrix_IN );
		virtual void execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT ) = 0;
		std::string name;
		void print_info();
	protected:
		Transform (const std::string &s) { name = s;}
		Transform (const char *s) { name = s;}
		Transform() {};
};

class EmptyTransform : public Transform {
	public:
		EmptyTransform ();
		virtual void execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT );
};


class FourierTransform : public Transform {
	public:
		FourierTransform();
		virtual void execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT );
};

class ChebyshevTransform: public Transform {
	public:
		ChebyshevTransform();
		virtual void execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT );
};

class WaveletTransform : public Transform {
	public:
		WaveletTransform();
		virtual void execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT );
};

class EdgeTransform : public Transform {
	public:
		EdgeTransform();
		virtual void execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT );
};

class ColorTransform : public Transform {
	public:
		ColorTransform();
		vector<double> histogram_vals;
		virtual void execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT );
};

class HueTransform : public Transform {
	public:
		HueTransform();
		virtual void execute( const SharedImageMatrix * matrix_IN, SharedImageMatrix * matrix_OUT );
};

/*	
#define WNDCHARM_REGISTER_TRANSFORM(tform_name) \
struct tform_name##TransformRegistrar \
{ \
  tform_name##TransformRegistrar() \
  { \
    FeatureNames *phonebook = FeatureNames::get_instance(); \
		tform_name *tform_instance = new tform_name; \
		int retval = phonebook->register_transform( tform_instance->name, dynamic_cast<Transform*>( tform_instance ) ); \
  } \
}; \
static tform_name##TransformRegistrar tform_name##TransformRegistrar_instance;
//std::cout << "call to register_transform " << #tform_name << " returned " << retval << std::endl; \
*/
#endif

