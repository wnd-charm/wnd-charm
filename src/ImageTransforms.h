#ifndef __TRANSFORMS_H_
#define __TRANSFORMS_H_

#include <string>
#include <vector>
#include "cmatrix.h"



/*! ImageTransform
 *  defines the interface for all inheriting transform classes
 *  Turns any class that inherits this interface into a singleton
 */
class ImageMatrix; // forward declaration
class ImageTransform {
	public:
		std::string name;
		virtual void execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const = 0;
		void print_info();
	protected:
		ImageTransform (const std::string &s) { name = s;}
		ImageTransform (const char *s) { name = s;}
		ImageTransform() {};
};

class EmptyTransform : public ImageTransform {
	public:
		EmptyTransform () : ImageTransform ("Empty") {};
		EmptyTransform (const std::string &s) : ImageTransform (s) {};
		EmptyTransform (const char *s) : ImageTransform (s) {};
		virtual void execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const;
};

class FourierTransform : public ImageTransform {
	public:
		FourierTransform();
		virtual void execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const;
};

class ChebyshevTransform: public ImageTransform {
	public:
		ChebyshevTransform();
		virtual void execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const;
};

class WaveletTransform : public ImageTransform {
	public:
		WaveletTransform();
		virtual void execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const;
};

class EdgeTransform : public ImageTransform {
	public:
		EdgeTransform();
		virtual void execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const;
};

class ColorTransform : public ImageTransform {
	public:
		ColorTransform();
		virtual void execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const;
};

class HueTransform : public ImageTransform {
	public:
		HueTransform();
		virtual void execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const;
};

#endif

