/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/* Copyright (C) 2007 Open Microscopy Environment                                */
/*       Massachusetts Institue of Technology,                                   */
/*       National Institutes of Health,                                          */
/*       University of Dundee                                                    */
/*                                                                               */
/*                                                                               */
/*                                                                               */
/*    This library is free software; you can redistribute it and/or              */
/*    modify it under the terms of the GNU Lesser General Public                 */
/*    License as published by the Free Software Foundation; either               */
/*    version 2.1 of the License, or (at your option) any later version.         */
/*                                                                               */
/*    This library is distributed in the hope that it will be useful,            */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          */
/*    Lesser General Public License for more details.                            */
/*                                                                               */
/*    You should have received a copy of the GNU Lesser General Public           */
/*    License along with this library; if not, write to the Free Software        */
/*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA  */
/*                                                                               */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Written by:                                                                   */
/*      Ilya G. Goldberg <goldbergil [at] mail [dot] nih [dot] gov>              */
/*      Christopher E. Coletta <colettace [at] mail [dot] nih [dot] gov>         */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifndef __TRANSFORMS_H_
#define __TRANSFORMS_H_

#include <string>
#include <vector>
#include "Tasks.h"



/*! ImageTransform
 *  defines the interface for all inheriting transform classes
 *  Turns any class that inherits this interface into a singleton
 */
class ImageMatrix; // forward declaration
class ImageTransform : public ComputationTask {
	public:
		virtual void execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const = 0;
		virtual bool register_task() const;
	protected:
		ImageTransform (const std::string &s) : ComputationTask (s, ImageTransformTask) {};
};

class EmptyTransform : public ImageTransform {
	public:
		EmptyTransform () : ImageTransform ("Empty") {};
		EmptyTransform (const std::string &s) : ImageTransform (s) {};
		EmptyTransform (const char *s) : ImageTransform (s) {};
		virtual void execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const;
};

class FourierTransform : public ImageTransform {
	public:
		FourierTransform();
		virtual void execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const;
};

class ChebyshevTransform: public ImageTransform {
	public:
		ChebyshevTransform();
		virtual void execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const;
};

class WaveletTransform : public ImageTransform {
	public:
		WaveletTransform();
		virtual void execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const;
};

class EdgeTransform : public ImageTransform {
	public:
		EdgeTransform();
		virtual void execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const;
};

class ColorTransform : public ImageTransform {
	public:
		ColorTransform();
		virtual void execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const;
};

class HueTransform : public ImageTransform {
	public:
		HueTransform();
		virtual void execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const;
};

#endif

