"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                               
 Copyright (C) 2015 National Institutes of Health 

    This library is free software; you can redistribute it and/or              
    modify it under the terms of the GNU Lesser General Public                 
    License as published by the Free Software Foundation; either               
    version 2.1 of the License, or (at your option) any later version.         
                                                                               
    This library is distributed in the hope that it will be useful,            
    but WITHOUT ANY WARRANTY; without even the implied warranty of             
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          
    Lesser General Public License for more details.                            
                                                                               
    You should have received a copy of the GNU Lesser General Public           
    License along with this library; if not, write to the Free Software        
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA  
                                                                               
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                               
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Written by:  Christopher Coletta (github.com/colettace)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""

import wndcharm
from wndcharm import ImageMatrix
import numpy as np
import ctypes

class PyImageMatrix (ImageMatrix):
	bytes_per_double = np.dtype(np.double).itemsize
	def __init__(self):
		super(PyImageMatrix, self).__init__()
	def as_ndarray(self):
		# self.data_ptr() is None unless allocate() has been called
		if (self.data_ptr() is None):
			return None
		# works, No control over strides
		# return np.ctypeslib.as_array(ctypes.cast(int(self.data_ptr()), ctypes.POINTER(ctypes.c_double)),shape=(self.height,self.width))
		#
		# works, full control over ndarray constructor
		n_mem_doubles = self.width * self.height
		ap = ctypes.cast(int(self.data_ptr()), ctypes.POINTER(ctypes.c_double * n_mem_doubles))
# strides aren't quite right yet...
# 		stride_bytes = (self.width * PyImageMatrix.bytes_per_double, self.height * PyImageMatrix.bytes_per_double)
# 		return np.ndarray(shape=(self.height,self.width), dtype = np.double, strides=stride_bytes, buffer=ap.contents)
		return np.ndarray(shape=(self.height,self.width), dtype = np.double, buffer=ap.contents)
if __name__ == "__main__":
	im = PyImageMatrix()
	im.allocate (200,200)
	mat = im.as_ndarray()
	print "mat = im.as_ndarray()"
	print mat
	print "mat.ctypes.data: {0:x}".format(mat.ctypes.data)
	mat[:] = 0.
	print "mat.ctypes.data: {0:x} after mat[:] = 0".format(mat.ctypes.data)
	print mat
	print "deallocating mat"
	mat = None
	print "deallocating im"
	im = None
	print "last line"
