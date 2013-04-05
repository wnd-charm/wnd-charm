import pychrm
from pychrm import ImageMatrix
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
