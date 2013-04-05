from pychrm.FeatureSet import *
import numpy as np
import ctypes

class PyImageMatrix (pychrm.ImageMatrix):
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
		stride_bytes = (self.width * PyImageMatrix.bytes_per_double, self.height * PyImageMatrix.bytes_per_double)
		ap = ctypes.cast(int(self.data_ptr()), ctypes.POINTER(ctypes.c_double * n_mem_doubles))
		return np.ndarray(shape=(self.height,self.width), dtype = np.double, strides=stride_bytes, buffer=ap.contents)
if __name__ == "__main__":
	foo = PyImageMatrix()
	foo.allocate (100,100)
	bar = foo.as_ndarray()
	print "The numpy returned by as_ndarray()"
	print dir (bar)
	foo = None
	print "last line"
