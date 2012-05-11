#!/usr/bin/env python
import stringformat
import pychrm

#================================================================
def main():

	feature_lists = LoadFeatureLists()
	full_list = "\n"
	full_list = full_list.join( feature_lists )
	original = pychrm.ImageMatrix()
	original.OpenImage("test.tif", 0, None, 0, 0)
	
	im_cache = {}
	im_cache[ '' ] = original

	global global_AllAlgorithms
	global global_AllTransforms

	global_AllAlgorithms = LoadFeatureAlgorithms()
	global_AllTransforms = LoadFeatureTransforms()

	feature_groups = []

	for feature_group_str in full_list.splitlines():
		feature_groups.append( ParseFeatureString( feature_group_str ) )

	print "{} feature groups".format( len( feature_groups ) )

	feature_names = []
	signatures = []
	for fg in feature_groups:
		print "Group {}".format(fg.Name)
		feature_vector = fg.CalculateFeatures( im_cache )
		count = 0
		for value in feature_vector:
			feature_names.append( fg.Name + " [{}]".format( count ) )
			signatures.append( value )	
			count += 1
	
	WriteFeaturesToSigFile( "pychrm_calculated.sig", feature_names, signatures )


#================================================================
class FeatureGroup:
	"""
	Attributes Name, Alg and Tforms are string representations, not links to the object
	"""

	Name = ""
	Alg = None
	Tforms = []
	def __init__( self, name_str = "", algorithm = None, tform_list = [] ):
		#print "Creating new FeatureGroup for string {}:".format(name_str)
		#print "\talgorithm: {}, transform list: {}".format( algorithm, tform_list )
		self.Name = name_str 
		self.Alg = algorithm
		self.Tforms = tform_list
	def CalculateFeatures( self, cached_pixel_planes ):
		"""Returns a tuple with the features"""
		pixel_plane = None
		try:
			#print "transforms: {}".format( self.Tforms )
			pixel_plane = RetrievePixelPlane( cached_pixel_planes, self.Tforms )
		except:
			raise
		return global_AllAlgorithms[ self.Alg ].calculate( pixel_plane )


#================================================================
def RetrievePixelPlane( image_matrix_cache, tform_list ):
	"""
	Returns the image matrix prescribed in tform_list
	If it already exists in cache, just return.
	If it doesn't exist calculates it
	Recurses through the compound transform chain in tform_list
	"""
	#print "passed in: {}".format( tform_list )
	requested_transform = " ".join( tform_list )
	#print "requesting pixel plane: '{}'".format( requested_transform )
	if requested_transform in image_matrix_cache:
		return image_matrix_cache[ requested_transform ]
	
	# Required transform isn't in the cache, gotta make it
	# Pop transforms off the end sequentially and check to see if
	# lower-level transforms have already been calculated and stored in cache

	# Can't begin if there isn't at least the raw (untransformed) pixel plane
	# already stored in the cache
	if image_matrix_cache is None or len(image_matrix_cache) == 0:
		raise ValueError( "Can't calculate features: couldn't find the original pixel plane" +\
		                  "to calculate features {}.".format( self.Name ) )

	sublist = tform_list[:]
	sublist.reverse()
	top_level_transform_name = sublist.pop()
	intermediate_pixel_plane = RetrievePixelPlane( image_matrix_cache, sublist )

	
	tformed_pp = global_AllTransforms[ top_level_transform_name ].transform( intermediate_pixel_plane )
	#assert( intermediate_pixel_plane ), "Pixel Plane returned from transform() was NULL"
	image_matrix_cache[ requested_transform ] = tformed_pp
	return tformed_pp


#================================================================
def ParseFeatureString( name ):
	"""Takes a string input, parses, and returns an instance of a FeatureGroup class"""
	#TBD: make a member function of the FeatureGroup
	# get the algorithm
	string_rep = name.rstrip( ")" )
	parsed = string_rep.split( ' (' )
	
	alg = parsed[0]
	if alg not in global_AllAlgorithms:
		raise KeyError( "Don't know about a feature algorithm with the name {}".format(alg) )
	
	tform_list = parsed[1:]
	try:
		tform_list.remove( "" )
	except ValueError:
		pass
	if len(tform_list) != 0:
		for tform in tform_list:
			if tform not in global_AllTransforms:
				raise KeyError( "Don't know about a transform named {}".format(tform) )
	return FeatureGroup( name, alg, tform_list )


#================================================================
def	WriteFeaturesToSigFile( file_name, feature_names, signatures ):
	"""Write a sig file."""
	with open( file_name, "w" ) as out_file:
		if len( feature_names ) != len( signatures ):
			raise ValueError( "Internal error: Can't write sig file: # of feature names " +\
			"({}) doesn't match # signatures ({})".format( len( feature_names ), len( signatures ) ) )
		for i in range( 0, len( feature_names ) ):
			out_file.write( "{val:0.6f} {name}\n".format( val=signatures[i], name=feature_names[i] ) )

#================================================================
def LoadFeatureAlgorithms():
	out = {}
	out[ 'Chebyshev Coefficients' ] = pychrm.ChebyshevCoefficients()
	out[ 'Chebyshev-Fourier Coefficients' ] = pychrm.ChebyshevFourierCoefficients()
	out[ 'Zernike Coefficients' ] = pychrm.ZernikeCoefficients()
	out[ 'Haralick Textures' ] = pychrm.HaralickTextures()
	out[ 'Multiscale Histograms' ] = pychrm.MultiscaleHistograms()
	out[ 'Tamura Textures' ] = pychrm.TamuraTextures()
	out[ 'Comb Moments' ] = pychrm.CombFirstFourMoments()
	out[ 'Radon Coefficients' ] = pychrm.RadonCoefficients()
	out[ 'Fractal Features' ] = pychrm.FractalFeatures()
	out[ 'Pixel Intensity Statistics' ] = pychrm.PixelIntensityStatistics()
	out[ 'Edge Features' ] = pychrm.EdgeFeatures()
	out[ 'Object Features' ] = pychrm.ObjectFeatures()
	out[ 'Gabor Textures' ] = pychrm.GaborTextures()
	out[ 'Gini Coefficient' ] = pychrm.GiniCoefficient()
	return out


#================================================================
def LoadFeatureTransforms():
	out = {}
	out[ 'Fourier' ] = pychrm.FourierTransform()
	out[ 'Chebyshev' ] = pychrm.ChebyshevTransform()
	out[ 'Wavelet' ] = pychrm.WaveletTransform()
	out[ 'Edge' ] = pychrm.EdgeTransform()
	out[ 'Color' ] = pychrm.ColorTransform()
	out[ 'Hue' ] = pychrm.HueTransform()
	return out


#================================================================
def LoadFeatureLists():

	small = """Chebyshev Coefficients ()
Chebyshev Coefficients (Fourier ())
Chebyshev-Fourier Coefficients ()
Chebyshev-Fourier Coefficients (Fourier ())
Comb Moments ()
Comb Moments (Chebyshev ())
Comb Moments (Chebyshev (Fourier ()))
Comb Moments (Fourier ())
Comb Moments (Wavelet ())
Comb Moments (Wavelet (Fourier ()))
Edge Features ()
Gabor Textures ()
Haralick Textures ()
Haralick Textures (Chebyshev ())
Haralick Textures (Chebyshev (Fourier ()))
Haralick Textures (Fourier ())
Haralick Textures (Wavelet ())
Haralick Textures (Wavelet (Fourier ()))
Multiscale Histograms ()
Multiscale Histograms (Chebyshev ())
Multiscale Histograms (Chebyshev (Fourier ()))
Multiscale Histograms (Fourier ())
Multiscale Histograms (Wavelet ())
Multiscale Histograms (Wavelet (Fourier ()))
Object Features ()
Radon Coefficients ()
Radon Coefficients (Chebyshev ())
Radon Coefficients (Chebyshev (Fourier ()))
Radon Coefficients (Fourier ())
Tamura Textures ()
Tamura Textures (Chebyshev ())
Tamura Textures (Chebyshev (Fourier ()))
Tamura Textures (Fourier ())
Tamura Textures (Wavelet ())
Tamura Textures (Wavelet (Fourier ()))
Zernike Coefficients ()
Zernike Coefficients (Fourier ())"""

	large = """Chebyshev Coefficients (Chebyshev ())
Chebyshev Coefficients (Edge ())
Chebyshev Coefficients (Fourier (Edge ()))
Chebyshev Coefficients (Fourier (Wavelet ()))
Chebyshev Coefficients (Wavelet ())
Chebyshev Coefficients (Wavelet (Edge ()))
Chebyshev-Fourier Coefficients (Chebyshev ())
Chebyshev-Fourier Coefficients (Edge ())
Chebyshev-Fourier Coefficients (Fourier (Edge ()))
Chebyshev-Fourier Coefficients (Fourier (Wavelet ()))
Chebyshev-Fourier Coefficients (Wavelet ())
Chebyshev-Fourier Coefficients (Wavelet (Edge ()))
Comb Moments (Chebyshev (Wavelet ()))
Comb Moments (Edge ())
Comb Moments (Fourier (Chebyshev ()))
Comb Moments (Fourier (Edge ()))
Comb Moments (Fourier (Wavelet ()))
Comb Moments (Wavelet (Edge ()))
Fractal Features ()
Fractal Features (Chebyshev ())
Fractal Features (Chebyshev (Fourier ()))
Fractal Features (Chebyshev (Wavelet ()))
Fractal Features (Edge ())
Fractal Features (Fourier ())
Fractal Features (Fourier (Chebyshev ()))
Fractal Features (Fourier (Edge ()))
Fractal Features (Fourier (Wavelet ()))
Fractal Features (Wavelet ())
Fractal Features (Wavelet (Edge ()))
Fractal Features (Wavelet (Fourier ()))
Gini Coefficient ()
Gini Coefficient (Chebyshev ())
Gini Coefficient (Chebyshev (Fourier ()))
Gini Coefficient (Chebyshev (Wavelet ()))
Gini Coefficient (Edge ())
Gini Coefficient (Fourier ())
Gini Coefficient (Fourier (Chebyshev ()))
Gini Coefficient (Fourier (Edge ()))
Gini Coefficient (Fourier (Wavelet ()))
Gini Coefficient (Wavelet ())
Gini Coefficient (Wavelet (Edge ()))
Gini Coefficient (Wavelet (Fourier ()))
Haralick Textures (Chebyshev (Wavelet ()))
Haralick Textures (Edge ())
Haralick Textures (Fourier (Chebyshev ()))
Haralick Textures (Fourier (Edge ()))
Haralick Textures (Fourier (Wavelet ()))
Haralick Textures (Wavelet (Edge ()))
Multiscale Histograms (Chebyshev (Wavelet ()))
Multiscale Histograms (Edge ())
Multiscale Histograms (Fourier (Chebyshev ()))
Multiscale Histograms (Fourier (Edge ()))
Multiscale Histograms (Fourier (Wavelet ()))
Multiscale Histograms (Wavelet (Edge ()))
Pixel Intensity Statistics ()
Pixel Intensity Statistics (Chebyshev ())
Pixel Intensity Statistics (Chebyshev (Fourier ()))
Pixel Intensity Statistics (Chebyshev (Wavelet ()))
Pixel Intensity Statistics (Edge ())
Pixel Intensity Statistics (Fourier ())
Pixel Intensity Statistics (Fourier (Chebyshev ()))
Pixel Intensity Statistics (Fourier (Edge ()))
Pixel Intensity Statistics (Fourier (Wavelet ()))
Pixel Intensity Statistics (Wavelet ())
Pixel Intensity Statistics (Wavelet (Edge ()))
Pixel Intensity Statistics (Wavelet (Fourier ()))
Radon Coefficients (Chebyshev (Wavelet ()))
Radon Coefficients (Edge ())
Radon Coefficients (Fourier (Chebyshev ()))
Radon Coefficients (Fourier (Edge ()))
Radon Coefficients (Fourier (Wavelet ()))
Radon Coefficients (Wavelet ())
Radon Coefficients (Wavelet (Edge ()))
Radon Coefficients (Wavelet (Fourier ()))
Tamura Textures (Chebyshev (Wavelet ()))
Tamura Textures (Edge ())
Tamura Textures (Fourier (Chebyshev ()))
Tamura Textures (Fourier (Edge ()))
Tamura Textures (Fourier (Wavelet ()))
Tamura Textures (Wavelet (Edge ()))
Zernike Coefficients (Chebyshev ())
Zernike Coefficients (Edge ())
Zernike Coefficients (Fourier (Edge ()))
Zernike Coefficients (Fourier (Wavelet ()))
Zernike Coefficients (Wavelet ())
Zernike Coefficients (Wavelet (Edge ()))"""

	return small, large


#================================================================
if __name__=="__main__":
	main()


