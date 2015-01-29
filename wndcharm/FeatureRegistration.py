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

try:
	from . import wndcharm
except:
	import wndcharm
#================================================================
def LoadFeatureAlgorithms():
	out = {}
	out[ 'Chebyshev Coefficients' ] = wndcharm.ChebyshevCoefficients()
	out[ 'Chebyshev-Fourier Coefficients' ] = wndcharm.ChebyshevFourierCoefficients()
	out[ 'Zernike Coefficients' ] = wndcharm.ZernikeCoefficients()
	out[ 'Haralick Textures' ] = wndcharm.HaralickTextures()
	out[ 'Multiscale Histograms' ] = wndcharm.MultiscaleHistograms()
	out[ 'Tamura Textures' ] = wndcharm.TamuraTextures()
	out[ 'Comb Moments' ] = wndcharm.CombFirstFourMoments()
	out[ 'Radon Coefficients' ] = wndcharm.RadonCoefficients()
	out[ 'Fractal Features' ] = wndcharm.FractalFeatures()
	out[ 'Pixel Intensity Statistics' ] = wndcharm.PixelIntensityStatistics()
	out[ 'Edge Features' ] = wndcharm.EdgeFeatures()
	out[ 'Object Features' ] = wndcharm.ObjectFeatures()
	out[ 'Gabor Textures' ] = wndcharm.GaborTextures()
	out[ 'Gini Coefficient' ] = wndcharm.GiniCoefficient()
	return out


#================================================================
def LoadFeatureTransforms():
	out = {}
	out[ 'Fourier' ] = wndcharm.FourierTransform()
	out[ 'Chebyshev' ] = wndcharm.ChebyshevTransform()
	out[ 'Wavelet' ] = wndcharm.WaveletTransform()
	out[ 'Edge' ] = wndcharm.EdgeTransform()
	out[ 'Color' ] = wndcharm.ColorTransform()
	out[ 'Hue' ] = wndcharm.HueTransform()
	return out


#================================================================
def LoadSmallAndLargeFeatureSetStringLists():

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
#Gini Coefficient ()
#Gini Coefficient (Chebyshev ())
#Gini Coefficient (Chebyshev (Fourier ()))
#Gini Coefficient (Chebyshev (Wavelet ()))
#Gini Coefficient (Edge ())
#Gini Coefficient (Fourier ())
#Gini Coefficient (Fourier (Chebyshev ()))
#Gini Coefficient (Fourier (Edge ()))
#Gini Coefficient (Fourier (Wavelet ()))
#Gini Coefficient (Wavelet ())
#Gini Coefficient (Wavelet (Edge ()))
#Gini Coefficient (Wavelet (Fourier ()))
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

