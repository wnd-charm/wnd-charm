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
import numpy as np
from . import feature_vector_major_version
from . import feature_vector_minor_version_from_num_features
from .utils import normalize_by_columns

class WrongFeatureSetVersionError( Exception ):
    pass

class IncompleteFeatureSetError( Exception ):
    pass

# Couldn't get this "Python singleton inherited from swig-wrapped C++ object" to work:
#*** NotImplementedError: Wrong number or type of arguments for overloaded function 'FeatureComputationPlan_add'.
#  Possible C/C++ prototypes are:
#    FeatureComputationPlan::add(std::string const &)
#    FeatureComputationPlan::add(FeatureGroup const *)
# "self" below was of type "<wndcharm.FeatureVector.PyFeatureComputationPlan;  >"
# when what was required was a SWIG proxy object to translate native python strings
# into std::string
# "<wndcharm.wndcharm.FeatureComputationPlan; proxy of <Swig Object of type 'FeatureComputationPlan *' at 0x111263cf0> >"
##================================================================
#class PyFeatureComputationPlan( wndcharm.FeatureComputationPlan ):
#    """Contains a cache to save memory, as there may be tens of thousands of samples
#    and therefore the same number of of redundant instances of the same computation plan."""
#
#    plan_cache = {}
#
#    def __new__( cls, feature_list, name='custom' ):
#        """Takes list of feature strings and chops off bin number at the first
#        space on right, e.g., "feature alg (transform()) [bin]" """
#
#        feature_groups = frozenset( [ feat.rsplit(" ",1)[0] for feat in feature_list ] )
#
#        if feature_groups in cls.plan_cache:
#            return cls.plan_cache[ feature_groups ]
#
#        self = super( PyFeatureComputationPlan, cls ).__new__( cls, name )
#        [ self.add( family ) for family in feature_groups ]
#
#        cls.plan_cache[ feature_groups ] = self
#        return self

# instead implement with a global dict to serve as feature plan cache

plan_cache = {}

def GenerateFeatureComputationPlan( feature_list, name='custom' ):
    """Takes list of feature strings and chops off bin number at the first
    space on right, e.g., "feature alg (transform()) [bin]" """

    global plan_cache
    feature_groups = frozenset( [ feat.rsplit(" ",1)[0] for feat in feature_list ] )

    if feature_groups in plan_cache:
        return plan_cache[ feature_groups ]

    obj = wndcharm.FeatureComputationPlan( name )
    [ obj.add( family ) for family in feature_groups ]

    plan_cache[ feature_groups ] = obj
    return obj


#############################################################################
# class definition of FeatureVector
#############################################################################
class FeatureVector( object ):
    """
    FeatureVector contains features for a single image or ROI, as well as all the sampling
    options. It is the values contained here that are grouped together to form FeatureSpaces."""

    import re
    sig_filename_parser = re.compile( ''.join( [
    # basename:
    r'(?P<basename>.+?)',
    # ROI:
    r'(?P<roi>-B(?P<x>\d+)_(?P<y>\d+)_(?P<w>\d+)_(?P<h>\d+))?',
    # downsampling:
    r'(?:-d(?P<downsample>\d+))?',
    # pixel intensity adjustment:
    r'(?:-S(?P<pixel_intensity_mean>\d+)(?:_(?P<pixel_intensity_stddev>\d+))?)?',
    # rotations:
    r'(?:-R_(?P<rot>\d))?',
    # tiling info:
    r'(?P<tiling_scheme>-t(?P<tile_num_rows>\d+)(?:x(?P<tile_num_cols>\d+))?_(?P<tile_row_index>\d+)_(?P<tile_col_index>\d+))?',
    # color features:
    r'(?P<color>-c)?',
    # long feature set:
    r'(?P<long>-l)?',
    # extension: .sig or .pysig
    r'\.(?:py)?sig$' ] ) )

    #==============================================================
    def __init__( self, **kwargs ):

        #: Row name, common across all channels
        self.name = None

        #: feature names
        self.feature_names = None
        #: feature vector
        self.values = None
        #: feature_maxima, feature_minima, and normalized_against are members
        #: this object shares with the FeatureSpaceObject, as well as implementations
        #: of FeatureReduce() and Normalize(), and CompatibleFeatureSetVersion()
        self.feature_maxima = None
        self.feature_minima = None
        self.normalized_against = None
        #: the prefix string to which sampling options will be appended to form .sig filepath
        self.basename = None
        #: Can also be a reference to a wndchrm.ImageMatrix object
        self.source_filepath = None
        #: Path to .sig file, in future hdf/sql file
        self.auxiliary_feature_storage = None
        #: label is stringified ground truth
        self.label = None
        self.ground_truth = None
        self.samplegroupid = None
        self.channel = None
        self.time_index = None
        self.tiling_scheme = None
        #: If no ROI image subsample, whole image is tile 1 of 1 in this sample group.
        self.tile_num_rows = 1
        self.tile_num_cols = 1
        #: indices count from 0
        self.tile_row_index = 0
        self.tile_col_index = 0
        self.samplesequenceid = 0
        #: downsample (in percents)
        self.downsample = 0
        self.pixel_intensity_mean = None
        self.pixel_intensity_stddev = None
        self.roi = None
        self.h = None
        self.w = None
        self.x = None
        self.y = None
        self.z = None
        self.z_delta = None
        self.rot = None
        self.fs_col = 0

        #: self.num_features should always be == len( self.feature_names ) == len( self.values )
        self.num_features = None

        # WND-CHARM feature bank-specific params:
        self.color = None
        self.long = None
        self.feature_set_version = None
        self.feature_computation_plan = None

        self.Update( **kwargs )
    #==============================================================
    def __len__( self ):
        try:
            length = len( self.feature_names )
        except:
            length = 0
        return length

    #==============================================================
    def __str__( self ):
        outstr = '<' + self.__class__.__name__
        if self.name is not None:
            if len(self.name) > 30:
                name = '...' + self.name[-30:]
            else:
                name = self.name
            outstr += ' "' + name + '"'
        if self.label is not None:
            outstr += ' label="' + self.label + '"'
        if self.feature_names is not None:
            outstr += ' n_features=' + str( len( self ) )
        if self.samplegroupid is not None:
            outstr += ' grp=' + str( self.samplegroupid )
        if self.samplesequenceid is not None:
            outstr += ' seq=' + str( self.samplesequenceid )
        if self.fs_col is not None:
            outstr += ' fs_col=' + str( self.fs_col )
        return outstr + '>'

    #==============================================================
    def __repr__( self ):
        return str(self)

    #==============================================================
    def Update( self, force=False, **kwargs ):
        """force - if True then kwargs coming in with None values will overwrite
        members in self, even if those members are non-None before calling Update."""

        self_namespace = vars( self )
        for key, val in kwargs.iteritems():
            # don't bother setting a val if it's None unless force
            if val is not None or force:  
                if key in self_namespace:
                    #if self_namespace[ key ] != None and self_namespace[ key ] != val:
                    #    from warnings import warn
                    #    warn( "Overwriting attrib {0} old val {1} new val {2}".format( key, self_namespace[ key ], val ) )
                    self_namespace[ key ] = val
                else:
                    raise AttributeError( 'No instance variable named "{0}" in class {1}'.format(
                        key, self.__class__.__name__ ) )

        def ReturnNumFeaturesBasedOnMinorFeatureVectorVersion( minor_fvv ):
            if major == 1:
                num_feats_dict = feature_vector_minor_version_from_num_features_v1
            else:
                num_feats_dict = feature_vector_minor_version_from_num_features

            for num_feats, version in num_feats_dict.iteritems():
                if version == minor_fvv:
                    return num_feats
            return None

        # FIXME: feature_set_version refers to WND-CHARM specific feature set specifications.
        # Want to be able to handle other feature sets from other labs in the future.
        if self.feature_set_version == None:
            major = feature_vector_major_version

            # The feature_set_version helps describe what features contained in the set.
            # Major version has to do with fixing bugs in the WND_CHARM algorithm code base.
            # The minor version describes the composition of features in the feature set.
            # Minor versions 1-4 have specific combination of WND-CHARM features.
            # Minor version 0 refers to user-defined combination of features.

            # Check to see if there is a user-defined set of features for this feature vector:
            if self.feature_computation_plan or self.feature_names:
                # set num_features
                if self.feature_names:
                    self.num_features = len( self.feature_names )
                else:
                    self.num_features = self.feature_computation_plan.n_features

                if self.num_features not in feature_vector_minor_version_from_num_features:
                    minor = 0
                else:
                    # FIXME: If features are out of order, should have a minor version of 0
                    minor = feature_vector_minor_version_from_num_features[ len( self.feature_names ) ]
            else:
                if not self.long:
                    if not self.color:
                        minor = 1
                    else:
                        minor = 3
                else:
                    if not self.color:
                        minor = 2
                    else:
                        minor = 4
                self.num_features = ReturnNumFeaturesBasedOnMinorFeatureVectorVersion( minor )
            self.feature_set_version = '{0}.{1}'.format( major, minor )
        else:
            major, minor = [ int( val ) for val in self.feature_set_version.split('.') ]
            self.num_features = ReturnNumFeaturesBasedOnMinorFeatureVectorVersion( minor )

        # When reading in sampling opts from the path, they get pulled out as strings
        # instead of ints:
        if self.tile_row_index is not None and type( self.tile_row_index ) != int:
            self.tile_row_index = int( self.tile_row_index )
        if self.tile_col_index is not None and type( self.tile_col_index ) != int:
            self.tile_col_index = int( self.tile_col_index )
        if self.tile_num_rows is not None and type( self.tile_num_rows ) != int:
            self.tile_num_rows = int( self.tile_num_rows )
        if self.tile_num_cols is not None and type( self.tile_num_cols ) != int:
            self.tile_num_cols = int( self.tile_num_cols )
        if self.samplegroupid is not None and type( self.samplegroupid ) != int:
            self.samplegroupid = int( self.tile_num_cols )

        # sequence order:
        # index 0 = position row 0, col 0
        # index 1 = position row 0, col 1
        # index 2 = position row 1, col 0, etc...
        self.samplesequenceid = (self.tile_row_index * self.tile_num_cols) + self.tile_col_index
        return self

    #==============================================================
    def Derive( self, **kwargs ):
        """Make a copy of this FeatureVector, except members passed as kwargs"""

        from copy import deepcopy
        new_obj = self.__class__()
        self_namespace = vars( self )
        new_obj_namespace = vars( new_obj )

        # skip these (if any):
        convenience_view_members = []

        # Are all keys in kwargs valid instance attribute names?
        invalid_kwargs = set( kwargs.keys() ) - set( self_namespace.keys() )
        if len( invalid_kwargs ) > 0:
            raise ValueError( "Invalid keyword arg(s) to Derive: {0}".format( invalid_kwargs ) )

        # Go through all of self's members and copy them to new_fs
        # unless a key-val pair was passed in as kwargs
        for key in self_namespace:
            if key in convenience_view_members:
                continue
            if key in kwargs:
                new_obj_namespace[key] = kwargs[key]
            else:
                new_obj_namespace[key] = deepcopy( self_namespace[key] )
        return new_obj

    #==============================================================
    def __deepcopy__( self, memo ):
        """Make a deepcopy of this FeatureVector"""
        return self.Derive()

    #==============================================================
    def GenerateSigFilepath( self ):
        """The C implementation of wndchrm placed feature metadata
        in the filename in a specific order, recreated here."""

        from os.path import splitext

        # FIXME: sigpaths for FeatureVectors with different channels
        # may have sig file names that will collide/overwrite each other.
        if self.basename:
            base = self.basename
        elif isinstance( self.source_filepath, wndcharm.ImageMatrix ) and \
                self.source_filepath.source:
            base, ext = splitext( self.source_filepath.source )
            self.basename = base
        elif self.source_filepath:
            base, ext = splitext( self.source_filepath )
            self.basename = base
        elif self.name:
            # ext may be nothing, that's ok
            base, ext = splitext( self.name )
            self.basename = base
        else:
            raise ValueError( 'Need for "basename" or "source_filepath" or "name" attribute in FeatureVector object to be set to generate sig filepath.')

        self_namespace = vars(self)
        
        if self.x is not None and self.y is not None and self.w is not None and self.h is not None:
            base += "-B{0}_{1}_{2}_{3}".format( self.x, self.y, self.w, self.h )
        if self.downsample:
            base += "-d" + str(self.downsample)
        if self.pixel_intensity_mean is not None:
            base += "-S" + str(self.pixel_intensity_mean)
            if self.pixel_intensity_stddev is not None:
                base += "-S" + str(self.pixel_intensity_stddev)
        if self.rot is not None:
            base += "-R_" + str(self.rot)
        if self.tile_num_rows and self.tile_num_rows != 1:
            base += "-t" + str(self.tile_num_rows)
            if self.tile_num_cols and self.tile_num_cols != 1:
                base +="x" + str(self.tile_num_cols)
            if self.tile_row_index is not None and self.tile_col_index is not None:
                base += "_{0}_{1}".format( self.tile_row_index, self.tile_col_index )
            else:
                raise ValueError('Need to specify tile_row_index and tile_col_index in self for tiling params')
        if self.color:
            base += '-c'
        if self.long:
            base += '-l'

        return base + '.sig'

    #================================================================
    def GenerateFeatures( self, write_to_disk=True, quiet=True ):
        """@brief Loads precalculated features, or calculates new ones, based on which instance
        attributes have been set, and what their values are.

        write_to_disk (bool) - save features to text file which by convention has extension ".sig"
        
        Returns self for convenience."""

        # 0: What features does the user want?
        # 1: are there features already calculated somewhere?
        # 2: if so are they current/complete/expected/correct?
        # 3: if not, what's left to calculate?
        # 4: Calculate the rest
        # 5: Reduce the features down to what the user asked for

        if self.values is not None and len( self.values ) != 0:
            return self

        # Make sure Feature Vector version string is correct, etc:
        #self.Update()

        partial_load = False
        try:
            self.LoadSigFile( quiet=quiet )
            # FIXME: Here's where you'd calculate a small subset of features
            # and see if they match what was loaded from file. The file could be corrupted
            # incomplete, or calculated with different options, e.g., -S1441
            return self
        except IOError:
            # File doesn't exist
            pass
        except WrongFeatureSetVersionError:
            # File has different feature version than desired
            pass
        except IncompleteFeatureSetError:
            # LoadSigFile should create a FeatureComputationPlan
            if not quiet:
                print 'Loaded {0} features from disk for sample "{1}"'.format(
                        len( self.temp_names ), self.name )
            partial_load = True
            pass

        # All hope is lost, calculate features.

        # Use user-assigned feature computation plan, if provided:
        if self.feature_computation_plan != None:
            comp_plan = self.feature_computation_plan

            # I Commented the following out because the computation plan may only reflect
            # the subset of features that haven't been calculated yet:
            # comp_plan.feature_vec_type seems to only contain the minor version
            # i.e., number after the '.'. Assume major version is the latest.
            #self.feature_set_version = '{0}.{1}'.format( 
            #        feature_vector_major_version, comp_plan.feature_vec_type )
        else:
            major, minor = self.feature_set_version.split('.')
            if minor == '0':
                comp_plan = GenerateFeatureComputationPlan( self.feature_names )
            elif minor == '1':
                comp_plan = wndcharm.StdFeatureComputationPlans.getFeatureSet()
            elif minor == '2':
                comp_plan = wndcharm.StdFeatureComputationPlans.getFeatureSetLong()
            elif minor == '3':
                comp_plan = wndcharm.StdFeatureComputationPlans.getFeatureSetColor()
            elif minor == '4':
                comp_plan = wndcharm.StdFeatureComputationPlans.getFeatureSetColorLong()
            else:
                raise ValueError( "Not sure which features you want." )
            self.feature_computation_plan = comp_plan

        # Here are the ImageMatrix API calls:
        # void normalize(double min, double max, long range, double mean, double stddev);
        # int OpenImage(char *image_file_name, int downsample, rect *bounding_rect, double mean, double stddev);
        # void Rotate (const ImageMatrix &matrix_IN, double angle);

        if self.rot is not None:
            raise NotImplementedError( "FIXME: Implement rotations." )

        if self.x is not None and self.y is not None and self.w is not None and self.h is not None:
            bb = wndcharm.rect()
            bb.x = self.x
            bb.y = self.y
            bb.w = self.w
            bb.h = self.h
        else:
            bb = None

        if self.pixel_intensity_mean:
            mean = self.pixel_intensity_mean
            # stddev arg only used in ImageMatrix::OpenImage() if mean is set
            stddev = self.pixel_intensity_stddev
        else:
            # setting mean = 0 is flag to not use mean in ImageMatrix::OpenImage()
            mean = 0
            stddev = 0

        if isinstance( self.source_filepath, str ):
            the_tiff = wndcharm.ImageMatrix()
            if 1 != the_tiff.OpenImage( self.source_filepath, self.downsample, bb, mean, stddev ):
                raise ValueError( 'Could not build an ImageMatrix from {0}, check the path.'.\
                    format( self.source_filepath ) )
        elif isinstance( self.source_filepath, wndcharm.ImageMatrix ):
            if self.downsample or mean:
                raise NotImplementedError( 'still need to implement modifying open pixel plane with downsample, mean or stddev' )
            if not bb:
                the_tiff = self.source_filepath
            else:
                # API calls for copying desired pixels into empty ImageMatrix instance:
                # the_tiff is garbage collected on return
                the_tiff = wndcharm.ImageMatrix()
                # bb only used when calling OpenImage
                the_tiff.submatrix( self.source_filepath, self.x, self.y, self.w, self.h ) # no retval
        else:
            raise ValueError("image parameter 'image_path_or_mat' is not a string or a wndcharm.ImageMatrix")

        # pre-allocate space where the features will be stored (C++ std::vector<double>)
        tmp_vec = wndcharm.DoubleVector( comp_plan.n_features )

        # Get an executor for this plan and run it
        plan_exec = wndcharm.FeatureComputationPlanExecutor( comp_plan )
        plan_exec.run( the_tiff, tmp_vec, 0 )

        # get the feature names from the plan
        comp_names = [ comp_plan.getFeatureNameByIndex(i) for i in xrange( comp_plan.n_features ) ]

        # convert std::vector<double> to native python list of floats
        comp_vals = list( tmp_vec )

        # Feature Reduction/Reorder step:
        # Feature computation may give more features than are asked for by user, or out of order.
        if self.feature_names:
            if self.feature_names != comp_names:
                if partial_load:
                    # If we're here, we've already loaded some but not all of the features
                    # we need. Take what we've already loaded and slap it at the end 
                    # of what was calculated.  Doesn't matter if some of the features are
                    # redundant, because the .index() method returns the first item it finds.
                    # FIXME: if there is overlap between what was loaded and what was 
                    # calculated, check to see that they match.
                    comp_names.extend( self.temp_names )
                    comp_vals.extend( self.temp_values )
                    del self.temp_names
                    del self.temp_values
                self.values = np.array( [ comp_vals[ comp_names.index( name ) ] for name in self.feature_names ] )
        else:
            self.feature_names = comp_names
            self.values = comp_vals

        if not quiet:
            print str( self ), '(calculated ' + str(len(comp_vals)) + ' features)'

        # FIXME: write to disk BEFORE feature reduce
        if write_to_disk:
            self.ToSigFile( quiet=quiet )

        # Feature names need to be modified for their sampling options.
        # Base case is that channel goes in the innermost parentheses, but really it's not
        # just channel, but all sampling options.
        # For now, let the FeatureSpace constructor code handle the modification of feature names
        # for its own self.feature_names
        return self

    #==============================================================
    def CompatibleFeatureSetVersion( self, version ):
        """Note that if either minor version is 0 (i.e not a standard feature vector)
        we return true while in fact, the compatibility is unknown"""

        try:
            version = version.feature_set_version
        except AttributeError:
            # Assume version is the version string then.
            pass

        if self.feature_set_version is None or version is None:
            err_str = "Can't tell if FeatureSpace {0} is compatible with version {1} because "
            if self.feature_set_version is None and version is None:
                err_str += "both are null."
            elif self.feature_set_version is None:
                err_str += "the FS instance's version string is null."
            else:
                err_str += "input version string is null."
            raise AttributeError( err_str.format( self.name, version ) )

        their_major, their_minor = [ int(v) for v in version.split('.',1) ]
        our_major, our_minor = [ int(v) for v in self.feature_set_version.split('.',1) ]

        if( their_major != our_major ):
            return False
        if our_minor and their_minor and our_minor != their_minor:
            return False

        return True

    #==============================================================
    def Normalize( self, reference_features=None, inplace=True, quiet=False ):
        """By convention, the range of feature values in the WND-CHARM algorithm are
        normalized on the interval [0,100]. Normalizing is useful in making the variation 
        of features human readable. Normalized samples are only comprable if they've been 
        normalized against the same feature maxima/minima."""

        if self.normalized_against:
            # I've already been normalized, and you want to normalize me again?
            raise ValueError( "{0} \"{1}\" has already been normalized against {2}.".format (
                self.__class__.__name__, self.name, self.normalized_against ) )

        newdata = {}

        if not reference_features:
            # Specific to FeatureVector implementation:
            # Doesn't make sense to Normalize a 1-D FeatureVector against itself
            # The FeatureSpace implementation of this function has stuff in this block
            err = "Can't normalize {0} \"{1}\" against itself (Normalize() called with blank arg)."
            raise ValueError( err.format( self.__class__.__name__, self.name ) )
        else:
            # Recalculate my feature space according to maxima/minima in reference_features
            if reference_features.feature_names != self.feature_names:
                err_str = "Can't normalize {0} \"{1}\" against {2} \"{3}\": Features don't match.".format(
                  self.__class__.__name__, self.name,
                    reference_features.__class__.__name__, reference_features.name )
                raise ValueError( err_str )
            if not self.CompatibleFeatureSetVersion( reference_features ):
                err_str = 'Incompatible feature versions: "{0}" ({1}) and "{2}" ({3})'
                raise ValueError( err_str.format( self.name, self.feature_set_version,
                    reference_features.name, reference_features.feature_set_version ) )

            if not quiet:
                # Specific to FeatureVector implementation:
                # no num_samples member:
                print 'Normalizing {0} "{1}" ({2} features) against {3} "{4}"'.format(
                    self.__class__.__name__, self.name, len( self.feature_names),
                    reference_features.__class__.__name__, reference_features.name )

            # Need to make sure there are feature minima/maxima to normalize against:
            if not reference_features.normalized_against:
                reference_features.Normalize( quiet=quiet )

            mins = reference_features.feature_minima
            maxs = reference_features.feature_maxima
            newdata['normalized_against'] = reference_features

        newdata['values'] = np.copy( self.values )
        newdata['feature_minima'], newdata['feature_maxima'] = \
            normalize_by_columns( newdata['values'], mins, maxs )

        if inplace:
            return self.Update( **newdata )
        return self.Derive( **newdata )

    #==============================================================
    def FeatureReduce( self, requested_features, inplace=False, quiet=False ):
        """Returns a new FeatureVector that contains a subset of the data by dropping
        features (columns), and/or rearranging columns.

        requested_features := an object with a "feature_names" member
            (FeatureVector/FeatureSpace/FeatureWeights) or an iterable containing
            strings that are feature names.

        Implementation detail: compares input "requested_features" to self.feature_names,
        and "requested_features" becomes the self.feature_names of the returned FeatureVector."""

        try:
            requested_features = requested_features.feature_names
        except AttributeError:
            # assume it's already a list then
            pass

        # Check that self's featurelist contains all the features in requested_features
        selfs_features = set( self.feature_names )
        their_features = set( requested_features )
        if not their_features <= selfs_features:
            missing_features_from_req = their_features - selfs_features
            err_str = "Feature Reduction error:\n"
            err_str += '{0} "{1}" is missing '.format( self.__class__.__name__, self.name )
            err_str += "{0}/{1} features that were requested in the feature reduction list.".format(\
                    len( missing_features_from_req ), len( requested_features ) )
            err_str += "\nDid you forget to convert the feature names into their modern counterparts?"
            raise IncompleteFeatureSetError( err_str )

        # The implementation of FeatureReduce here is similar to FeatureSpace.FeatureReduce
        # Here is where the implementations diverge"
        num_features = len( requested_features )

        if not quiet:
            orig_len = len( self )

        newdata = {}
        newdata[ 'name' ] = self.name + "(feature reduced)"
        newdata[ 'feature_names' ] = requested_features
        newdata[ 'num_features' ] = num_features

        new_order = [ self.feature_names.index( name ) for name in requested_features ]

        # N.B. 1-D version used here, contrast with FeatureSpace.FeatureReduce() implementation.
        newdata[ 'values' ] = self.values[ new_order ]

        if self.feature_maxima is not None:
            newdata[ 'feature_maxima' ] = self.feature_maxima[ new_order ]
        if self.feature_minima is not None:
            newdata[ 'feature_minima' ] = self.feature_minima[ new_order ]

        # If the feature vectors sizes changed then they are no longer standard feature vectors.
        if self.feature_set_version is not None and num_features != self.num_features:
            newdata[ 'feature_set_version' ] = \
                    "{0}.0".format( self.feature_set_version.split('.',1)[0] )

        if inplace:
            newfv = self.Update( **newdata )
        else:
            newfv = self.Derive( **newdata )

        if not quiet:
            print newfv, 'features reduced/reordered from orig len {0}'.format( orig_len )
        return newfv

    #================================================================
    def LoadSigFile( self, sigfile_path=None, quiet=False ):
        """Load computed features from a sig file.

        Desired features indicated by strings currently in self.feature_names.
        Desired feature set version indicated self.feature_set_version.

        Compare what got loaded from file with desired."""

        import re

        if sigfile_path:
            path = sigfile_path
            update_sampling_opts = True
        elif self.auxiliary_feature_storage:
            path = self.auxiliary_feature_storage
            update_sampling_opts = True
        else:
            path = self.GenerateSigFilepath()
            update_sampling_opts = False

        with open( path ) as infile:

            # First, check to see feature set versions match:
            firstline = infile.readline()
            m = re.match( '^(\S+)\s*(\S+)?$', firstline )
            if not m:
                # Deprecate old-style naming support anyway, those features are pretty buggy
                # -CEC 20150104
                raise ValueError( "Can't read a WND-CHARM feature set version from file {0}. File my be corrupted or calculated by an unsupported version of WND-CHARM. Recalculate features and try again.".format( path ) )
                #input_major = 1
                # For ANCIENT sig files, with features calculated YEARS ago
                # Cleanup for legacy edge case:
                # Set the minor version to the vector type based on # of features
                # The minor versions should always specify vector types, but for
                # version 1 vectors, the version is not written to the file.
                #self.feature_set_version = "1." + str(
                #feature_vector_minor_version_from_num_features_v1.get( len( self.values ),0 ) )
                # This is really slow:
                #for i, name in enumerate( names ):
                #retval = wndcharm.FeatureNames.getFeatureInfoByName( name )
                #if retval:
                #    self.feature_names[i] = retval.name
                #else:
                # self.feature_names[i] = name
                # Use pure Python for old-style name translation
                #from wndcharm import FeatureNameMap
                #self.feature_names = FeatureNameMap.TranslateToNewStyle( feature_names )
            else:
                class_id, input_fs_version = m.group( 1, 2 )
                input_fs_major_ver, input_fs_minor_ver = input_fs_version.split('.')
            if self.feature_set_version:
                desired_fs_major_ver, desired_fs_minor_ver = self.feature_set_version.split('.')
                if desired_fs_major_ver != input_fs_major_ver:
                    errstr = 'Desired feature set version "{0}" different from "{1}" in file {2}'
                    raise WrongFeatureSetVersionError(
                            errstr.format( desired_fs_major_ver, input_fs_major_ver, path ) )

            # 2nd line is path to original tiff file, which may be nonsense
            # if sig file was moved post-feature calculation.
            orig_source_tiff_path = infile.readline()
            if self.source_filepath is None:
                from os.path import exists
                # FIXME: Maybe try a few directories?
                if exists( orig_source_tiff_path ):
                    self.source_filepath = orig_source_tiff_path

            # Load data into local variables:
            values, names = \
                zip( *[ line.split( None, 1 ) for line in infile.read().splitlines() ] )

        # Re: converting read-in text to numpy array of floats, np.fromstring is a 3x PIG:
        # %timeit out = np.array( [ float(val) for val in thing ] )
        # 10 loops, best of 3: 38.3 ms per loop
        # %timeit out = np.fromstring( " ".join( thing ), sep=" " )
        # 10 loops, best of 3: 98.1 ms per loop

        # By now we would know by know if there was a sigfile processing error,
        # e.g., file doesn't exist.
        # Safe to set this member now if not already set
        if not self.auxiliary_feature_storage:
            self.auxiliary_feature_storage = path

        # Check to see that the sig file contains all of the desired features:
        if self.feature_names:
            if self.feature_names == names:
                # Perfect! Do nothing.
                pass
            else:
                features_we_want = set( self.feature_names )
                features_we_have = set( names )
                if not features_we_want <= features_we_have:
                    # Need to calculate more features
                    missing_features = features_we_want - features_we_have
                    # create a feature computation plan based on missing features only:
                    self.feature_computation_plan = GenerateFeatureComputationPlan( missing_features )
                    # temporarily store loaded features in temp members to be used by 
                    # self.GenerateFeatures to create the final feature vector.
                    self.temp_names = names
                    self.temp_values = [ float( val ) for val in values ]
                    raise IncompleteFeatureSetError
                else:
                    # If you get to here, we loaded MORE features than asked for,
                    # or the features are out of desired order, or both.
                    values = [ values[ names.index( name ) ] for name in self.feature_names ]
        else:
            # User didn't indicate what features they wanted.
            # It's a pretty dangerous assumption to make that the user just "got 
            # what they wanted" by loading the file, but danger is my ... middle name ;-)
            self.feature_names = list( names )

        self.values = np.array( [ float( val ) for val in values ] )

        # Subtract path so that path part doesn't become part of name
        from os.path import basename
        # Pull sampling options from filename
        path_removed = basename( path )
        self.name = path_removed
        if update_sampling_opts:
            result = self.sig_filename_parser.search( path_removed )
            if result:
                self.Update( **result.groupdict() )

        if not quiet:
            print "Loaded features from file {0}".format( path )
        return self

    #================================================================
    @classmethod
    def NewFromSigFile( cls, sigfile_path, image_path=None, quiet=False ):
        """@return  - An instantiated FeatureVector class with feature names translated from
        the old naming convention, if applicable."""
        return cls( source_filepath=image_path ).LoadSigFile( sigfile_path, quiet )

    #================================================================
    def ToSigFile( self, path=None, quiet=False ):
        """Write features C-WND-CHARM .sig file format

        If filepath is specified, you get to name it whatever you want and put it
        wherever you want. Otherwise, it's named according to convention and placed 
        next to the image file in its directory."""
        from os.path import exists
        if path:
            self.auxiliary_feature_storage = path
        elif self.auxiliary_feature_storage is not None:
            path = self.auxiliary_feature_storage
        else:
            path = self.auxiliary_feature_storage = self.GenerateSigFilepath()

        if not quiet:
            if exists( path ):
                print "Overwriting {0}".format( path )
            else:
                print 'Writing signature file "{0}"'.format( path )
        
        with open( path, "w" ) as out:
            # FIXME: line 1 contains class membership and version
            # Just hardcode the class membership for now.
            out.write( "0\t{0}\n".format( self.feature_set_version ) )
            out.write( "{0}\n".format( self.source_filepath ) )
            for val, name in zip( self.values, self.feature_names ):
                out.write( "{0:0.6g} {1}\n".format( val, name ) )

# end definition class FeatureVector
