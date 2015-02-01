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
    r'(?P<tiling_scheme>-t(?P<tile_num_rows>\d+)(?:x(?P<tile_num_cols>\d+))?_(?P<tile_col_index>\d+)_(?P<tile_row_index>\d+))?',
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
        self.featurenames_list = None
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

        #: self.num_features should always be == len( self.featurenames_list ) == len( self.values )
        self.num_features = None

        # WND-CHARM feature bank-specific params:
        self.color = None
        self.long = None
        self.feature_set_version = None
        self.feature_computation_plan = None

        self.Update( **kwargs )

    #==============================================================
    def __repr__( self ):
        return '"{0}" ({1}) grp {2}, idx {3}, col {4}'.format( self.name, self.label,
                self.samplegroupid, self.samplesequenceid, self.fs_col )

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
            if self.featurenames_list:
                if len( self.featurenames_list ) not in feature_vector_minor_version_from_num_features:
                    minor = 0
                else:
                    # FIXME: If features are out of order, should have a minor version of 0
                    minor = feature_vector_minor_version_from_num_features[ len( self.featurenames_list ) ]
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
            self.feature_set_version = '{0}.{1}'.format( major, minor )
        else:
            major, minor = [ int( val ) for val in self.feature_set_version.split('.') ]

        # reset number of features
        if self.featurenames_list:
            self.num_features = len( self.featurenames_list )
        else:
            if major == 1:
                num_feats_dict = feature_vector_minor_version_from_num_features_v1
            else:
                num_feats_dict = feature_vector_minor_version_from_num_features
            for num_feats, version in num_feats_dict.iteritems():
                if version == minor:
                    self.num_features = num_feats

        # When reading in sampling opts from the path, they get pulled out as strings
        # instead of ints:
        if self.tile_row_index and type( self.tile_row_index ) != int:
            self.tile_row_index = int( self.tile_row_index )
        if self.tile_col_index and type( self.tile_col_index ) != int:
            self.tile_col_index = int( self.tile_col_index )
        if self.tile_num_rows and type( self.tile_num_rows ) != int:
            self.tile_num_rows = int( self.tile_num_rows )
        if self.tile_num_cols and type( self.tile_num_cols ) != int:
            self.tile_num_cols = int( self.tile_num_cols )

        # sequence order:
        # index 0 = position row 0, col 0
        # index 1 = position row 0, col 1
        # index 2 = position row 1, col 0, etc...
        self.samplesequenceid = self.tile_row_index + ( self.tile_num_cols * self.tile_col_index )
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

        # FIXME: sigpaths for FeatureVectors with different channels
        # may have sig file names that will collide/overwrite each other.
        if self.basename:
            base = self.basename
        elif self.source_filepath:
            from os.path import splitext
            base, ext = splitext( self.source_filepath )
        elif self.name:
            base = self.name
            self.basename = base
        else:
            raise ValueError( 'Need for "basename" or "source_filepath" or "name" attribute in FeatureVector object to be set to generate sig filepath.')

        self_namespace = vars(self)
        # ROI:
        roi_params = 'x', 'y', 'w', 'h',
        bb = [self_namespace[key] for key in roi_params]

        if all( [ val is not None for val in bb ] ):
            base += "-B{0}_{1}_{2}_{3}".format( *bb )
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
    def GenerateFeatures( self, write_sig_files_to_disk=True, quiet=True ):
        """@brief Loads precalculated sigs, or calculates new ones, based on which instance
        attributes have been set, and what their values are.
        
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
        self.Update()

        try:
            self.LoadSigFile( quiet=quiet )
            # FIXME: Here's where you'd calculate a small subset of features
            # and see if they match what was loaded from file. The file could be corrupted
            # incomplete, or calculated with different options, e.g., -S1441
            return self
        except IOError:
            pass

        # All hope is lost, calculate features.

        # Use user-assigned feature computation plan, if provided:
        if self.feature_computation_plan != None:
            comp_plan = self.feature_computation_plan
            self.feature_set_version = comp_plan.feature_vec_type
        else:
            major, minor = self.feature_set_version.split('.')
            if minor == '0':
                comp_plan = GenerateComputationPlanFromListOfFeatureStrings( self.featurenames_list )
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

        self_namespace = vars( self )
        bb_members = 'x', 'y', 'w', 'h',
        bb_vals = tuple( [ self_namespace[ var ] for var in bb_members ] )
        if all( [ val is not None for val in bb_vals ] ):
            bb = wndcharm.rect( *bb_vals )
        else:
            bb = None

        if self.pixel_intensity_mean:
            mean = self.pixel_intensity_mean
        else:
            mean = 0

        if self.pixel_intensity_stddev:
            stddev = self.pixel_intensity_stddev
        else:
            stddev = 0

        if isinstance( self.source_filepath, str ):
            the_tiff = wndcharm.ImageMatrix()
            if 1 != the_tiff.OpenImage( self.source_filepath, self.downsample, bb, mean, stddev ):
                raise ValueError( 'Could not build an ImageMatrix from {0}, check the path.'.\
                    format( self.source_filepath ) )
        elif isinstance( self.source_filepath, wndcharm.ImageMatrix ):
            the_tiff = self.source_filepath
            raise NotImplementedError( "FIXME: Still haven't implemented downsample, bounding box, pixel intensity mean and stddev for an already open instance of ImageMatrix." )
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
        if self.featurenames_list:
            if self.featurenames_list != comp_names:
                self.values = [ comp_vals[ comp_names.index( name ) ] for name in self.featurenames_list ]
        else:
            self.featurenames_list = comp_names
            self.values = comp_vals

        if write_sig_files_to_disk:
            self.ToSigFile( quiet=quiet )

        # Feature names need to be modified for their sampling options.
        # Base case is that channel goes in the innermost parentheses, but really it's not
        # just channel, but all sampling options.
        # For now, let the FeatureSpace constructor code handle the modification of feature names
        # for its own self.featurenames_list
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
    def Normalize( self, input_feat_container=None, inplace=True, quiet=False ):
        """By convention, the range of feature values in the WND-CHARM algorithm are
        normalized on the interval [0,100]. Normalizing is useful in making the variation 
        of features human readable. Normalized samples are only comprable if they've been 
        normalized against the same feature maxima/minima."""

        if self.normalized_against:
            # I've already been normalized, and you want to normalize me again?
            raise ValueError( "{0} \"{1}\" has already been normalized against {2}.".format (
                self.__class__.__name__, self.name, self.normalized_against ) )

        newdata = {}

        if not input_feat_container:
            # Specific to FeatureVector implementation:
            # Doesn't make sense to Normalize a 1-D FeatureVector against itself
            # The FeatureSpace implementation of this function has stuff in this block
            err = "Can't normalize {0} \"{1}\" against itself (Normalize() called with blank arg)."
            raise ValueError( err.format( self.__class__.__name__, self.name ) )
        else:
            # Recalculate my feature space according to maxima/minima in input_feat_container
            if input_feat_container.featurenames_list != self.featurenames_list:
                err_str = "Can't normalize {0} \"{1}\" against {2} \"{3}\": Features don't match.".format(
                  self.__class__.__name__, self.name,
                    input_feat_container.__class__.__name__, input_feat_container.name )
                raise ValueError( err_str )
            if not self.CompatibleFeatureSetVersion( input_feat_container ):
                err_str = 'Incompatible feature versions: "{0}" ({1}) and "{2}" ({3})'
                raise ValueError( err_str.format( self.name, self.feature_set_version,
                    input_feat_container.name, input_feat_container.feature_set_version ) )

            if not quiet:
                # Specific to FeatureVector implementation:
                # no num_samples member:
                print 'Normalizing {0} "{1}" ({2} features) against {3} "{4}"'.format(
                    self.__class__.__name__, self.name, len( self.featurenames_list),
                    input_feat_container.__class__.__name__, input_feat_container.name )

            # Need to make sure there are feature minima/maxima to normalize against:
            if not input_feat_container.normalized_against:
                input_feat_container.Normalize( quiet=quiet )

            mins = input_feat_container.feature_minima
            maxs = input_feat_container.feature_maxima
            newdata['normalized_against'] = input_feat_container

        newdata['values'] = np.copy( self.values )
        newdata['feature_minima'], newdata['feature_maxima'] = \
            normalize_by_columns( newdata['values'], mins, maxs )

        if inplace:
            return self.Update( **newdata )
        return self.Derive( **newdata )

    #==============================================================
    def FeatureReduce( self, requested_features, inplace=False ):
        """Returns a new FeatureVector that contains a subset of the data by dropping
        features (columns), and/or rearranging columns.

        requested_features := an object with a "featurenames_list" member
            (FeatureVector/FeatureSpace/FeatureWeights) or an iterable containing
            strings that are feature names.

        Implementation detail: compares input "requested_features" to self.featurenames_list,
        and "requested_features" becomes the self.featurenames_list of the returned FeatureVector."""

        try:
            requested_features = requested_features.featurenames_list
        except AttributeError:
            # assume it's already a list then
            pass

        # Check that self's faturelist contains all the features in requested_features
        selfs_features = set( self.featurenames_list )
        their_features = set( requested_features )
        if not their_features <= selfs_features:
            missing_features_from_req = their_features - selfs_features
            err_str = "Feature Reduction error:\n"
            err_str += '{0} "{1}" is missing '.format( self.__class__.__name__, self.name )
            err_str += "{0}/{1} features that were requested in the feature reduction list.".format(\
                    len( missing_features_from_req ), len( requested_features ) )
            err_str += "\nDid you forget to convert the feature names into their modern counterparts?"
            raise ValueError( err_str )

        # The implementation of FeatureReduce here is similar to FeatureSpace.FeatureReduce
        # Here is where the implementations diverge"
        num_features = len( requested_features )

        newdata = {}
        newdata[ 'name' ] = self.name + "(feature reduced)"
        newdata[ 'featurenames_list' ] = requested_features
        newdata[ 'num_features' ] = num_features

        new_order = [ self.featurenames_list.index( name ) for name in requested_features ]

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
            return self.Update( **newdata )
        return self.Derive( **newdata )
    #================================================================
    def LoadSigFile( self, sigfile_path=None, quiet=False ):

        if sigfile_path:
            path = sigfile_path
        elif self.auxiliary_feature_storage:
            path = self.auxiliary_feature_storage
        else:
            path = self.GenerateSigFilepath()

        with open( path ) as infile:
            lines = infile.read().splitlines()

        # First line is metadata, 2nd line is path to original tiff file (skip that).
        # Process features first, then deal with metadata
        values, names = zip( *[ line.split( None, 1 ) for line in lines[2:] ] )

        # np.fromstring is a 3x PIG:
        # %timeit out = np.array( [ float(val) for val in thing ] )
        # 10 loops, best of 3: 38.3 ms per loop
        # %timeit out = np.fromstring( " ".join( thing ), sep=" " )
        # 10 loops, best of 3: 98.1 ms per loop
        self.values = np.array( [ float( val ) for val in values ] )

        # Now that we know howmany values there are, deal with metadata.
        import re
        self.class_id, self.feature_set_version = \
                re.match( '^(\S+)\s*(\S+)?$' , lines[0] ).group( 1, 2 )
        if self.feature_set_version is None:
            # Cleanup for legacy edge case:
            # Set the minor version to the vector type based on # of features
            # The minor versions should always specify vector types, but for
            # version 1 vectors, the version is not written to the file.
            self.feature_set_version = "1." + str(
                feature_vector_minor_version_from_num_features_v1.get( len( self.values ),0 ) )

        # We would know by know if there was a sigfile processing error,
        # e.g., file doesn't exist.
        # Safe to set this member now if not already set
        if not self.auxiliary_feature_storage:
            self.auxiliary_feature_storage = path

        # This is really slow:
        #for i, name in enumerate( names ):
            #retval = wndcharm.FeatureNames.getFeatureInfoByName( name )
            #if retval:
            #    self.featurenames_list[i] = retval.name
            #else:
            # self.featurenames_list[i] = name

        # Use pure Python for old-style name translation
        #from wndcharm import FeatureNameMap
        #self.featurenames_list = FeatureNameMap.TranslateToNewStyle( featurenames_list )

        # Deprecate old-style naming support anyway, those features are pretty buggy
        # -CEC 20150104
        self.featurenames_list = list( names )

        # Subtract path so that path part doesn't become part of name
        from os.path import basename
        # Pull sampling options from filename
        path_removed = basename( path )
        self.name = path_removed
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
        if not path:
            path = self.GenerateSigFilepath()
        if not quiet:
            if exists( path ):
                print "Overwriting {0}".format( path )
            else:
                print 'Writing signature file "{0}"'.format( path )
        self.auxiliary_feature_storage = path
        with open( path, "w" ) as out:
            # FIXME: line 1 contains class membership and version
            # Just hardcode the class membership for now.
            out.write( "0\t{0}\n".format( self.feature_set_version ) )
            out.write( "{0}\n".format( self.source_filepath ) )
            for val, name in zip( self.values, self.featurenames_list ):
                out.write( "{0:0.6g} {1}\n".format( val, name ) )

# end definition class FeatureVector
