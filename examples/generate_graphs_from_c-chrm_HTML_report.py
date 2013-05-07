#!/usr/bin/env python
"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                               
 Copyright (C) 2013 National Institutes of Health 

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
 Written by:  Christopher Coletta <christopher.coletta [at] nih [dot] gov>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Takes a C-chrm-generated HTML report and generates a kernel-smoothed probability
density estimate graph or a rank-ordered predicted values graph"""


# import pychrm
from pychrm.FeatureSet import *
from pychrm import __version__ as pychrm_version
print "pychrm "+pychrm_version

import argparse
import os

parser = argparse.ArgumentParser( description="Takes a C-chrm-generated HTML report and generates a kernel-smoothed probability density estimate graph or a rank-ordered predicted values graph")
parser.add_argument( 'html_report_filepath', help='path WND-CHARM HTML report',
                     nargs=1 )
parser.add_argument( 'graph_output_filepath', help='Base string used for naming outputted .png graphs',
                     nargs='?', default = None )
parser.add_argument( 'debug', help='Name of file into which HTML parsing output and aggregation output will go.',
                     nargs='?', default = None )
args = parser.parse_args()

html_file = args.html_report_filepath[0]

print "parsing HTML file {0}".format( html_file )

if args.graph_output_filepath is None:
	graph_basename = html_file
else:
	graph_basename = args.graph_output_filepath

graph_basename, extension = os.path.splitext( graph_basename )
graph_basename = graph_basename.replace( '.', '_' )

print "using graph basename {0}".format( graph_basename )

if args.debug is None:
	text_out = os.devnull
else:
	text_out = args.debug

graph = PredictedValuesGraph.NewFromHTMLFile( html_file, output_filepath=text_out )
graph.KernelSmoothedDensityGraph( chart_title=html_file )
graph.SaveToFile( graph_basename + '-ks_density' )
graph.RankOrderedPredictedValuesGraph( chart_title=html_file )
graph.SaveToFile( graph_basename + '-rank-ordered' )






