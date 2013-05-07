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
parser.add_argument( '--figure_basename', metavar='<string>', help='Base string used for naming outputted .png graphs',
                     nargs='?', default = None )
parser.add_argument( '--debug', metavar='<optional path>', help='Specify where to send parsing output. Switch specified with argument directs output to filepath specified by argument; no argument means output to STDOUT',
                     nargs='?', default = 'unset' )
args = parser.parse_args()

html_file = args.html_report_filepath[0]

print "parsing HTML file {0}".format( html_file )

if args.figure_basename is None:
	figure_basename = html_file
else:
	figure_basename = args.figure_basename

figure_basename, extension = os.path.splitext( figure_basename )
figure_basename = figure_basename.replace( '.', '_' )

print "using figure basename {0}".format( figure_basename )

text_out = None
if args.debug is 'unset':
	# unset means user doesn't need debug out
	print "User doesn't want to see debug output, sending it to /dev/null"
	text_out = os.devnull
elif args.debug is not None:
	text_out = args.debug
	print "Sending debug output to file {0}".format( text_out )
else:
	print "Sending debug output to STDOUT"


if text_out:
	graph = PredictedValuesGraph.NewFromHTMLFile( html_file, output_filepath=text_out )
else:
	graph = PredictedValuesGraph.NewFromHTMLFile( html_file )

graph.KernelSmoothedDensityGraph( chart_title=html_file )
graph.SaveToFile( figure_basename + '-ks_density' )
graph.RankOrderedPredictedValuesGraph( chart_title=html_file )
graph.SaveToFile( figure_basename + '-rank-ordered' )






