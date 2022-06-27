import astropy
import astropy.units as u
import astroquery
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import numpy as np

# def compare_distances(input_filename, output_filename, dist=___): 
#     """Compares distances of Gaia-Pulsar "matches".

#     Takes in a text file of Gaia-pulsar matches, and filters out the ones with a difference in distances 
#     greater than ___(smt). Returns a text file of the matches within that distance range and each of their 
#     respective distances

#     Args:
#         input_filename (str): Name of the input file with Gaia-Pulsar matches. This could be the output of 
#             get_matches()
#         output_filename (str): Name of the file to which the matches within the distance range are output
#         dist(:obj:'float', optional): Maximum allowed distance between the pulsar and the white dwarf

    
#     """
#     from astropy.table import Table, vstack

#     f = open(input_filename, 'r')

#     for line in f: 
#         values = line.split(',')

# not from parallax, only trust to a factor of three
# with parallax match, trust to 3 sigma 
# only restrict things with a parallax 
# distances should come later on 

# after positional coincidence 

#simple way of propogating errors: take most extreme positions and proper motions -- this will give 1 sigma bounds 