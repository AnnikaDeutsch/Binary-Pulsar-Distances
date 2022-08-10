import astropy
import astropy.units as u
import astroquery
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import quantity_support
from astropy.time import Time 
import os
from astropy.io.votable import parse_single_table
from astropy.time import Time
import pytest
from astropy.table import Table, vstack
import csv
import pandas as pd


# query gaia for one of the matches, and only retrieve the parameters we care about, in this case parallax
# pm, and g-band magnitude
job = Gaia.launch_job("select top 20 "
                      "gaia_source.source_id, gaia_source.parallax, gaia_source.parallax_error,"
                      " gaia_source.pmra,gaia_source.pmdec,gaia_source.phot_g_mean_mag "
                      "from gaiadr2.gaia_source "
                      "where (gaiadr2.gaia_source.source_id=44308738051547264)")

# store the results of the query in an astropy table
gaia_results = job.get_results()

# query atnf for the attributes we want of the pulsars 
# query = QueryATNF(params=['JName', 'PMRA', 'PMDec', 'BINARY', 'PX', 'DM'], psrs= ['J0337+1715'])
# atnf_results = query.table


# create a txt file to store the data in a nice format
g = open('psrJ0337+1715.txt', 'w')
g.write('PSR JNAME: ' + 'J0337+1715\n')
g.write('\n----Gaia Match Attributes----\n')
g.write('Gaia Source ID: ' + str(gaia_results[0][0]) + '\n')
g.write('Gaia PMRA: ' + str(gaia_results[0][3]) + ' mas/yr\n')
g.write('Gaia PMDEC: ' + str(gaia_results[0][4]) + ' mas/yr\n')
g.write('Gaia Parallax: ' + str(gaia_results[0][1]) + ' \u00B1 ' + str(gaia_results[0][2]) + ' mas')
g.write('\nGaia G-band-mag: ' + str(gaia_results[0][5]) + '\n')
g.write('\n----Pulsar Attributes----\n')
g.write('ATNF DM: ' + str(atnf_results['DM'][0]) + ' \u00B1 ' + str(atnf_results['DM_ERR'][0]) + ' pc / cm^{3}\n')
g.write('ATNF PMRA: ' + str(atnf_results['PMRA'][0]) + ' \u00B1 ' + str(atnf_results['PMRA_ERR'][0]) + 
        ' mas/yr\n')
g.write('ATNF PMDEC: ' + str(atnf_results['PMDEC'][0]) + ' \u00B1 ' + str(atnf_results['PMDEC_ERR'][0]) + 
        ' mas/yr\n')
g.write('ATNF Parallax: ' + str(atnf_results['PX'][0]) + ' \u00B1 ' + str(atnf_results['PX_ERR'][0]) + 
        ' mas\n')
g.close()

# now write this into a function to do it for any given pulsar
def pretty_print(pulsar, source_id, filename, no_glob):
    """Prints attributes of the pulsar/match pair in a readable format.

    Args:
        pulsar (str): JName of the pulsar as it is in the list of matched files from matching_pipeline()
        source_id (str): Source ID of the Gaia match as given from matching_pipeline()
    
    """

    # query gaia for one of the matches, and only retrieve the parameters we care about, in this case parallax
    # pm, and g-band magnitude
    job = Gaia.launch_job("select top 20 "
                          "gaia_source.source_id, gaia_source.parallax, gaia_source.parallax_error,"
                          " gaia_source.pmra, gaia_source.pmdec, gaia_source.phot_g_mean_mag "
                          "from gaiadr3.gaia_source "
                          "where (gaiadr3.gaia_source.source_id=" + source_id + ")")
    print(source_id)
    # store the results of the query in an astropy table
    gaia_results = job.get_results()
    print(gaia_results)

    # query atnf for the attributes we want of the pulsars 
    l = open(no_glob, 'r')
    for line in l:
        values = line.split(';')
        if values[1] == pulsar:
                dm = values[9]
                pmra = values[4]
                pmra_err = values[5]
                pmdec = values[6]
                pmdec_err = values[9]

    # create a txt file to store the data in a nice format
    g = open(filename+'.txt', 'a')
    g.write('PSR JNAME: ' + pulsar + '\n')
    g.write('\n----Gaia Match Attributes----\n')
    g.write('Gaia Source ID: ' + str(gaia_results[0][0]) + '\n')
    g.write('Gaia PMRA: ' + str(gaia_results[0][3]) + ' mas/yr\n')
    g.write('Gaia PMDEC: ' + str(gaia_results[0][4]) + ' mas/yr\n')
    g.write('Gaia Parallax: ' + str(gaia_results[0][1]) + ' \u00B1 ' + str(gaia_results[0][2]) + ' mas')
    g.write('\nGaia G-band-mag: ' + str(gaia_results[0][5]) + '\n')
    g.write('\n----Pulsar Attributes----\n')
    g.write('ATNF DM: ' + dm + 'pc / cm^{3}\n')
    g.write('ATNF PMRA: ' + pmra + ' \u00B1 ' + pmra_err + ' mas/yr\n')
    g.write('ATNF PMDEC: ' + pmdec + ' \u00B1 ' + pmdec_err + 
            ' mas/yr\n')
    g.write('\n')
    g.write('-------------------------------------------------------------')
    g.write('\n')
    g.close()

pretty_print('J0348+0432', '3273288485744249344', 'J0348+0432') # cool it works, time to put it in the matching pipeline

