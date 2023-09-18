"""
Created June 2022

@author: Annika Deutsch
@date: 06/2022
@title: eliminating.py
@description: This python script holds all the function definitions of the package PSRmatch.
"""


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


def psr_to_gaia(jname, raj, decj,  pmra, pmdec, posepoch, binary, bincomp, radius):
    """Searches Gaia for possible companion to any given pulsar

    Given input ATNF parameters for a single pulsar, 
    query Gaia to find Gaia objects near that pulsar based on RA and DEC. Returns 
    Gaia query results in a table and writes them to a text file

    Args:
        jname (str): Name of the pulsar being checked for matches 
        raj (str): Right ascension of the pulsar in hh:mm:ss.ss format -- is this actually what it is?
        decj (str): Declination of the pulsar in degrees:mm:ss.ss format -- is this actually what it is?
        pmra (str): Proper motion in ra of the pulsar in string format and mas/yr units
        pmdec (str): proper motion in dec of the pulsar in string format and mas/yr units
        posepoch (str): epoch that the data was taken in string format and mjd units
        height (float): Height of Gaia box search in arcminutes
        width (float): Width of Gaia box search in arcminutes 
        radius (float): Radius of the Gaia cone search in arcminutes
    
    Returns:
        results (Table): results of the Gaia query in an astropy Table  
    """

    from astropy.time import Time 
    p_pmra = u.Quantity(pmra, u.mas/u.yr) # comes in as a string of units mas/yr
    p_pmdec = u.Quantity(pmdec, u.mas/u.yr) # comes in as a string of units mas/yr
    p_epoch = Time(posepoch, format='mjd').jyear # comes in in units mjd, is immediately converted to jyear tcb

    psr = SkyCoord(ra= raj, dec= decj, unit= (u.hourangle, u.deg), frame = 'icrs', pm_ra_cosdec= p_pmra, pm_dec= p_pmdec)

    # calculate epoch difference
    gaia_epoch = Time('2016.0', format='jyear').jyear
    year_diff = (gaia_epoch.tolist() * u.yr) - (p_epoch.tolist() * u.yr) # difference b/w epochs in years

    # pm times the time diff
    pmtransra = psr.pm_ra_cosdec * year_diff
    pmtransdec = psr.pm_dec * year_diff

    ra_dr3 = psr.ra.to_value(u.mas) * u.mas + (pmtransra)
    dec_dr3 = psr.dec.to_value(u.mas) * u.mas + (pmtransdec)

    # new SkyCoord object with the propagated positions
    psr_dr3 = SkyCoord(ra= ra_dr3, dec= dec_dr3, unit= (u.mas, u.mas), frame= 'icrs', pm_ra_cosdec= p_pmra, pm_dec= p_pmdec)

    # Query Gaia within the range of the given pulsar 
    Gaia.ROW_LIMIT = 2000
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source" # Select eDR3
    j = Gaia.cone_search_async(coordinate=psr_dr3, radius=u.Quantity(radius, u.arcsec))
    results = j.get_results()

    # use python sort function
    if len(results) == 0:
        # with open('/home/annika_deutsch/Binary-Pulsar-Distances/Binary_Pulsar_Distances/matches_10arcsec.csv', 'a') as fd:
        #     fd.write(jname + ',\n')
        #     fd.close()
        return results
    else:
        results.add_column(jname, name='Companion Pulsar', index=0)
        results.add_column(raj, name='Pulsar RA', index=1)
        results.add_column(decj, name='Pulsar DEC', index=2)
        results.add_column(binary, name='Binary', index=3)
        results.add_column(bincomp, name='Binary Companion', index=4)
        results.write('temp.csv', overwrite=True) #writes the results of a single query to a csv file
        # f = open('temp.csv', 'r')
        # for line in f:
        #     with open('/home/annika_deutsch/Binary-Pulsar-Distances/Binary_Pulsar_Distances/matches_10arcsec.csv', 'a') as fd:
        #         fd.write(line)
        #         fd.close()
        return results
    

def psr_to_gaia_nominal(jname, raj, decj, radius):
    """
        Looks for Gaia matches to a single pulsar within a specified radius, and returns the table resulting 
        from the Gaia query with the psr name, ra and dec. A simpler version of psr_to_gaia().

        Args:
            jname (str): Name of the pulsar being checked for matches 
            raj (str): Right ascension of the pulsar in hh:mm:ss.ss format -- is this actually what it is?
            decj (str): Declination of the pulsar in degrees:mm:ss.ss format -- is this actually what it is?
            radius (float): Radius of the Gaia cone search in arcminutes
    
        Returns:
            Table: results of the Gaia query in an astropy Table 
    """
    
    # Query Gaia within the range of the given pulsar 
    Gaia.ROW_LIMIT = 2000
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source" # Select early Data Release 3
    sc = SkyCoord(ra= raj, dec= decj, frame= 'icrs', unit=(u.hourangle, u.deg))
    j = Gaia.cone_search_async(coordinate=sc, radius=u.Quantity(radius, u.arcsec))
    results = j.get_results()

    # use python sort function
    
    if len(results) == 0:
        return results
    else:
        results.add_column(jname, name='Companion Pulsar', index=0)
        results.add_column(raj, name='Pulsar RA', index=1)
        results.add_column(decj, name='Pulsar DEC', index=2)
        return results



def get_matches(input_file, radius=1.):
    """Takes a list of pulsars and returns a list of those with potential matches.

    Takes as input a text file containing parameters of a list of pulsars, and produces all of the gaia 
    matches of ra and dec to within a certain range.

    Args: 
        input_file (str): Name of the text file (csv) containing each pulsar with the parameters 'index', 'jname', 
            'ra', 'dec', 'pmra', 'pmdec', 'posepoch' row by row for each object.
        radius (:obj:'float', optional): Radius of the circle Gaia will query in. 

    Returns:
        results (Table): Results of Gaia queries for all pulsars in input list, as an astropy table
        hits (int): Total number of Gaia objects returned for the pulsar list cross-matched to a certain radius
        skipped (int): Total number of pulsars skipped in the input list due to lack of information

    """
    from astropy.table import Table, vstack

    f = open(input_file, "r")
    results = Table()
    first_time = True
    skipped = 0

    # Loop through file of ATNF data and combine tables of Gaia matches into one supertable
    for line in f:

      # Parse input
      values = line.split(';')
      
      jname = values[1]
      ra = values[3]
      ra_err = values[4]
      dec = values[6]
      dec_err = values[7]
      pmra = values[9]
      pmra_err = values[10]
      pmdec = values[12]
      pmdec_err = values[13]
      posepoch = values[15]
      binary = values[17]
      bincomp = values[19]

      if ra == '*' or ra_err == '*' or dec == '*' or dec_err == '*' or pmra == '*' or pmra_err == '*' or pmdec == '*' or pmdec_err == '*' or posepoch == '*':
        skipped += 1
        # with open('/home/annika_deutsch/Binary-Pulsar-Distances/Binary_Pulsar_Distances/matches_10arcsec.csv', 'a') as fd:
        #   fd.write(values[1] + ',\n')
        #   fd.close()
        continue

      # write a condition that will perform the query if that pulsar is not in the table, and will skip if it is
      query = True
      g = open('/home/billee/Binary-Pulsar-Distances/Binary_Pulsar_Distances/matches_10arcsec.csv', 'r')
      for line in g:
        name = line.split(',')
        if name[0] == jname:
          query = False
          break
      
      if query:
        # Add result to supertable
        search_result = psr_to_gaia(jname,ra,dec,pmra,pmdec,posepoch,binary,bincomp,radius)
        if (len(search_result) == 0):
          continue
        if first_time:
          results = search_result
          first_time = False
        else:  
          results = vstack([results, search_result])

    hits = len(results)
    # results.write(output_file, format='csv', overwrite=True)
    return results, hits, skipped


def check_binary(input_file, output_file):
    """Removes isolated pulsars from input file.

    Given an input file, checks that each pulsar in the file is in a known binary, and 
    removes those that are not. Checks the 'type' parameter for each entry, and removes 
    those objects with an asterisk as the value in that field.

    Args:
        input_file (str): Name of the text file (csv) containing each pulsar, with the 
            parameters 'index', 'jname', 'ra', 'dec', 'pmra', 'pmdec', 'posepoch' and 
            'binary' as listed in the ATNF catalouge.
        output_file (str): Name of the new text file (csv) created with only the binary 
            pulsars.
    

    """
    f = open(input_file, 'r')
    first_time = True
    new = []
    count = 1
    for line in f:
        values = line.split(';')
        while len(values) - 1 > 11:
            values.pop()
        if values[11] != "*":
            values[0] = count
            values[11] = values[11].replace('\n',';')
            if first_time:
                new.insert(0, values)
                first_time = False
            else:
                new.append(values)
            count += 1
    
    with open(output_file, 'w') as g:
        write = csv.writer(g, delimiter=';')
        write.writerows(new)


def check_pos_uncertainty(input_file, output_file):
    """Removes pulsars with position uncertainty > 1".

    Given an input file, checks the uncertainty in ra and dec of each pulsar, and removes 
    any object that has an uncertainty greater than >1", to exclude pulsars that cannot
    be confidently matched against Gaia astrometry. The input file MUST be of the ATNF 
    output style "long with errors", and then copied directly into a csv for the function 
    to work correctly.

    Args:
        input_file (str): Name of the text file (csv) containing each pulsar, with the 
            parameters 'index', 'jname', 'ra', 'dec', 'pmra', 'pmdec', 'posepoch', 'DM', and 
            'binary' as listed in the ATNF catalouge.
        output_file (str): Name of the new text file (csv) created with only the binary 
            pulsars.
    
    """

    import csv

    f = open(input_file, 'r')
    first_time = True
    new = []
    count = 1
    index = 0
    for line in f:
        # if (index) % 5 == 0:
        #     continue
        if (index%5) != 0 or index == 0:
            values = line.split()
            # print(values)
            # while len(values) - 1 > 7:
            #     values.pop()
            if values[4] != '0':
                ra_err = values[4].split('e')
            else:
                ra_err = [0,-1]
            if values[7] != '0':
                dec_err = values[7].split('e')
            else:
                ra_err = [0,-1]
            print(ra_err)
            if int(ra_err[1]) < 0:
                if int(dec_err[1]) < 0:
                    values.pop(21)
                    values.pop(19)
                    values.pop(16)
                    values.pop(14)
                    values.pop(11)
                    values.pop(8)
                    values.pop(7)
                    values.pop(5)
                    values.pop(4)
                    values.pop(2)
                    values[0] = count 
                    # values[9] = values[9].replace('\n',';')
                    values.append('')
                    if first_time:
                        new.insert(0, values)
                        first_time = False
                    else:
                        new.append(values)
                    count += 1
            index += 1
        else:
            index = 0
    

    
    with open(output_file, 'w') as g:
        write = csv.writer(g, delimiter=';')
        write.writerows(new)


def check_in_globular(input_file, output_file):
    """Removes pulsars in globular clusters.

    Given an input file, checks the pulsar names against a list of pulsars in globular clusters, stored in 
    a file called 'gc_pulsar_names.csv', and removes those that match with any of the names in the list. Make
    sure the file 'gc_pulsar_names.csv' from the GitHub repo is in the directory in which you are running
    this function!

    Args:
        input_file (str): Name of the text file (csv) containing each pulsar on a separate line, with the 
            parameters 'index', 'jname', 'ra', 'dec', 'pmra', 'pmdec', 'posepoch', 'binary' as listed in ATNF
        output_file (str): Name of the new text file (csv) created with only pulsars not in globular 
            clusters.
    
    """
    import csv

    f = open(input_file, 'r')
    first_time = True
    new = []
    count = 1
    for line in f:
        not_match = True
        values = line.split(';')
        while len(values) - 1 > 9:
            values.pop()
        g = open('/home/annika_deutsch/Binary-Pulsar-Distances/Binary_Pulsar_Distances/gc_pulsar_names.csv', 'r')
        for line in g:
            check = line.split()
            if values[1] == check[0]:
                not_match = False
        if not_match:
            values[0] = count
            #values[9] = values[9].replace('\n','')
            if first_time:
                new.insert(0,values)
                first_time = False
            else:
                new.append(values)
            count += 1
    with open(output_file, 'w') as h:
        write = csv.writer(h, delimiter=';')
        write.writerows(new)


def gc_psr_names():
    """
    Gets the names of the pulsars in globular clusters from 'gc_pulsars.csv' and
    saves them in the csv file 'gc_pulsar_names.csv'.
    """
    j = open('gc_pulsars.csv','r')
    keep = []
    first_time = True
    for line in j:
        values = line.split()
        if first_time:
            keep.insert(0,values[0])
        else:
            keep.append(values[0])

    df = pd.DataFrame(keep)
    df.to_csv('gc_pulsar_names.csv', index=False)

    with open('gc_pulsar_names.csv', 'w') as g:
        write = csv.writer(g)
        write.writerows(keep)



def pretty_print(pulsar, source_id, filename, no_glob):
    """Prints attributes of the pulsar/match pair in a readable format.

    Args:
        pulsar (str): JName of the pulsar as it is in the list of matched files from matching_pipeline().
        source_id (str): Source ID of the Gaia match as given from matching_pipeline().
        filename (str): Name of the (txt) file to which the results will be output.
        no_glob (str): name of the file from the matching pipeline that contains all the pulsars from 
            after the step of removing pulsars in globular clusters. 
    
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



def matching_pipeline(input_file, output_file, no_pos, no_bin, no_glob, match_all_params, radius=1., 
                      pretty_print= False, remove_bin= False):
    """Carries out the identification of sources from a list of pulsars to a list potential Gaia matches.

    Takes as input a text file of the ATNF format "long with errors", and runs through each step of the source
    identification to produce a list of potential Gaia matches and some or all of their associated Gaia 
    parameters. Pretty print allows the output to be printed in a more readable 
    and user-friendly format. 

    Args:
        input_file (str): Path to the text file with the original list of pulsars with the paramters 
            'index', 'jname', 'ra', 'dec', 'pmra', 'pmdec', 'posepoch', 'DM', and 
            'binary' as listed in the ATNF catalouge.
        output_file (str): Path to the final text file (csv) with the list of potential Gaia matches and 
            a smaller amount of their associated parameters.
        no_pos (str): Path to the intermediate file containing only pulsars with position uncertainties 
            <1" as it comes out after being run through check_pos_uncertainty().
        no_bin (str): Path to the intermediate file containing only pulsars in binaries
            as it comes out after being run through check_binary().
        no_glob (str): Path to the intermediate file containing only pulsars not in globular clusters
            as it comes out after being run through check_in_globular().
        match_all_params (str): Path to the intermediate file containing the list of all potential 
            Gaia matches and all of their associated parameters.
        radius (:obj:'float', optional): Radius of the circle Gaia will query in. 
        pretty_print (:obj:'bool', optional): Will print an additional file of the list in a more user-
            friendly format if set to True.
    
    """
    check_pos_uncertainty(input_file, no_pos)
    if remove_bin: 
        check_binary(no_pos, no_bin)
        check_in_globular(no_bin, no_glob)
        get_matches(no_glob, match_all_params, radius= radius)
    else: 
        check_in_globular(no_pos, no_glob)
        get_matches(no_glob, match_all_params, radius= radius)

    # script to create an updated version of all_final.csv with only the parameters we wanna look at

    c = open(match_all_params, 'r')
    new = []
    first_time = True
    for line in c:
        values = line.split(',')
        values.pop(95)
        count = 0 
        while count < 40:
            values.pop(91-count)
            count += 1
        index = 0
        while index < 30:
            values.pop(50-index)
            index += 1
        values.pop(12)
        values.pop(5)
        values.pop(4)
        values.pop(2)
        values.pop(1)
        values.pop()
        if first_time:
            new.insert(0,values)
            first_time = False
        else:
            new.append(values)
    print(new)
    import csv
    with open(output_file, 'w') as h:
        write = csv.writer(h, delimiter=';')
        write.writerows(new)

    if pretty_print:
        f = open(match_all_params, 'w')
        first_time = True
        filename = output_file + 'prettyprinted'
        no_glob = '/home/annika_deutsch/Binary-Pulsar-Distances/text_files_test/matching_pipeline_test/t1_noglob.csv'
        for line in f:
            values = line.split(';')
            if first_time:
                first_time = False
                continue
            else:
                pretty_print(values[0], values[1], filename, no_glob)
