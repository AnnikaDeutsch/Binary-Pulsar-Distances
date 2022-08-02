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

def psr_to_gaia(jname, raj, decj,  pmra, pmdec, posepoch, height, width, radius):
    """Search Gaia for Possible Companion to Pulsar

    Given input parameters read in from a text file following the guidelines of ATNF parameters, 
    queries Gaia DR2 to find matches (nearby objects from Gaia) based on RA and Dec for each object 
    from the text file to within a certain range, the default being 1 arcmin in both ra and dec.

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
        Table: results of the Gaia query in an astropy Table  
    """

    from astropy.time import Time 
    p_ra = raj # comes in as a string of units hh:mm:ss.ss
    p_dec = decj # comes in as a string of units dd:mm:ss.s
    p_pmra = pmra # comes in as a string of units mas/yr
    p_pmdec = pmdec # comes in as a string of units mas/yr
    p_epoch = Time(posepoch, format='mjd').jyear # comes in in units mjd, is immediately converted to jyear tcb

    p_ra_ang = Angle(p_ra + ' hours') # stores the ra in hms as an Angle object
    p_dec_ang = Angle(p_dec + ' degrees') # stores the dec in dms as an Angle object

    p_ra_deg = p_ra_ang.degree * u.deg # stores the ra converted to degrees 
    p_dec_deg = p_dec_ang.degree * u.deg # stores the dec converted to degrees

    # first create variables for the pm units as they come in from atnf, mas/yr
    pmra_masyr = float(p_pmra) * u.mas / u.yr
    pmdec_masyr = float(p_pmdec) * u.mas /u.yr

    # then convert to variable that represent the pms in deg/yr
    pmra_degyr = pmra_masyr.to(u.deg / u.yr)
    pmdec_degyr = pmdec_masyr.to(u.deg / u.yr)

    # to propogate location of pulsar up to gaia time, must calculate epoch difference
    gaia_epoch = Time('2015.5', format='jyear').jyear
    year_diff = (gaia_epoch.tolist() * u.yr) - (p_epoch.tolist() * u.yr) # difference b/w epochs in years

    # get the new ra and dec for the pulsar by updating to gaia epoch 
    p_new_ra = p_ra_deg + (pmra_degyr * year_diff)
    p_new_dec = p_dec_deg + (pmdec_degyr * year_diff)

    print(p_new_ra)
    print(p_new_dec)

    # ra_ext_pos = raj + rajerr + (pmra + pmraerr)*year_diff
    # ra_ext_neg = raj - rajerr + (pmra - pmraerr)*year_diff

    # dec_ext_pos = decj + decjerr + (pmdec + pmdecerr)*year_diff
    # dec_ext_neg = decj - decjerr + (pmdec - pmdecerr)*year_diff


    # Query Gaia within the range of the given pulsar 
    Gaia.ROW_LIMIT = 2000
    coord=SkyCoord(ra=p_new_ra, dec=p_new_dec, unit=(u.degree, u.degree), frame='icrs')
    radius = u.Quantity(radius, u.arcsec)
    width_gaia = u.Quantity(1., u.arcmin) # by default, queries in 1 arcmin range
    height_gaia = u.Quantity(1., u.arcmin) # by default, queries in 1 arcmin range
    j = Gaia.cone_search_async(coordinate=coord, radius=radius)
    results = j.get_results()

    # use python sort function
    
    if len(results) == 0:
        return results
    else:
        results.add_column(jname, name='Companion Pulsar', index=0)
        return results

def get_matches(input_file, output_file, height=1., width=1., radius=1.):
    """Give Gaia matches to Pulsars 

    Takes as input a text file (.csv file) with index number, name, ra, dec, proper
    motion ra, proper motion dec and posepoch of a list of pulsars and produces all of the gaia 
    matches of ra and dec to within a certain range.

    Args: 
        input_file (str): Name of the text file (csv) containing each pulsar with the parameters 'index', 'jname', 
            'ra', 'dec', 'pmra', 'pmdec', 'posepoch' row by row for each object.
        output_file (str): Name of the text file which the pulsar-gaia matches will be output to.
        height (:obj:'float', optional): Height of the rectangle Gaia will query in.
        width (:obj:'float', optional): Width of the rectangle Gaia will query in.
        radius (:obj:'float', optional): Radius of the ractangle Gaia will query in. 


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

      if values[2] == '*' or values[3] == '*' or values[4] == '*' or values[5] == '*' or values[6] == '*':
        skipped += 1
        continue

      # Add result to supertable
      search_result = psr_to_gaia(values[1],values[2],values[3],values[4],values[5],values[6],height,width,radius)
      if (len(search_result) == 0):
        continue
      if first_time:
        results = search_result
        first_time = False
      else:  
        results = vstack([results, search_result])

    results.write(output_file, format='csv', overwrite=True)
    return skipped

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

    from astropy.table import Table, vstack
    import csv

    f = open(input_file, 'r')
    first_time = True
    new = []
    count = 1
    for line in f:
        values = line.split(';')
        while len(values) - 1 > 7:
            values.pop()
        if values[7] != "*":
            values[0] = count
            values[7] = values[7].replace('\n',';')
            if first_time:
                new.insert(0, values)
                first_time = False
            else:
                new.append(values)
            count += 1

# == 'MS' or values[7] == 'NS' or values[7] == 'CO' or values[7] == 'He' or values[7] == 'UL' or values[7] == 'ELL1'
    
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
            parameters 'index', 'jname', 'ra', 'dec', 'pmra', 'pmdec', 'posepoch' and 
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
            if int(ra_err[1]) < 0:
                if int(dec_err[1]) < 0:
                    values.pop(18)
                    values.pop(16)
                    values.pop(14)
                    values.pop(13)
                    values.pop(11)
                    values.pop(10)
                    values.pop(8)
                    values.pop(7)
                    values.pop(5)
                    values.pop(4)
                    values.pop(2)
                    values[0] = count 
                    values[7] = values[7].replace('\n',';')
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
        while len(values) - 1 > 7:
            values.pop()
        g = open('gc_pulsar_names.csv', 'r')
        for line in g:
            check = line.split()
            if values[1] == check[0]:
                not_match = False
        if not_match:
            values[0] = count
            values[7] = values[7].replace('\n','')
            if first_time:
                new.insert(0,values)
                first_time = False
            else:
                new.append(values)
            count += 1
    with open(output_file, 'w') as h:
        write = csv.writer(h, delimiter=';')
        write.writerows(new)


class TestBinaryCheck:

    def test_eliminates_all_not_binaries(self):
        """
        Tests that, given a small file with a known number of pulsars not in binaries,
        the function check_binary() eliminates all of the isolated pulsars and keeps
        all of the binary pulsars.
        """

        input_file = 'e1_input.csv'
        output_file = 'e1_output.csv'
        check_binary(input_file, output_file)

        count = 0
        f = open(output_file,'r')
        for line in f:
            count += 1

        assert count == 2

class TestPosUncertaintyCheck:

    def test_pos_uncerainty_check(self):
        """
        Tests that, given a small file (e2_input.csv) with 2 acceptable and 2 unacceptable position
        uncertainty, the resulting output file has only 2 entries.
        """

        input_file = 'e2_input.csv'
        output_file = 'e2_output.csv'
        check_pos_uncertainty(input_file, output_file)

        count = 0
        f = open(output_file,'r')
        for line in f:
            count += 1
        
        assert count == 2

class TestInGlobularCheck:

    def test_in_globular_check(self):
        """
        Tests that, given a file with 2 pulsars not in globular clusters and 1 pulsar in a 
        globular cluster, the 1 pulsar in the globular is removed and a file with the 2 remaining 
        pulsars is created.
        """

        input_file = 'e3_input.csv'
        output_file = 'e3_output.csv'

        check_in_globular(input_file, output_file)

        count = 0
        f = open(output_file, 'r')
        for line in f:
            count += 1

        assert count == 2


# check_pos_uncertainty('initial_input.csv', 'i1.csv')
# check_binary('i1.csv', 'i2.csv')
# check_in_globular('i2.csv', 'almost_final.csv')
# get_matches('almost_final.csv', 'final.csv')

# check_pos_uncertainty('all_atnf.csv', 'a1.csv')
# check_binary('a1.csv', 'a2.csv')
# check_in_globular('a2.csv', 'all_almost_final.csv')
# get_matches('all_almost_final.csv', 'all_final.csv')


# script to compare pmra and pmdec of the identified matches 
# many of the confirmed literature matches do not have great agreement in their proper motions, so 
# I realized narrowing them down that way is not particularly helpful...

# a = open('all_final.csv', 'r') # read output file of get_mathces() that has pms of Gaia objects
# first_time = True
# new = []
# top = True

# for line in a:
#     values = line.split(',')
#     if top:
#         top = False
#         continue
#     b = open('all_almost_final.csv', 'r') # read input file of get_matches() that has pms of pulsars
#     gaia_pmra = float(values[13])
#     gaia_pmdec = float(values[15])
#     for line in b:
#         check = line.split(';')
#         if values[0] == check[1]:
#             psr_pmra = float(check[4])
#             psr_pmdec = float(check[5])
#             if int(gaia_pmra) == int(psr_pmra) and int(gaia_pmdec) == int(psr_pmdec):
#                 values.pop()
#                 if first_time:
#                     new.insert(0,values)
#                     first_time = False
#                 else:
#                     new.append(values)


# import csv
# with open('all_binary_matches.csv', 'w') as h:
#     write = csv.writer(h, delimiter=';')
#     write.writerows(new)


# script to create an updated version of all_final.csv with only the parameters we wanna look at

c = open('all_final.csv', 'r')
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

import csv
with open('all_final_short.csv', 'w') as h:
    write = csv.writer(h, delimiter=';')
    write.writerows(new)




# script to make a csv of just the NAMES of the pulsars in Globular Clusters, from a file that had a lot
# more info

# import csv
# import pandas as pd
# j = open('gc_pulsars.csv','r')
# keep = []
# first_time = True
# for line in j:
#     values = line.split()
#     if first_time:
#         keep.insert(0,values[0])
#     else:
#         keep.append(values[0])

# df = pd.DataFrame(keep)
# df.to_csv('gc_pulsar_names.csv', index=False)

# with open('gc_pulsar_names.csv', 'w') as g:
#     write = csv.writer(g)
#     write.writerows(keep)