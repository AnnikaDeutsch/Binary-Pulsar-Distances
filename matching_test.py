# %%
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

def imoprting():
    """
    Imports every important package needed to run this package.
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

# %%
# inputs: ra, dec, height_of_rectangle, width_of_rectangle, radius_of_circle

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
    radius = u.Quantity(radius, u.arcmin)
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

# %%

def get_matches(input_file, output_file, height=1., width=1., radius=1.):
    """Give Gaia matches to Pulsars 

    Takes as input a text file (.csv file) with index number, name, ra, dec, proper
    motion ra, proper motion dec and posepoch of a list of pulsars and produces all of the gaia 
    matches of ra and dec to within a certain range.

    Args: 
        input_file (str): Name of the text file (csv) containing each pulsar with the parameters 'index', 'name', 
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

    # Loop through file of ATNF data and combine tables of Gaia matches into one supertable
    for line in f:

      # Parse input
      values = line.split(';')

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
# %% [markdown]
# Below is a number of unit tests on psr_to_gaia(), testing edge cases, cases that should throw an error, and 
# that cases with already known outcomes return the correct values

class TestPsrToGaia:
    def test_returns_empty_table_when_range_is_zero(self):
        """psr_to_gaia() returns empty table when height and width equal zero

        Tests that, when given an otherwise valid dataset, if height and width are zero, an empty astropy Table
        is returned.
        """

        jname = 'J1012+5307'
        raj = '10:12:33.4'
        decj = '+53:07:02.2'
        pmra = '9.240'
        pmdec = '1.770'
        posepoch = '56000.00'
        height = 0.
        width = 0.
        radius = 0.

        check = psr_to_gaia(jname, raj, decj, pmra, pmdec, posepoch, height, width, radius)
        assert len(check) == 0

    def test_minus_sign_in_decj_and_pm_strings(self):
        """
        Tests that the function psr_to_gaia() still behaves properly with string inputs that are intended to be 
        interpreted as negative values, such as decj and pmdec
        """

        jname = 'J0024-7204J'
        raj = '00:23:59.4'
        decj = '-72:03:58.7'
        pmra = '5.270'
        pmdec = '-3.590'
        posepoch = '51600.00'
        height = 1.
        width = 1.
        radius = 1.

        check = psr_to_gaia(jname, raj, decj, pmra, pmdec, posepoch, height, width, radius)
        astropy_table = parse_single_table('update_test_J0024-7204_matches.vot.gz').to_table(use_names_over_ids=True)
        condition = False
        for val in check['source_id']:
            if val == astropy_table['source_id'][0]:
                condition = True
                break
        assert len(check['ra']) == len(astropy_table['ra'])
        assert condition

    def test_dec_between_zero_and_minus_one(self):
        """
        Tests that the function psr_to_gaia() still gives the correct output when the declination is between 
        0 degrees and -1 degrees, where the declination conversions are often messed up; this should not be the 
        case for me as I used astropy quantities, however it is important to check.
        """

        jname = 'J1607-0032'
        raj = '16:07:12.0'
        decj = '-00:32:41.5'
        pmra = '-26.470'
        pmdec = '-27.500'
        posepoch = '56000.00'
        height = 1.
        width = 1. 
        radius = 1. 

        check = psr_to_gaia(jname, raj, decj, pmra, pmdec, posepoch, height, width, radius)
        check.write('check', format='csv', overwrite=True)
        astropy_table = parse_single_table('J1607-0032_matches.vot.gz').to_table(use_names_over_ids=True)
        assert len(check['ra']) == len(astropy_table['ra'])

    def test_error_thrown_when_invalid_ra(self):
        """
        Tests that the function psr_to_gaia() throws an error when an invalid right ascension 
        is input.
        """

        jname = 'J0348+0432'
        raj = '25:48:43.6'
        decj = '+04:32:11.4'
        pmra = '4.040'
        pmdec = '3.500'
        posepoch = '56000.00'
        height = 1.
        width = 1.
        radius = 1.

        with pytest.raises(Exception):
            psr_to_gaia(jname, raj, decj, pmra, pmdec, posepoch, height, width, radius)

class TestGetMatches:

    def test_more_than_twenty(self):
        """
        Tests whether the function get_matches returns a csv of all the pulsars in a given
        file. In particular, if there is an input file with more than 20 pulsars, it 
        still returns the matches for each pulsar 

        Note: Do not change or remove the file 'get_match_check_input.csv' as it is 
        important to the proper functionality of this test
        """

        input_file = 'get_match_check_input.csv'
        output_file = 'get_match_check_output.csv'

        get_matches(input_file, output_file)

        f = open(output_file, 'r')

        for line in f:
            pass
        last_line_output = line.split(',')

        g = open(input_file, 'r')

        for line in g:
            pass
        last_line_input = line.split(';')

        assert last_line_output[0] == last_line_input[1]


    def test_gaia_and_gm_return_same_results(self):
        """
        Tests that the function get_matches() returns the same exact table of results
        when only one pulsar is input, as does a direct Gaia query
        """

        input_file = 'gm1_input.csv'
        output_file = 'gm1_output.csv'

        get_matches(input_file, output_file)

        f = open(output_file, 'r')
        count = 0
        for line in f:
            count += 1

        g = open('gm1_fromgaia_csv.csv', 'r')
        index = 0
        for line in g:
            index += 1

        assert count == index


    def test_blank_file_returned_when_blank_file_input(self):
        """
        Test that a blank csv is returned when a blank csv is input by get_matches()
        """
        
        input_file = 'gm2_input.csv'
        output_file = 'gm2_output.csv'

        get_matches(input_file, output_file)

        f = open(output_file, 'r')
        count = 0
        for line in f:
            values = line.split(',')
            if values[0] != '':
                count += 1
        
        assert count == 1


    def test_error_thrown_when_wrong_input_params(self):
        """
        Tests that get_matches() throws an error when the csv input file has parameters
        that are incorrect
        """

        input_file = 'gm3_input.csv'
        output_file = 'gm3_output.csv'

        with pytest.raises(Exception):
            get_matches(input_file, output_file)

    
    # def test_when_missing_input_params(self):
    #     """
    #     Explores how get_matches() responds when objects with missing parameters are input
    #     """

    #     input_file = 'gm4_input.csv'
    #     output_file = 'gm4_output.csv'

    #     get_matches(input_file, output_file)


