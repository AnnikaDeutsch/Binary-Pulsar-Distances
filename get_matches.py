import astropy
import astropy.units as u
import astroquery
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import numpy as np

# inputs: ra, dec, height_of_rectangle, width_of_rectangle, radius_of_circle
def psr_to_gaia(jname, raj, rajerr, decj, decjerr,  pmra, pmraerr, pmdec, pmdecerr, posepoch, height, width):
    """Search Gaia for Possible Companion to Pulsar

    Given input parameters read in from a text file, queries Gaia DR2 to find matches 
    based on RA and Dec for each object from the text file

    Args:
        jname (str): Name of the pulsar being checked for matches 
        raj (str): Right ascension of the pulsar in hh:mm:ss.ss format
        decj (str): Declination of the pulsar in degrees:mm:ss.ss format
        radius (str): Radius of ATNF search in degrees
        height (str): Height of Gaia search in mas 
        width (str): Width of Gaia search in mas 
    
    Returns:
        Table: results of the Gaia query in an astropy Table  
    """

    #Import things

    # Given an item (via specifying ra/dec range or optionally inputing a pulsar name) from ATNF catalouge, 
    # find nearby Gaia objects
    from astropy.time import Time 
    p_ra = raj
    p_dec = decj
    p_pmra = pmra
    p_pmdec = pmdec
    p_epoch = Time(posepoch, format='mjd').jyear

    # convert ra and dec to angle objects that will know to behave as floats 
    p_ra_ang = Angle(p_ra, u.degree)
    p_dec_ang = Angle(p_dec, u.degree)

    # convert pmra and pmdec from mas/yr to deg/yr
    p_pmra_deg = (float(p_pmra) * u.mas).to(u.deg) / u.yr
    p_pmdec_deg = (float(p_pmdec) * u.mas).to(u.deg) / u.yr

    # update location of pulsar based on difference from gaia epoch and pmra/pmdec
    gaia_epoch = 2015.5 * u.yr
    year_diff = gaia_epoch - p_epoch.tolist() * u.yr

    # get the new ra and dec for the pulsar by updating to gaia epoch 
    p_new_ra = p_ra_ang + (p_pmra_deg * year_diff)
    p_new_dec = p_dec_ang + (p_pmdec_deg * year_diff)

    ra_ext_pos = raj + rajerr + (pmra + pmraerr)*year_diff
    ra_ext_neg = raj - rajerr + (pmra - pmraerr)*year_diff

    dec_ext_pos = decj + decjerr + (pmdec + pmdecerr)*year_diff
    dec_ext_neg = decj - decjerr + (pmdec - pmdecerr)*year_diff


    # Query Gaia within the range of the given pulsar 
    coord=SkyCoord(ra=p_new_ra, dec=p_new_dec, unit=(u.degree, u.degree), frame='icrs')
    width_gaia = u.Quantity(width, u.mas)
    height_gaia = u.Quantity(height, u.mas)
    results = Gaia.query_object_async(coordinate=coord, width=width_gaia, height=height_gaia)
    
    if len(results) == 0:
        return results
    else:
        results.add_column(jname, name='Companion Pulsar', index=0)
        return results

def get_matches(input_file, output_file, height=1.*u.arcmin, width=1.*u.arcmin):
    """Give Gaia matches to Pulsars 

    Takes as input a text file (.csv file) with index number, name, ra, dec, proper
    motion ra, proper motion dec and posepoch and produces all of the gaia matches
    of ra and dec to within a certain range.

    Args: 
        input_file (str): Name of the text file (csv) containing each pulsar with the parameters 'index', 'name', 
          'ra', 'dec', 'pmra', 'pmdec', 'posepoch' row by row for each object.
        output_file (str): Name of the text file which the pulsar-gaia matches will be output to.
        height (:obj:'float', optional): Height of the rectangle Gaia will query in.
        width (:obj:'float', optional): Width of the rectangle Gaia will query in.

  
    """
    from astropy.table import Table, vstack

    f = open(input_file, "r")
    results = Table()
    first_time = True

    counter = 0

    # Loop through file of ATNF data and combine tables of Gaia matches into one supertable
    for line in f:
      # Loop 20 times
      counter += 1
      if (counter == 20):
        break

      # Parse input
      values = line.split(';')

      # Add result to supertable
      search_result = psr_to_gaia(values[1],values[2],values[3],values[4],values[5],values[6], height, width)
      if (len(search_result) == 0):
        continue
      if first_time:
        results = search_result
        first_time = False
      else:  
        results = vstack([results, search_result])

    results.write(output_file, format='csv', overwrite=True)


