import astropy
import astropy.units as u
import astroquery
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import numpy as np

# inputs: ra, dec, height_of_rectangle, width_of_rectangle, radius_of_circle
def psr_to_gaia(jname, raj, decj, pmra, pmdec, posepoch, height=1000000*u.mas, width=1000000*u.mas, name=""):
    """Search Gaia for Possible Companion to Pulsar

    Given input parameters read in from a text file, queries Gaia DR2 to find matches 
    based on RA and Dec for each object from the text file

    Args:
        jname (str): Name of the pulsar being checked for matches 
        raj (:obj:`str`, optional): Right ascension of the pulsar in hh:mm:ss.ss format
        decj (:obj:`str`, optional): Declination of the pulsar in degrees:mm:ss.ss format
        radius (:obj:`float`, optional): Radius of ATNF search in degrees
        height (:obj:`float`, optional): Height of Gaia search
        width (:obj:`float`, optional): Width of Gaia search
        name (:obj:`str`, optional)): Name of pulsar

    Raises:
        Exception: If name is given but not all three of raj, decj, and radius are given

    # """

    # if name == "":
    #     if (raj == None or decj == None or radius == None):
    #         raise Exception("If no name given, must provide a right ascension, declination, and radius")

    #Import things

    # Given an item (via specifying ra/dec range or optionally inputing a pulsar name) from ATNF catalouge, find nearby
    # Gaia objects
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

    # print('for '+table['JNAME'][pulsar]+' with Right Ascension'+p_new_ra+' and Declination'
    # +p_new_dec+' in the Gaia epoch, the following Gaia objects are at the same RA and Dec to within a'+
    # width+' by '+height+' milliacrsecond range:')

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


from astropy.table import Table, vstack

f = open("name_input.csv", "r")
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
  search_result = psr_to_gaia(values[1],values[2],values[3],values[4],values[5],values[6])
  if (len(search_result) == 0):
    continue
  if first_time:
    results = search_result
    first_time = False
  else:  
    results = vstack([results, search_result])

results.write('matches2.ecsv', format='csv', overwrite=True)