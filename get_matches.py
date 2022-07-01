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

    # ra_ext_pos = raj + rajerr + (pmra + pmraerr)*year_diff
    # ra_ext_neg = raj - rajerr + (pmra - pmraerr)*year_diff

    # dec_ext_pos = decj + decjerr + (pmdec + pmdecerr)*year_diff
    # dec_ext_neg = decj - decjerr + (pmdec - pmdecerr)*year_diff


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


def compare_pm_param_space(psr_name, pmra, pmra_err, pmdec, pmdec_err, gaia_matches_filename):
  """Plots pm parameter space for Gaia matches of one pulsar in cartesian coordinates 
  
  """
  # instantiate a figure and axis object where we will plot everything relating to the pm parameter space 
  fig, ax = plt.subplots(figsize=(16,4)) 

  # create variables for pmra, pmdec, and their errors
  pmra = pmra * u.mas / u.yr
  pmra_err = pmra_err * u.mas / u.yr
  pmdec = pmdec * u.mas / u.yr
  pmdec_err = pmdec_err * u.mas / u.yr

  # plot the \mu in \alpha and \delta of the actual pulsar 
  with quantity_support(): # to make sure astropy quantities get plotted properly 
      ax.errorbar(pmra, pmdec, pmdec_err, pmra_err, 'bo') 

  # read the csv file containing the gaia "matches" identified by get_matches()
  f = open(gaia_matches_filename, 'r')
  not_zero = 0 # counter that skips the first line in the text file, which is just a header 
  for line in f:
      values = line.split(',') # create a table containing the values of a given line of the file
      if not_zero == 0:
          not_zero+=1 # skips the first line of the file 
      elif values[0] == psr_name: # only looks at matches to the pulsar we are concerned about 
          if values[14] != '': # avoids entries that don't have a pm measurement 
              with quantity_support(): # plot the pm's of a gaia object with error bars 
                  ax.errorbar(float(values[14]), float(values[16]), float(values[15]), float(values[17]), 'ro')
      else:
        break # break once it has gone through every match for the given pulsar; we will have to modify this to work 
              # generically once we put it in a function 

  ax.set_title('Proper motions for possible Gaia matches to PSR {}'.format(psr_name), fontsize=16)
  ax.set_xlabel(r'$\mu_{\alpha}$ (mas)', fontsize=16)
  ax.set_ylabel(r'$\mu_{\delta}$ (mas)', fontsize=16)


def compare_pos_param_space(psr_name, ra, ra_err, dec, dec_err, pmra, pmra_err, pmdec, pmdec_err, pos_epoch, 
  gaia_matches_filename):
  """plots position parameter space for Gaia matches of one pulsar in cartesian coordinates 
  
  """
  fig1, ax1 = plt.subplots(figsize=(16,4)) 

  ra = ra 
  dec = dec 
  ra_err = ra_err*u.deg # must be floats
  dec_err = dec_err*u.deg # must be floats

  ra_ang = Angle(ra, u.deg)
  dec_ang = Angle(dec, u.deg)

  # this part is to plot the region within the propogated error after updating to the gaia epoch

  posepoch = pos_epoch

  p_epoch = Time(posepoch, format='mjd').jyear
  gaia_epoch = 2015.5 * u.yr
  year_diff = gaia_epoch - p_epoch.tolist() * u.yr

  bound1 = (ra_ang + ra_err) + ((pmra.to(u.deg/u.yr) + pmra_err.to(u.deg/u.yr))*year_diff) # right x err
  bound2 = (dec_ang + dec_err) + ((pmdec.to(u.deg/u.yr) + pmdec_err.to(u.deg/u.yr))*year_diff) # top y err

  bound3 = (ra_ang - ra_err) + ((pmra.to(u.deg/u.yr) - pmra_err.to(u.deg/u.yr))*year_diff) # left x err
  bound4 = (dec_ang - dec_err) + ((pmdec.to(u.deg/u.yr) - pmdec_err.to(u.deg/u.yr))*year_diff) # bottom y err
  bounds = [bound1, bound3, bound2, bound4]

  asym_err_x = [[bound3], [bound1]]
  asym_err_y = [[bound4], [bound2]]

  largest_err = 0
  for bound in bounds:
      if bound > largest_err:
          largest_err = bound 

  new_ra = ra_ang + (pmra.to(u.deg/u.yr)*year_diff)
  new_dec = dec_ang + (pmdec.to(u.deg/u.yr)*year_diff)

  region = plt.Circle((new_ra, new_dec), largest_err-new_ra, color='g', alpha=0.2)

  # plot the \mu in \alpha and \delta of the actual pulsar 
  with quantity_support():
      ax1.errorbar(ra_ang, dec_ang, dec_err, ra_err, 'bo') 
      ax1.plot(new_ra, new_dec, 'bo')
      ax1.add_patch(region)
  # plt.xlim([-1,2])
  # plt.ylim([-75,70])
  f = open(gaia_matches_filename, 'r')
  not_zero = 0
  for line in f:
      values = line.split(',')
      if not_zero == 0:
          not_zero+=1
      elif values[0] == 'J0024-7204Z':
          # if float(values[8]) <= 0.1 and float(values[10]) <= 0.1:
              with quantity_support():
                  ax1.errorbar(float(values[7])*u.deg, float(values[9])*u.deg, float(values[8])*u.mas,
                  float(values[10])*u.mas, 'ro')
      else:
          break

  ax1.set_title('Angular offsets of Gaia matches from PSR J0024-7204Z', fontsize=16)
  ax1.set_xlabel(r'$\theta_{\alpha}$', fontsize=16)
  ax1.set_ylabel(r'$\theta_{\delta}$', fontsize=16)

  plt.tight_layout()


