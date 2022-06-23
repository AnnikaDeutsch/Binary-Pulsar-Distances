import astropy
import astropy.units as u
# inputs: ra, dec, height_of_rectangle, width_of_rectangle, radius_of_circle
def psr_to_gaia(raj, decj, pmra, pmdec, posepoch, height=10*u.mas, width=10*u.mas, name=""):
    """Search Gaia for Possible Companion to Pulsar

    Can search for a pulsar with the given name, or if no name is given, will search for pulsar
    with the provided raj, decj, and radius.

    Args:
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
    import astropy
    import astroquery
    from astroquery.gaia import Gaia
    from astropy.coordinates import SkyCoord
    from astropy.coordinates import Angle
    import astropy.units as u
    import numpy as np

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
    results.show_in_notebook()


f = open("input.csv", "r")
for line in f:
  values = line.split(';')
  psr_to_gaia(values[1], values[2], values[3], values[4], values[5])