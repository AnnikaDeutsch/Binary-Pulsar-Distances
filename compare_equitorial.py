# %%
import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
import matplotlib.cm as cmap
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from mpl_toolkits.axisartist import SubplotHost

from mpl_toolkits.axisartist import GridHelperCurveLinear

import astropy.units as u

from astropy.coordinates import Angle

from astropy.time import Time 

from astropy.visualization import quantity_support
quantity_support()

# %%
def get_ra_matches_deg(matchfile, psr, plot=False):
    """Access right ascensions from a table of gaia matches in degrees

    From a table of Gaia matches, access the right ascensions for a given pulsar, either as astropy quantities
    to keep track of units, or if the table is for plotting, as dimensionless floats or ints so numpy accepts
    them.

    Args:
        matchfile (str): Full filename of the csv containing all the gaia matches for a list of pulsars,
            including all of their parameters as passed over from the output of get_matches().
        psr (str): Name of the pulsar being whose right ascensions should be accessed 
        plot (:obj:'bool', optional): Tells the function whether to record the quantities with or without 
            units depending on whether or not the table will be used to plot data

    Returns:
        Array: Array of right ascensions in the same order as the matches from the input file
    """

    f = open(matchfile, 'r')
    ras = []

    if plot:
        for line in f:
            values = line.split(',')
            if values[0] == psr:
                ras.append(float(values[7]))

        return ras
    else:
        for line in f:
            values = line.split(',')
            if values[0] == psr:
                ras.append(Angle(float(values[7]), u.deg))

        return ras  

def get_ra_matches_hms(matchfile, psr, plot=False): 
    """Access right ascensions from a table of gaia matches in hms 

    From a table of Gaia matches, access the right ascensions for a given pulsar, either as astropy quantities
    to keep track of units, or if the table is for plotting, as dimensionless floats or ints so numpy accepts
    them.

    Args:
        matchfile (str): Full filename of the csv containing all the gaia matches for a list of pulsars,
            including all of their parameters as passed over from the output of get_matches().
        psr (str): Name of the pulsar being whose right ascensions should be accessed 
        plot (:obj:'bool', optional): Tells the function whether to record the quantities with or without 
            units depending on whether or not the table will be used to plot data

    Returns:
        Array: Numpy array of right ascensions in the same order as the matches from the input file
    """
    f = open(matchfile, 'r')
    ras = []

    if plot:
        for line in f:
            values = line.split(',')
            if values[0] == psr:
                ras.append(float(values[7]))

        return np.array(ras)
    else:
        for line in f:
            values = line.split(',')
            if values[0] == psr:
                ras.append(Angle(float(values[7]), u.deg).hour)

        return np.array(ras)

# %%
def get_dec_matches_deg(matchfile, psr, plot=False):
    """Access declinations from a table of gaia matches in degrees 

    From a table of Gaia matches, access the declinations for a given pulsar, either as astropy quantities
    to keep track of units, or if the table is for plotting, as dimensionless floats or ints so numpy accepts
    them.

    Args:
        matchfile (str): Full filename of the csv containing all the gaia matches for a list of pulsars,
            including all of their parameters as passed over from the output of get_matches().
        psr (str): Name of the pulsar being whose declinations should be accessed 
        plot (:obj:'bool', optional): Tells the function whether to record the quantities with or without 
            units depending on whether or not the table will be used to plot data

    Returns:
        Array: Array of declinations in the same order as the matches from the input file
    """

    f = open(matchfile, 'r')
    decs = []

    if plot:
        for line in f:
            values = line.split(',')
            if values[0] == psr:
                decs.append(float(values[9]))

        return decs
    else:
        for line in f:
            values = line.split(',')
            if values[0] == psr:
                decs.append(Angle(float(values[9]), u.deg))

        return decs

def get_dec_matches_dms(matchfile, psr, plot=False):
    """Access declinations from a table of gaia matches in dms

    From a table of Gaia matches, access the declinations for a given pulsar, either as astropy quantities
    to keep track of units, or if the table is for plotting, as dimensionless floats or ints so numpy accepts
    them.

    Args:
        matchfile (str): Full filename of the csv containing all the gaia matches for a list of pulsars,
            including all of their parameters as passed over from the output of get_matches().
        psr (str): Name of the pulsar being whose declinations should be accessed 
        plot (:obj:'bool', optional): Tells the function whether to record the quantities with or without 
            units depending on whether or not the table will be used to plot data

    Returns:
        Array: Array of declinations in the same order as the matches from the input file
    """
    f = open(matchfile, 'r')
    decs = []

    if plot:
        for line in f:
            values = line.split(',')
            if values[0] == psr:
                decs.append(Angle(float(values[9]), u.deg).deg)
    #.to_string(unit=u.deg, sep=':')

        return np.array(decs)
    else:
        for line in f:
            values = line.split(',')
            if values[0] == psr:
                decs.append(Angle(float(values[9]), u.deg).dms)

        return np.array(decs)

# %%
def get_pmras(matchfile, psrname, plot=False):
    """Access proper motions in right ascension from a table of gaia matches

    Access a table of the proper motions of the gaia matches to a given pulsar from an output file of matches
    to a list of pulsars. Proper motions can either be given as astropy quantities with units of u.mas/u.yr, or as
    unitless quantities to be accepted by numpy for plotting.

    Args:
        matchfile (str): Full filename of the csv containing all the gaia matches for a list of pulsars,
            including all of their parameters as passed over from the output of get_matches().
        psr (str): Name of the pulsar being whose proper motions should be accessed 
        plot (:obj:'bool', optional): Tells the function whether to record the quantities with or without 
            units depending on whether or not the table will be used to plot data
    
    Returns:
        Array: Numpy array of proper motions in right ascension in the same order as the matches from 
            the input file
    """
    # all pms are in mas/yr
    f = open(matchfile, 'r')
    pmras = []

    if plot:
        for line in f:
            values = line.split(',')
            if values[0] == psrname:
                if values[14] == '':
                    pmras.append(0)
                else:
                    pmras.append(float(values[14]))

        return np.array(pmras)
    else:
        for line in f:
            values = line.split(',')
            if values[0] == psrname:
                if values[14] == '':
                    pmras.append(0*u.mas/u.yr)
                else:
                    pmras.append(float(values[14])*u.mas/u.yr)

        return np.array(pmras)


# %%
def get_pmdecs(matchfile, psrname, plot=False):
    """Access proper motions in declination from a table of gaia matches

    Access a table of the proper motions of the gaia matches to a given pulsar from an output file of matches
    to a list of pulsars. Proper motions can either be given as astropy quantities with units of mas/yr, or as
    unitless quantities to be accepted by numpy for plotting.

    Args:
        matchfile (str): Full filename of the csv containing all the gaia matches for a list of pulsars,
            including all of their parameters as passed over from the output of get_matches().
        psr (str): Name of the pulsar being whose proper motions should be accessed 
        plot (:obj:'bool', optional): Tells the function whether to record the quantities with or without 
            units depending on whether or not the table will be used to plot data
    
    Returns:
        Array: Numpy array of proper motions in declination in the same order as the matches from 
            the input file
    """
    # all pms are in mas/yr
    f = open(matchfile, 'r')
    pmdecs = []

    if plot:
        for line in f:
            values = line.split(',')
            if values[0] == psrname:
                if values[16] == '':
                    pmdecs.append(0)
                else:
                    pmdecs.append(float(values[16]))

        return np.array(pmdecs)
    else: # this part is going to throw an error having to do with numpy not liking quantities with units
        for line in f:
            values = line.split(',')
            if values[0] == psrname:
                if values[16] == '':
                    pmdecs.append(0*u.mas/u.yr)
                else:
                    pmdecs.append(float(values[16])*u.mas/u.yr)

        return np.array(pmdecs)

# %%
def get_distances(matchfile, psrname, plot=False):
    """Access distances from a table of gaia matches

    Access a table of the distances in Mpc of the gaia matches to a given pulsar from an output file of matches
    to a list of pulsars. Distances can either be given as astropy quantities with units of u.Mpc, or as
    unitless quantities to be accepted by numpy for plotting.

    Args:
        matchfile (str): Full filename of the csv containing all the gaia matches for a list of pulsars,
            including all of their parameters as passed over from the output of get_matches().
        psr (str): Name of the pulsar being whose distances should be accessed 
        plot (:obj:'bool', optional): Tells the function whether to record the quantities with or without 
            units depending on whether or not the table will be used to plot data
    
    Returns:
        Array: Array of distances in the same order as the matches from the input file
    """
    f = open(matchfile, 'r')
    dists = []

    if plot:
        for line in f:
            values = line.split(',')
            if values[0] == psrname:
                if values[11] == '':
                    dists.append(0)
                else:
                    dists.append(1/float(values[11])) # 1/parallax gives Mpc in this case as parallax is in mas
        return dists
    else:
        for line in f:
            values = line.split(',')
            if values[0] == psrname:
                if values[11] == '':
                    dists.append(0*u.Mpc)
                else:
                    dists.append(1/float(values[11])*u.Mpc) # 1/parallax gives Mpc in this case as parallax is in mas
        return dists

# %%
def curvelinear_test2(fig):
    """ 
    
    """

    global ax1

    tr = Affine2D().scale(np.pi/180., 1.) + PolarAxes.PolarTransform()

    extreme_finder = angle_helper.ExtremeFinderCycle(10, 60, 
                                                    lon_cycle= 360,
                                                    lat_cycle= None,
                                                    lon_minmax= None,
                                                    lat_minmax= (0, np.inf),
                                                    )
    
    grid_locator1 = angle_helper.LocatorHMS(12)
    tick_formatter1 = angle_helper.FormatterHMS()

    grid_locator2 = angle_helper.LocatorDMS(6)
    tick_formatter2 = angle_helper.FormatterDMS()

    grid_helper = GridHelperCurveLinear(tr,
                                        extreme_finder= extreme_finder,
                                        grid_locator1= grid_locator1,
                                        tick_formatter1= tick_formatter1,
                                        grid_locator2= grid_locator2,
                                        tick_formatter2= tick_formatter2)

    ax1 = SubplotHost(fig, 1, 1, 1, grid_helper= grid_helper)

    ax1.axis['right'].major_ticklabels.set_visible(True)
    ax1.axis['top'].major_ticklabels.set_visible(True)
    ax1.axis['bottom'].major_ticklabels.set_visible(True)

    ax1.axis['right'].get_helper().nth_coord_ticks=0
    ax1.axis['bottom'].get_helper().nth_coord_ticks=0

    fig.add_subplot(ax1)

    grid_helper = ax1.get_grid_helper()

    ax1.set_aspect(1.)
    ax1.set_xlim(-1,1)
    ax1.set_ylim(-1,1)

    ax1.set_xlabel('Right Ascension')
    ax1.set_ylabel('90$^\circ$ + Declination')
    ax1.grid(True)

    return tr

# %%
import matplotlib.pyplot as plt

fig = plt.figure(1,figsize= (5,5), constrained_layout=True)
plt.rcParams['axes.titley'] = 1.0
plt.rcParams['axes.titlepad'] = 25
fig.clf()

tr = curvelinear_test2(fig) # tr


decs = get_dec_matches_dms('testing_output.csv', 'J0024-7204Z', plot=True)
ras = get_ra_matches_hms('testing_output.csv', 'J0024-7204Z', plot=True)


# out_test = tr.transform(zip(ras, decs))
both = np.vstack((decs, ras))
out_test = tr.transform(both)

dist = get_distances('testing_output.csv', 'J0024-7204Z', plot=True)
print(dist)

cm = plt.cm.get_cmap('RdYlBu_r')
z = dist # for z to actually represent distances, write a function that gets
# distances for the matches in an array and set z to that, rather than random

with quantity_support():
    SC = ax1.scatter(out_test[:,0],
        out_test[:,1],
        c = z,
        cmap = cm,
        zorder = 9)

cbar = plt.colorbar(SC, shrink = 1., pad=0.2)
cbar.ax.tick_params(labelsize=8)
cbar.set_label('distance (Mpc)', fontsize=8)

plt.title('J0024-7204Z')
plt.show
plt.savefig('J0024-7204Z_equitorial_pos.pdf')

# %%
def plot_pos_eqt(inputfile, psrname, pmra=False, pmdec=False, dist=False):
    """Plots the positions of matches in polar coordinates.

    Plots in polar coordinates the positions of the gaia matches to a given pulsar. Specifying one of the 
    optional parameters to be true will add a colorbar of that parameter to the side of the graph, to show
    a thrid dimension of the parameters. One of the three optional boolean parameters must be set True for the 
    function to plot anything.

    Args:
        inputfile (str): Full filename of the csv file holding all the gaia matches that will be plotted. 
            Equivalent to the output file return from get_matches().
        psrname (str): Name of the pulsar for which we are plotting the match positions for comparison.
        pmra (:obj:'bool', optional): Also include a colorbar of the pmras for the given pulsar. Default 
            value is False.
        pmdec (:obj:'bool', optional): Also include a colorbar of the pmdecs for the given pulsar. Default 
            value is False.
        dist (:obj:'bool', optional): Also include a colorbar of the distances for the given pulsar. Default 
            value is False.
    """
    fig = plt.figure(1,figsize= (5,5), constrained_layout=True)
    plt.rcParams['axes.titley'] = 1.0
    plt.rcParams['axes.titlepad'] = 25
    fig.clf()

    tr = curvelinear_test2(fig) # tr

    decs = get_dec_matches_dms(inputfile, psrname, plot=True)
    ras = get_ra_matches_hms(inputfile, psrname, plot=True)

    print('raw ras data')
    print(ras)
    print('raw decs data')
    print(decs)
    

    # out_test = tr.transform(zip(ras, decs))
    both = np.vstack((decs, ras))
    out_test = tr.transform(both)

    print('values of out_test')
    print(out_test)
    print('ras that are being plotted:')
    print(out_test[:,0])
    print('decs that are being plotted:')
    print(out_test[:,1])

    if pmra:
        pmras = get_pmras(inputfile, psrname, plot=True)

        cm = plt.cm.get_cmap('RdYlBu_r')
        z = pmras
        name = 'pmra (mas/yr)'
        file = '_pmra'

    if pmdec:
        pmdecs = get_pmdecs(inputfile, psrname, plot=True)

        cm = plt.cm.get_cmap('RdYlBu_r')
        z = pmdecs
        name = '$\mu_{\delta}$ (mas/yr)'
        file = '_pmdec'

    if dist:
        dists = get_distances(inputfile, psrname, plot=True)

        cm = plt.cm.get_cmap('RdYlBu_r')
        z = dists
        name = 'distance (Mpc)'
        file = '_dist'

    with quantity_support():
        SC = ax1.scatter(ras, 
            decs,
            c = z,
            cmap = cm,
            zorder = 9)

            #out_test[:,0], out_test[:,1],

    cbar = plt.colorbar(SC, shrink = 1., pad=0.2)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(name, fontsize=8)

    plt.title(psrname)
    plt.show
    plt.savefig(psrname + '_eqt_pos' + file + '.pdf')

# %%
plot_pos_eqt('cross_check.csv', 'J1012+5307', pmra=True)

# %%
get_pmdecs('actual_match.csv', 'J1816+4510', plot=True)

plot_pos_eqt('actual_match.csv', 'J1816+4510', pmdec=True)
# plot_pos_eqt('actual_match.csv', 'J1816+4510', pmdec=True)
# plot_pos_eqt('actual_match.csv', 'J1816+4510', dist=True)

# %%
plot_pos_eqt('testing_output.csv', 'J0024-7204Z', pmra=True)

# %%
def compare_pos_param_space(psr_name, ra, ra_err, dec, dec_err, pmra, pmra_err, pmdec, pmdec_err, pos_epoch, 
                            gaia_matches_filename):
    """Plots position parameter space for gaia matches of one pulsar in cartesian coordinates

    Plots the right ascension vs declination for all the matches of a given pulsar in cartesian coordinates, 
    just to see the comparison between the position of the pulsar and the identified matches in the parameter 
    space

    Args:
        psr_name (str): Name of the pulsar whose positions we are plotting and comparing to the position of 
            the pulsar.
        ra (str): Right ascension of the pulsar in hms exactly as it is in the input file of pulsars 
            from ATNF, as a string. (hms string separated by :)
        ra_err (float): Right ascension error of the pulsar as a float. (deg)
        dec (str): Declinations of the pulsar in hms exactlt as it is in the input file of pulsars 
            from ATNF, as a string. (hms string separated by :)
        dec_err (float): Declination error of the pulsar as a float. (deg)
        pmra (float): Proper motion in ra of the pulsar as a float (mas/yr)
        pmra_err (float): Proper motion in ra error of the pulsar as a float (mas/yr)
        pmdec (float): Proper motion in dec of the pulsar as a float (mas/yr)
        pmdec_err (float): Proper motion in dec of the pulsar as a float (mas/yr)
        pos_epoch (int): Epoch that the pulsar data was taken as an int (MJD)
        gaia_matches_filename (str): Full filename of the csv file that holds the gaia matches for a list of 
            pulsars, equivalent to the output of get_matches. This is the file that will be parsed to plot the
            positions of the matches 
    """
    fig1, ax1 = plt.subplots(figsize=(4,1)) 

    ra = ra 
    dec = dec 
    ra_err = ra_err*u.deg # must be floats
    dec_err = dec_err*u.deg # must be floats

    pmra = pmra*u.mas/u.yr
    pmra_err = pmra_err * u.mas / u.yr
    pmdec = pmdec * u.mas / u.yr
    pmdec_err = pmdec_err * u.mas / u.yr

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
    temp_zero = 0
    for line in f:
        values = line.split(',')
        if temp_zero == 0:
            temp_zero+=1
        elif values[0] == psr_name:
            # if float(values[8]) <= 0.1 and float(values[10]) <= 0.1:
                with quantity_support():
                    ax1.errorbar(float(values[7])*u.deg, float(values[9])*u.deg, float(values[8])*u.mas,
                    float(values[10])*u.mas, 'ro')
        else:
            break

    ax1.set_title('Angular offsets of Gaia matches from PSR' + psr_name, fontsize=16)
    ax1.set_xlabel(r'$\theta_{\alpha}$', fontsize=16)
    ax1.set_ylabel(r'$\theta_{\delta}$', fontsize=16)

    plt.xlim([ra_ang - 0.01*u.deg, ra_ang + 0.01*u.deg])
    plt.ylim([dec_ang - 0.01*u.deg, dec_ang + 0.01*u.deg])

    plt.tight_layout()

    plt.savefig(psr_name + 'pos_cartesian.pdf')

# %%
compare_pos_param_space('J1816+4510', '18:16:35.93436', 0.00007, '45:10:33.8618', 0.0004, 5.3, 0.8, -3., 1., 56047, 'actual_match.csv')

# %%
def compare_pos_param_space2(inputfile, psrname, pmra=False, pmdec=False, dist=False): 
    """Plots position parameter space for gaia matches of one pulsar in a different way

    A second (and hopefully easier way) of plotting the right ascension versus declination of the gaia matches 
    to a given pulsar. This function can also allow the option of including a colorbar to show a third 
    parameter.

    Args:
        inputfile (str): Full filename of the csv file that holds all of the gaia matches and their 
            parameters. Equivalent to the output of get_matches()
        psrname (str): Name of the pulsar whose matches right ascensions vs. declinations are being plotted. 
        pmra (:obj:'bool', optional): Also include a colorbar of the pmras for the given pulsar. Default 
            value is False.
        pmdec (:obj:'bool', optional): Also include a colorbar of the pmdecs for the given pulsar. Default 
            value is False.
        dist (:obj:'bool', optional): Also include a colorbar of the distances for the given pulsar. Default 
            value is False.
    
    """
    fig = plt.figure(1,figsize= (3,3), constrained_layout=True)
    plt.rcParams['axes.titley'] = 1.0
    plt.rcParams['axes.titlepad'] = 25
    fig.clf()


    decs = get_dec_matches_dms(inputfile, psrname, plot=True)
    ras = get_ra_matches_hms(inputfile, psrname, plot=True)

    print(decs)
    print(ras)

    if pmra:
        pmras = get_pmras(inputfile, psrname, plot=True)

        cm = plt.cm.get_cmap('RdYlBu_r')
        z = pmras
        name = 'pmra (mas/yr)'
        file = '_pmra'

    if pmdec:
        pmdecs = get_pmdecs(inputfile, psrname, plot=True)

        cm = plt.cm.get_cmap('RdYlBu_r')
        z = pmdecs
        name = '$\mu_{\delta}$ (mas/yr)'
        file = '_pmdec'

    if dist:
        dists = get_distances(inputfile, psrname, plot=True)

        cm = plt.cm.get_cmap('RdYlBu_r')
        z = dists
        name = 'distance (Mpc)'
        file = '_dist'

    with quantity_support():
        plt.scatter(ras,
            decs,
            c = z,
            cmap = cm,
            zorder = 9)

    plt.xlim([np.amin(ras)-.01,np.amax(ras)+.01])
    plt.ylim([np.amin(decs)-.01,np.amax(decs)+.01])

    cbar = plt.colorbar(SC, shrink = 1., pad=0.2)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(name, fontsize=8)

    plt.title(psrname)
    plt.show
    plt.savefig(psrname + '_cartesian_pos' + file + '.pdf')

# %%
compare_pos_param_space2('actual_match.csv', 'J1816+4510', pmra=True)

# %%
compare_pos_param_space2('cross_check.csv', 'J1012+5307', pmra=True)

# %%
# get ra and dec arrays
ras = get_ra_matches_hms('actual_match.csv','J1816+4510', plot=True)
decs = get_dec_matches_dms('actual_match.csv','J1816+4510', plot=True)

print(ras)
print(decs)

pmras = get_pmras('actual_match.csv', 'J1816+4510', plot=True)

# convert ra from decimal degrees to radians -- will I want ras and decs to be in degrees?
ras = [x / 180.0 * np.pi for x in ras]

# actually make the plot
fig = plt.figure(figsize=(5,5))
gs = gridspec.GridSpec(4,2) # a grid layout to place subplots in a figure; I don't rlly understand it

# position plot in the figure using gridspec
ax = plt.subplot(gs[0], polar=True)
ax.set_ylim([90,30]) # this is dependent on the range of decs for the given psr matches

# set x, y ticks
angs = np.array([330., 345., 0., 15., 30., 45., 60.]) # angles in degrees 
plt.xticks(angs * np.pi / 180., fontsize= 8)
plt.yticks(np.arange(90,30,10), fontsize= 8)
ax.set_rlabel_position(120)
ax.set_xticklabels(['$22^h$', '$23^h$','$0^h$','$1^h$','$2^h$','$3^h$','$4^h$'], fontsize= 10)
ax.set_yticklabels(['$80^{\circ}$', '$70^{\circ}$', '$60^{\circ}$', '$50^{\circ}$', '$40^{\circ}$'])

cm = plt.cm.get_cmap('RdYlBu_r')

# Plot the points 
ax.scatter(ras, decs, c=pmras)

# Colorbar
cbar = plt.colorbar(SC, shrink=1., pad=0.05)
cbar.ax.tick_params(labelsize=8)
cbar.set_label('pmra (mas)', fontsize=8)







