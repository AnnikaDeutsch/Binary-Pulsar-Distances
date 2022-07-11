import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
import matplotlib.cm as cmap
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D
import matplotlib.pyplot as plt

from mpl_toolkits.axisartist import SubplotHost

from mpl_toolkits.axisartist import GridHelperCurveLinear

import astropy.units as u

from astropy.coordinates import Angle

from astropy.visualization import quantity_support
quantity_support()

def get_ra_matches_deg(matchfile, psr, plot=False):
    """
    
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
    f = open(matchfile, 'r')
    ras = []

    if plot:
        for line in f:
            values = line.split(',')
            if values[0] == psr:
                ras.append(Angle(float(values[7]), u.deg).hour)

        return np.array(ras)
    else:
        for line in f:
            values = line.split(',')
            if values[0] == psr:
                ras.append(Angle(float(values[7]), u.deg).hour)

        return np.array(ras)


def get_dec_matches_deg(matchfile, psr, plot=False):
    """
    
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


def get_pmras(matchfile, psrname, plot=False):
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


def get_pmdecs(matchfile, psrname, plot=False):
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
    else:
        for line in f:
            values = line.split(',')
            if values[0] == psrname:
                if values[16] == '':
                    pmdecs.append(0*u.mas/u.yr)
                else:
                    pmdecs.append(float(values[16])*u.mas/u.yr)

        return np.array(pmdecs)


def get_distances(matchfile, psrname, plot=False):
    """
    
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


def plot_pos_eqt(inputfile, psrname, pmra=False, pmdec=False, dist=False):
    """
    
    """
    fig = plt.figure(1,figsize= (5,5), constrained_layout=True)
    plt.rcParams['axes.titley'] = 1.0
    plt.rcParams['axes.titlepad'] = 25
    fig.clf()

    tr = curvelinear_test2(fig) # tr

    # this is the point at which I need to access the actual data. I am realizing that there is a standard way I go about
    # this for each function and so I may as well write them into functions themselves 

    decs = get_dec_matches_dms(inputfile, psrname, plot=True)
    ras = get_ra_matches_hms(inputfile, psrname, plot=True)


    # out_test = tr.transform(zip(ras, decs))
    both = np.vstack((decs, ras))
    out_test = tr.transform(both)

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
        SC = ax1.scatter(out_test[:,0],
            out_test[:,1],
            c = z,
            cmap = cm,
            zorder = 9)

    cbar = plt.colorbar(SC, shrink = 1., pad=0.2)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(name, fontsize=8)

    plt.title(psrname)
    plt.show
    plt.savefig(psrname + '_eqt_pos' + file + '.pdf')


plot_pos_eqt('cross_check.csv', 'J1012+5307', pmra=True)




