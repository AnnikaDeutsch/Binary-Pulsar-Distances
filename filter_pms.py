# script to compare pmra and pmdec of the identified matches 
# many of the confirmed literature matches do not have great agreement in their proper motions, so 
# I realized narrowing them down that way is not particularly helpful...

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

a = open('all_final.csv', 'r') # read output file of get_mathces() that has pms of Gaia objects
first_time = True
new = []
top = True

for line in a:
    values = line.split(',')
    if top:
        top = False
        continue
    b = open('all_almost_final.csv', 'r') # read input file of get_matches() that has pms of pulsars
    gaia_pmra = float(values[13])
    gaia_pmdec = float(values[15])
    for line in b:
        check = line.split(';')
        if values[0] == check[1]:
            psr_pmra = float(check[4])
            psr_pmdec = float(check[5])
            if int(gaia_pmra) == int(psr_pmra) and int(gaia_pmdec) == int(psr_pmdec):
                values.pop()
                if first_time:
                    new.insert(0,values)
                    first_time = False
                else:
                    new.append(values)


import csv
with open('all_binary_matches.csv', 'w') as h:
    write = csv.writer(h, delimiter=';')
    write.writerows(new)
