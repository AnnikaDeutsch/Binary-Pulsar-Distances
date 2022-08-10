# script to make a csv of just the NAMES of the pulsars in Globular Clusters, from a file that had a lot
# more info

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