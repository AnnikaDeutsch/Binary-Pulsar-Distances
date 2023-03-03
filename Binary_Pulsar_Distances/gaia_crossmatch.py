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

a=[1,2,3]
b=['a','b','c']
table = Table([a,b], names=['col1','col2'], meta={'meta':'first table'})
# Upload
Gaia.login()
Gaia.upload_table(upload_resource=table, table_name='table_test_from_astropy1')
Gaia.update_user_table(table_name="table_test_from_astropy1",
                       list_of_changes=[["col1", "flags", "Ra"], 
                                        ["col2", "flags", "Dec"]])