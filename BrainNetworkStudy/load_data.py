# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 23:42:01 2018

@author: Mattia
"""

import pyedflib
import mylib 
#import pandas as pd

f = pyedflib.EdfReader('../data/S001R01.edf')

# mylib.print_some_stuff(f)


data = mylib.edfToDataFrame(f)




f._close()
