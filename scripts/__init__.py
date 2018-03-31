# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 11:19:00 2017

@author: shooper
"""
import sys, os, warnings

this_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(this_dir)
for name in os.listdir(this_dir):
    path = os.path.join(this_dir, name)
    if os.path.isdir(path):
        sys.path.append(path)

try:
    from evaluation import *
except:
    warnings.warn('Could not import evaluation module', ImportWarning)
    pass

try:
    from utils import *
except:
    warnings.warn('Could not import utils module', ImportWarning)
    pass

try:
    from randomforest import *
except:
    warnings.warn('Could not import randomforest module', ImportWarning)
    pass

        