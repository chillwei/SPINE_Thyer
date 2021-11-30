#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 21:39:56 2021

@author: weiqiyao
"""

test = "AGCGGGAGACCGGGGTCTCTGAGCtaa"
#test_num = len(test)/3;
for i in range(3,len(test)-3,3):
    if i == len(test) - 3:
        test_v2 = test[0:i]
    else:
        test_L = test[0:i]
        test_R = test[(i+3):]
        test_v2 = test_L + test_R
    
    #table = test.maketrans(x, y)
    print(test_v2)