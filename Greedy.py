# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 12:23:09 2021

@author: d1p3y
"""
'''==========Solving Makespan Scheduling Problem (MS) by Gentic Algorithm in python======='''
# importing required modules
import pandas as pd
import numpy as np
import time
import copy
import random
import math

def runfile(n,nume,inst):

    ''' ================= initialization setting ======================'''    
    INF = 10**8
      
    input_data = open("data/Input_data_%s_%s_%s.txt" %(n,nume,inst), "r")
    list_of_lines = input_data.readlines()
    MaxTime = int(list_of_lines[0])
    Num_Mchn = int(list_of_lines[1])
    jobstr = list_of_lines[2]
    jobtemp = jobstr.split()
    Jobsizes = [int(s) for s in jobtemp]
    Num_Jobs = len(Jobsizes)
    input_data.close()    
    
    Indexlist = [i for i in range(Num_Jobs)]
    
    Jobtable = np.array([Indexlist, Jobsizes])
    #print(Jobtable)
    
    #Sort jobs according to size
    Orderedjobtable = (-Jobtable)[:, (-Jobtable)[1, :].argsort()]
    Orderedjobtable = -Orderedjobtable
    #print(Orderedjobtable)     
    
    Allocation = [[] for m in range(Num_Mchn) ] 
    for m in range(Num_Mchn):
        Allocation[m].append(Orderedjobtable[0][m])
        
    Time = [Orderedjobtable[1][m] for m in range(Num_Mchn) ]
        
    for j in range(Num_Mchn, Num_Jobs):
        shortest = INF
        earlist = 0
        for m in range(Num_Mchn):
            if Time[m] < shortest:
                shortest = Time[m]
                earlist = m
        Allocation[earlist].append(Orderedjobtable[0][j])
        Time[earlist] += Orderedjobtable[1][j]
    
    Maxspan = 0
    Maxmachin = 0
    for m in range(Num_Mchn):
        if Time[m] > Maxspan:
            Maxspan = Time[m]
            Maxmachin = m
        
    sequence = [0 for i in range(Num_Jobs)]
    for m in range(Num_Mchn):
        for j in Allocation[m]:
            sequence[j] = m
    
    #print(sequence)
    #print(Maxspan)
    
    return Maxspan, sequence 
