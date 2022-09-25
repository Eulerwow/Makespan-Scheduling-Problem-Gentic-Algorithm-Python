# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 10:32:34 2021

@author: d1p3y
"""
'''======= Main script to call greedy algorithm and GA ========='''
import csv

import Greedy as gd
import GAalgo as ga
import time

def call_and_write(x1,x2):
    for x in range(x1, x2+1):
        n = 100*x
        
        for nume in range(1, 10):
            m = round(n*nume/10)
            path = "result_data/Result_data_%d_%d.csv "%(n, nume)
            with open(path, 'w+') as result_data:
                csv.writer(result_data, delimiter='\t').writerow(['Input Data','Greedy Result', 'GA Result', 'Greedy Time', 'GA Time'])
            for instance in range(1, 11):
                print('running greedy on data n = %d, m = %d, instance %d' %(n,m,instance))
                start = time.time()            
                greedy_result, greedy_schedule = gd.runfile(n,nume,instance)
                end = time.time()
                greedy_time = end - start
                print('makespan from greedy algorithm is: %d' %greedy_result)
                
                print('running GA on data n = %d, m = %d, instance %d' %(n,m,instance))                
                start = time.time()
                GA_result, GA_schedule = ga.runfile(n,nume,instance)
                end = time.time()
                GA_time = end - start
                print('makespan from GA is: %d' %GA_result)
                
                input_data = str(n)+"-"+str(nume)+"-"+str(instance)
                with open(path, 'a+') as result_data:
                    csv.writer(result_data, delimiter='\t').writerow([input_data, greedy_result, GA_result, greedy_time, GA_time])
                print('finished running on data n = %d, m = %d, instance %d' %(n,m,instance))
 
                
def call_without_write(n,nume,instance):
    m = round(n*nume/10)
    print('running greedy on data n = %d, m = %d, instance %d' %(n,m,instance))
    m = round(n*nume/10)
    start = time.time()            
    greedy_result, greedy_schedule = gd.runfile(n,nume,instance)
    end = time.time()
    greedy_time = end - start
    print('makespan from greedy algorithm is: %d' %greedy_result)
    print('schedule from greedy algorithm is ', greedy_schedule)
    
    print('running GA on data n = %d, m = %d, instance %d' %(n,m,instance)) 
    start = time.time()
    GA_result, GA_schedule = ga.runfile(n,nume,instance)
    end = time.time()
    GA_time = end - start
    
    print('makespan from GA is: %d' %GA_result)
    print('schedule from GA is ', GA_schedule)
    print("n    m    intance ID    Greedy Result    GA Result   Greedy Time   GA Time")
    print("%d   %d     %02d          %d           %d          %.5f       %.2f" %(n,m,instance,greedy_result,GA_result,greedy_time,GA_time))


'''call both algorithm on data file: n-nume-instance'''    
call_without_write(300,2,1)

#call_and_write(1, 1)
    

