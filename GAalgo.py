# -*- coding: utf-8 -*-

"""
Created on Sun Oct 10 11:35:30 2021

@author: Pengyuan Ding
"""

'''==========Solving Makespan Scheduling Problem (MS) by Gentic Algorithm in python======='''
# importing required modules

import pandas as pd
import numpy as np
import time
import copy
import random
import math

import Greedy as gd

def runfile(n,nume,inst):
    
    ''' ================= initialization setting ======================'''
    print_each_generation = False #option to print best makespan of each generation
    use_greedy = False #option to put result from greedy algo into initial population
    
    #read data from txt file
    input_data = open("data/Input_data_%s_%s_%s.txt" %(n,nume,inst), "r")
    list_of_lines = input_data.readlines()
    MaxTime = int(list_of_lines[0])
    Num_Mchn = int(list_of_lines[1])
    
    jobstr = list_of_lines[2]
    jobtemp = jobstr.split()
    
    Jobsizes = [int(s) for s in jobtemp]
    total_processing_time = sum(Jobsizes)
    fantasy_makespan = int(total_processing_time/Num_Mchn)
    max_size = max(Jobsizes)
    Num_Jobs = len(Jobsizes)
    joblist = list(range(Num_Jobs))
    input_data.close()
    
    #set GA algorithm parameters
    
    INF = 10**8 #set a large number
    eps = 10**(-5) #set a small number for convergency check
    
    #set population size, larger the jobs quantity gets, smaller the population size need to be
    # population size = 150 for no. of jobs = 1000; size=50 for  no. of jobs = 100
    if Num_Jobs < 100:
        population_size = 50
    elif Num_Jobs <= 1000:
        population_size = int(50 + 2*int(35*(Num_Jobs-100)/900))
    else:
        population_size = 120
    crossover_rate = 0.95 #probablity of a pair of chromosomes crossover
    mutation_rate = 0.65 #probablity of a chromosome mutates
    
    num_mutation_jobs=2
    elite_rate = 0.05
    num_generation = 10000 #max generation of GA algorithm
    converge_gen = int(num_generation*0.4)
    
    cutpoint1 = int(1/3*Num_Jobs)
    cutpoint2 = int(2/3*Num_Jobs)
    
    '''=== functions that the GA would call ==='''
    
    #function that randomly generate a chromosome
    def random_chromosome(Num_Mchn, Num_Jobs):
        chromosome = [random.choice(range(Num_Mchn)) for i in range(Num_Jobs)]
        return chromosome
       
    def crossover(parent1, parent2):
        piece1_1 = parent1[0:cutpoint1]
        piece1_2 = parent1[cutpoint1:cutpoint2]
        piece1_3 = parent1[cutpoint2:]
        piece2_1 = parent2[0:cutpoint1]
        piece2_2 = parent2[cutpoint1:cutpoint2]
        piece2_3 = parent2[cutpoint2:]
        child1 = copy.deepcopy(piece1_1)
        child1.extend(piece2_2)
        child1.extend(piece1_3)
        child2 = copy.deepcopy(piece2_1)
        child2.extend(piece1_2)
        child2.extend(piece2_3)
        return child1, child2  
        
    def mutation(chromosome):
        mutated = copy.deepcopy(chromosome)
        m_chg = random.sample(joblist, num_mutation_jobs)
        while mutated == chromosome:
            for gene in m_chg:
                mutated[gene] = random.choice(range(Num_Mchn))
        return mutated
    
    #calculate the fitness value of a chromosome i.e. alpha*exp(-beta*makespan)
    def fitness(chromosome):
        machinespans = [0 for m in range(Num_Mchn)]
        for j in range(Num_Jobs):
            m = chromosome[j]
            machinespans[m] += Jobsizes[j]
        makespan = max(machinespans)           
        return makespan, fantasy_makespan /(makespan - fantasy_makespan)          
    
        
    '''==================== main code ==============================='''
    start_time = time.time()
    '''----- generate initial population -----'''
    Tbest = INF
    best_list,best_obj=[],[]
    population_list=[]
    makespan_record=[]
    
    gd_result, grdsqn = gd.runfile(n,nume,inst)
        
    #initial population    
    if use_greedy: 
        gdpop = round(population_size*0.5)
        #print(grdsqn)
        for i in range(gdpop):
            improved = copy.deepcopy(grdsqn)
            population_list.append(improved)
            improved = mutation(improved)
        randmgenert = population_size-gdpop
    else:
        randmgenert = population_size
    for i in range(randmgenert):
        population_list.append(random_chromosome(Num_Mchn, Num_Jobs))
        
    #start GA
    for n in range(num_generation):
        Tbest_now = INF           
        
        '''-------- crossover --------'''
        parent_list=copy.deepcopy(population_list)
        offspring_list=copy.deepcopy(population_list)
        S=list(np.random.permutation(population_size)) # generate a random sequence to select the parent chromosome to crossover
        
        for m in range(int(population_size/2)):
            crossover_prob=np.random.rand()
            if crossover_rate>=crossover_prob:
                parent_1= population_list[S[2*m]]
                parent_2= population_list[S[2*m+1]]
                children = crossover(parent_1, parent_2)
                offspring_list[S[2*m]]=children[0]
                offspring_list[S[2*m+1]]=children[1]
                
        '''--------mutation--------'''   
        for m in range(len(offspring_list)):
            mutation_prob=np.random.rand()
            if (mutation_rate >= mutation_prob and num_mutation_jobs > 1):
                offspring_list[m] = mutation(offspring_list[m])
                
        '''--------fitness value(calculate makespan)-------------'''
        total_chromosome=copy.deepcopy(parent_list)+copy.deepcopy(offspring_list)
        chrom_makespan = []
        chrom_fitness = []
        total_fitness = 0
        for m in range(population_size*2):
            chrom_makespan.append(fitness(total_chromosome[m])[0])
            chrom_fitness.append(fitness(total_chromosome[m])[1])
            total_fitness=total_fitness+chrom_fitness[m]
               
        '''----------selection(roulette wheel approach)----------'''
        pk,qk=[],[]
        
        for i in range(population_size*2):
            pk.append(chrom_fitness[i]/total_fitness)
        for i in range(population_size*2):
            cumulative=0
            for j in range(0,i+1):
                cumulative=cumulative+pk[j]
            qk.append(cumulative)
        
        selection_rand=[np.random.rand() for i in range(population_size)]
        
        for i in range(population_size):
            if selection_rand[i]<=qk[0]:
                population_list[i]=copy.deepcopy(total_chromosome[0])
            else:
                for j in range(0,population_size*2-1):
                    if selection_rand[i]>qk[j] and selection_rand[i]<=qk[j+1]:
                        population_list[i]=copy.deepcopy(total_chromosome[j+1])
                        break
                    
        '''----------comparison----------'''
        for i in range(population_size*2):
            if chrom_makespan[i] < Tbest_now:
                Tbest_now = chrom_makespan[i]
                sequence_now = copy.deepcopy(total_chromosome[i])
        if Tbest_now <= Tbest:
            Tbest = Tbest_now
            sequence_best = copy.deepcopy(sequence_now)
            
        makespan_record.append(Tbest)
        
        '''----------choose elite----------'''
        for i in range(round(elite_rate*population_size)):
            population_list[i*2]=copy.deepcopy(sequence_best)

        if print_each_generation:    
            print('generation %s, best makespan %s' %(n, Tbest))
        
        '''---------- stop criteria ----------'''
        if n > converge_gen:
            if abs(makespan_record[n] - makespan_record[n-converge_gen]) < eps:
                print("Stagnant Improvement")
                break
        
        cumutime = time.time() - start_time
        if cumutime >= MaxTime:
            print("Maximum allowable processing time reached")
            break
        
        if n == num_generation - 1:
            print("Maximum generations reached")
        
    '''----------result----------'''
    #print results    
    print("optimal value:%f"%Tbest)
    print('the elapsed time:%s'% (time.time() - start_time))
    
    #plot results
    import matplotlib.pyplot as plt
    #matplotlib inline
    plt.plot([i for i in range(len(makespan_record))],makespan_record,'b')
    plt.ylabel('makespan',fontsize=15)
    plt.xlabel('generation',fontsize=15)
    plt.show()    
    
    #print('fantasy_aid is: ', fantasy_makespan)
    
    return Tbest, sequence_best

        