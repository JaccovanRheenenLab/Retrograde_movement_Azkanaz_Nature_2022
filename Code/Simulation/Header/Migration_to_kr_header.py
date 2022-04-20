#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 21:20:43 2020

@author: bernat
"""
import random as ran
import scipy.integrate as integrate
import numpy as np

##################################################################
############### CLONE SURVIVAL PROBABILITY FUNCTIONS #############
##################################################################

def rho(x,k,t,krkd):
    E=-1/(2*krkd)*(((x-k*np.exp(t))**2)/(np.exp(2*t)-1))
    R=np.sqrt(1/((np.exp(2*t)-1))*2*krkd*np.pi)
    T=np.exp(t)
    
    return 0.5*R*T*np.exp(E)

def N_L(k,t,krkd,L):
    integral,res=integrate.quad(rho, -10,L, args=(k,t,krkd))
    
    return integral

def N_infinity(k,t,krkd):
    integral,res=integrate.quad(rho, 0,np.inf, args=(k,t,krkd))
    
    return integral
    
def R(k,t,krkd, L):
    return N_L(k,t,krkd,L)/N_infinity(k,t,krkd)

def psurvival(k,t,krkd, L):
    return 1.0-np.exp(-R(k,t,krkd,L))

def Rel_psurvival(k,t,krkd, L):
    return psurvival(k,t,krkd, L)/psurvival(0,t,krkd, L)
    
def Rel_psurvival_approx(k,t,krkd, L):
    return N_L(k,t,krkd,L)

def Rel_psurvival_Norm(k,t,krkd, L):
    Zdoub=0
    Zint=0
    for i in range(0,L):
        Zdoub=Zdoub+psurvival(i,t,krkd, L)
    Zint=int(Zdoub)

    return (Zint/Zdoub)*psurvival(k,t,krkd, L)
    
def Ns_surv_time(t,krkd, L):
    Zdoub=0
    for i in range(0,L):
        Zdoub=Zdoub+psurvival(i,t,krkd, L)
    if Zdoub >1:
        return Zdoub
    else:
        return 1

def p(x,k,t,krkd,L):
    integral,res=integrate.quad(rho, x-0.5,x+0.5, args=(k,t,krkd))
    Z=0
    for i in range(0,L):
        integralZ,res=integrate.quad(rho, x-0.5,x+0.5, args=(i,t,krkd))
        Z=Z+integralZ

    return integral/Z

def p_gauss_infinity(k,krkd):
    p=0
    p=np.exp(-1/(2*krkd)*k*k)

    return p

def p_fermi_infinity(k,krkd,L):
    E=0
    E=np.exp(-1/(2*krkd)*k*k)
    Q=np.sqrt(1/(2*krkd*np.pi))
    p=E*Q
    p1=np.exp(-p*L)

    return (1.0-p1)

def p_Taylor2nd_infinity(k,krkd,L):
    E=0
    E=np.exp(-1/(2*krkd)*k*k)
    Q=np.sqrt(1/(2*krkd*np.pi))
    p=E*Q

    return p-0.5*(p*L)*(p)+(1/6)*(p*L)*(p*L)*(p)
    
def Gauss(x,a,b):
    x0=0
    y=a*np.exp(-((x0-x)**2)/(2*b))

    return y

def Logistic(t, kd, cells):
    a_0=1/(1+np.exp(-(kd*(0)-np.log2(cells/2)+1)))## initial condition t=0
    result=1/(1-a_0)*(1/(1+np.exp(-(kd*(t-1)-np.log2(cells/2)+1)))-a_0)

    return result

####################################################
##### initializing a system with n cells labelled 0,1,...,n at random positions

def create_system_1D(cells, N):
    system=[-1]*N
    T=1
    cell=0
    n=0
    while T==1:
        r=ran.randint(0,N-1)
        if system[r]==-1:
            system[r]=cell
            cell=cell+1
            n=n+1
            if n==cells:
                T=0

    return(system)

def create_system_1D_Top(cells, N):
    system=[-1]*N
    system[cells]=0
    
    return(system)   

def create_system_1D_confluent(N):
    system=[]
    for i in range(0,N):
        system.append(i)
        
    return(system)   

def create_system_1D_ablation(N,cell):
    system=[]
    for i in range(0,cell):
        system.append(-1)

    for i in range(cell, N):
        system.append(i-cell)

    return(system)   

##### END initializing a system with n cells labelled 0,1,...,n at random positions
####################################################

######### Identifying the deepest cell of the system
######### Empty system = -1-1-1-1-1-1-1...
######### Full system = 0000000000...   
    
def first_cell(system):
    a=-1
    for i in range(0, len(system)):
        if system[i]==0:
            a=i
            break

    return a
    
######### END Identifying the deepest cell of the system
####################################################

####################################################
#### event of migration/swap    

def migration_swap(system, q): ### competition.Prob of swap if it is possible

    r=ran.randint(0,len(system)-1) #### chose a cell to perform a migration step towards position 0  
    hit_temp=system[r]
    if system[r]>-1 and r>0:
        prev_temp=system[r-1]
        if prev_temp==-1: #### free migration
            system[r-1]=hit_temp
            system[r]=prev_temp

        if prev_temp>-1:
            if ran.uniform(0, 1)<q: #### competition with the cell at the bottom
                system[r-1]=hit_temp
                system[r]=prev_temp

    return(system)
    
#### END event of migration/swap    
####################################################

####################################################
##### event of Duplication

def duplication_event(system):
    
    N=len(system)
    system_dup=[]

    for i in range(0, N):
        system_dup.append(system[i])
    
    r=ran.randint(0,N-1) #### chose a cell to duplicate
    if system[r]>-1: #### there is a cell in the position
        if r>0: #not at the bottom
            ##### cases of trivial movement
            if system[r-1]==-1: #### goes down by default if it can 
                system_dup[r-1]=system[r]
            if r<N-1: 
                if system[r-1]>-1 and system[r+1]==-1:
                    system_dup[r+1]=system[r]
                    
            ### neighbours in both sides
            
                if system[r-1]>-1 and system[r+1]>-1: 
                    up=0
                    down=0
                    cut_up=-1
                    cut_down=-1
                    
                    ### computing the number of successive cells without empty spaces
                    for c in range(r+1,N): ##### up cells
                        if system[c]>-1:
                            up=up+1
                        if system[c]==-1:
                            cut_up=c
                            break
                        
                    for c in range(0,r): ##### down cells
                        if system[r-1-c]>-1:
                            down=down+1
                        if system[r-1-c]==-1:
                            cut_down=r-1-c
                            break
                    
                    if cut_down==-1: #### full until the bottom!
                        if cut_up>-1:
                            for c in range(r+1, cut_up+1):
                                system_dup[c]=system[c-1]
                        if cut_up==-1: #### full until the end (total full)
                            for c in range(r+1, N):
                                system_dup[c]=system[c-1]
                    if cut_down>-1:
                        p=up/(down+up)
                        if ran.uniform(0,1)<p:
                            for c in range(cut_down,r):
                                system_dup[c]=system[c+1]
                        else:
                            for c in range(r+1, cut_up+1):
                                system_dup[c]=system[c-1]
        
        if r==0: ### if the chosen cell is at the bottom
            if system[r+1]==-1:
                system_dup[r+1]=system[r]
            if system[r+1]>-1:
                up=0
                cut_up=-1
                for c in range(r+1,N):
                    if system[c]>-1:
                        up=up+1
                    if system[c]==-1:
                        cut_up=c
                        break
                if cut_up==-1: #### if the system is full
                    for c in range(r+1, N):
                        system_dup[c]=system[c-1]
                if cut_up>-1: ### if there is a hole, advance 1 until the hole
                    for c in range(r+1, cut_up+1):
                        system_dup[c]=system[c-1]

    return(system_dup)

##### END event of Duplication
####################################################

####################################################
##### event of Duplication IF the above position is not empty

def duplication_event_connex(system):
    
    N=len(system)
    system_dup=[]

    for i in range(0, N):
        system_dup.append(system[i])
    
    r=ran.randint(0,N-1) #### chose a cell to duplicate
    if system[r]>-1: #### there is a cell in the position
        if r>0:
            ##### cases of trivial movement
            if r <N-1 and system[r+1]>-1:
                if system[r-1]==-1: #### goes down by default if it can 
                    system_dup[r-1]=system[r]
                if r<N-1: 
                    if system[r-1]>-1 and system[r+1]==-1:
                        system_dup[r+1]=system[r]
            if r==N-1:
                if system[r-1]==-1: #### goes down by default if it can 
                    system_dup[r-1]=system[r]
                if r<N-1: 
                    if system[r-1]>-1 and system[r+1]==-1:
                        system_dup[r+1]=system[r]
                    
            ### neighbours in both sides
            
                if system[r-1]>-1 and system[r+1]>-1: 
                    up=0
                    down=0
                    cut_up=-1
                    cut_down=-1
                    
                    ### computing the number of successive cells without empty spaces
                    for c in range(r+1,N): ##### up cells
                        if system[c]>-1:
                            up=up+1
                        if system[c]==-1:
                            cut_up=c
                            break
                        
                    for c in range(0,r): ##### down cells
                        if system[r-1-c]>-1:
                            down=down+1
                        if system[r-1-c]==-1:
                            cut_down=r-1-c
                            break
                    
                    if cut_down==-1: #### full until the bottom!
                        if cut_up>-1:
                            for c in range(r+1, cut_up+1):
                                system_dup[c]=system[c-1]
                        if cut_up==-1: #### full until the end (total full)
                            for c in range(r+1, N):
                                system_dup[c]=system[c-1]
                    if cut_down>-1:
                        p=up/(down+up)
                        if ran.uniform(0,1)<p:
                            for c in range(cut_down,r):
                                system_dup[c]=system[c+1]
                        else:
                            for c in range(r+1, cut_up+1):
                                system_dup[c]=system[c-1]
        
        if r==0: ### if the chosen cell is at the bottom
            if system[r+1]==-1:
                system_dup[r+1]=system[r]
            if system[r+1]>-1:
                up=0
                cut_up=-1
                for c in range(r+1,N):
                    if system[c]>-1:
                        up=up+1
                    if system[c]==-1:
                        cut_up=c
                        break
                if cut_up==-1: #### if the system is full
                    for c in range(r+1, N):
                        system_dup[c]=system[c-1]
                if cut_up>-1: ### if there is a hole, advance 1 until the hole
                    for c in range(r+1, cut_up+1):
                        system_dup[c]=system[c-1]

    return(system_dup)

##### END event of Duplication
####################################################
  

####################################################
#### event of migration DO NOT CARE WHO IS MOVING
    
def migration_system(system, v):

    r=ran.randint(0,len(system)-1) #### chose a cell to perform a migration step towards position 0  
    hit_temp=system[r]
    if system[r]>-1 and r>0:
        prev_temp=system[r-1]
        if prev_temp==-1: #### free migration
            if ran.uniform(0,1)<v:
                system[r-1]=hit_temp
                system[r]=prev_temp
            
    return(system)

####################################################
#### event of migration/swap WE CARE WHO IS MOVING

def migration_swap_trajectory(system, v, q,cell): ### competition.Prob of swap if it is possible
    mov=0
    r=ran.randint(0,len(system)-1) #### chose a cell to perform a migration step towards position 0  
    hit_temp=system[r]
    if system[r]>-1 and r>0:
        prev_temp=system[r-1]
        if prev_temp==-1: #### free migration
            if ran.uniform(0,1)<v:
                system[r-1]=hit_temp
                system[r]=prev_temp
                if r==cell:
                    mov=-1
        if prev_temp>-1:
            if ran.uniform(0, 1)<q: #### competition with the cell at the bottom
                system[r-1]=hit_temp
                system[r]=prev_temp
            if r==cell:
                mov=-1
            if r-1==cell:
                mov=1

    return(system,mov)
    
#### END event of migration/swap    
####################################################
    
####################################################
##### event of Duplication

def duplication_event_trajectory(system, cell):
    
    mov=0
    N=len(system)
    system_dup=[]

    for i in range(0, N):
        system_dup.append(system[i])
    
    r=ran.randint(0,N-1) #### chose a cell to duplicate
    if system[r]>-1: #### there is a cell in the position
        if r>0:
            ##### cases of trivial movement
            if system[r-1]==-1: #### goes down by default if it can 
                system_dup[r-1]=system[r]
            if r<N-1: 
                if system[r-1]>-1 and system[r+1]==-1:
                    system_dup[r+1]=system[r]
                    
            ### neighbours in both sides
            
                if system[r-1]>-1 and system[r+1]>-1: 
                    up=0
                    down=0
                    cut_up=-1
                    cut_down=-1
                    
                    ### computing the number of successive cells without empty spaces
                    for c in range(r+1,N): ##### up cells
                        if system[c]>-1:
                            up=up+1
                        if system[c]==-1:
                            cut_up=c
                            break
                        
                    for c in range(0,r): ##### down cells
                        if system[r-1-c]>-1:
                            down=down+1
                        if system[r-1-c]==-1:
                            cut_down=r-1-c
                            break
                    
                    if cut_down==-1: #### full until the bottom!
                        if cut_up>-1:
                            for c in range(r+1, cut_up+1):
                                system_dup[c]=system[c-1]
                                if c==cell:
                                    mov=1

                        if cut_up==-1: #### full until the end (total full)
                            for c in range(r+1, N):
                                system_dup[c]=system[c-1]
                                if c==cell:
                                    mov=1

                    if cut_down>-1:
                        p=up/(down+up)
                        if ran.uniform(0,1)<p:
                            for c in range(cut_down,r):
                                system_dup[c]=system[c+1]
                                if c==cell:
                                    mov=-1

                        else:
                            for c in range(r+1, cut_up+1):
                                system_dup[c]=system[c-1]
                                if c==cell:
                                    mov=1
        
        if r==0: ### if the chosen cell is at thew bottom
            if system[r+1]==-1:
                system_dup[r+1]=system[r]
            if system[r+1]>-1:
                up=0
                cut_up=-1
                for c in range(r+1,N):
                    if system[c]>-1:
                        up=up+1
                    if system[c]==-1:
                        cut_up=c
                        break
                if cut_up==-1: #### if the system is full
                    for c in range(r+1, N):
                        system_dup[c]=system[c-1]
                        if c==cell:
                            mov=1

                if cut_up>-1: ### if there is a hole, advance 1 until the hole
                    for c in range(r+1, cut_up+1):
                        system_dup[c]=system[c-1]
                        if c==cell:
                            mov=1

    return(system_dup, mov)

##### END event of Duplication
####################################################

####################################################
##### event of Duplication

def duplication_event_trajectory_T(system, cell):
    
    mov=0
    N=len(system)
    system_dup=[]

    for i in range(0, N):
        system_dup.append(system[i])
    
    r=ran.randint(0,N-1) #### chose a cell to duplicate
    if system[r]>-1: #### there is a cell in the position
        if r>0:
            ##### cases of trivial movement
            if system[r-1]==-1: #### goes down by default if it can 
                system_dup[r-1]=system[r]
                
            if r==N-1 and system[r-1]>-1:
                
                cut_down=-1
                for c in range(0,r): ##### Find if it is full 
                    if system[r-c]==-1:
                        cut_down=r-c
                        break
                    
                if cut_down>-1: #### Hole before the bottom!
                    for c in range(cut_down,r):
                        system_dup[c]=system[c+1]
                        if c==cell:
                            mov=-1

            ##### cases of trivial movement
            if r<N-1 and system[r-1]>-1: 
                if system[r+1]==-1:
                    system_dup[r+1]=system[r]
                    
            ### neighbours in both sides            
                if system[r+1]>-1: 
                    up=0
                    down=0
                    cut_up=-1
                    cut_down=-1
                    
                    ### computing the number of successive cells without empty spaces
                    for c in range(r+1,N): ##### up cells
                        if system[c]>-1:
                            up=up+1
                        if system[c]==-1:
                            cut_up=c
                            break
                        
                    for c in range(0,r): ##### down cells
                        if system[r-c]>-1:
                            down=down+1
                        if system[r-c]==-1:
                            cut_down=r-c
                            break
                    
                    if cut_down==-1: #### full until the bottom!
                        
                        if cut_up>-1: ### hole somewhere up
                            for c in range(r+1, cut_up+1):
                                system_dup[c]=system[c-1]
                                if c==cell:
                                    mov=1

                        if cut_up==-1: #### full until the top (total full)
                            for c in range(r+1, N):
                                system_dup[c]=system[c-1]
                                if c==cell:
                                    mov=1

                    #If there is a hole down, it goes down
                    if cut_down>-1: #holes somewhere down
                        if cut_up==-1: ### hole somewhere up
                            for c in range(cut_down, r):
                                system_dup[c]=system[c+1]
                                if c==cell:
                                    mov=-1
        
        if r==0: ### if the chosen cell is at the bottom
            if system[r+1]==-1:
                system_dup[r+1]=system[r]
            if system[r+1]>-1:
                up=0
                cut_up=-1
                for c in range(r+1,N):
                    if system[c]>-1:
                        up=up+1
                    if system[c]==-1:
                        cut_up=c
                        break
                if cut_up==-1: #### if the system is full
                    for c in range(r+1, N):
                        system_dup[c]=system[c-1]
                        if c==cell:
                            mov=1

                if cut_up>-1: ### if there is a hole, advance 1 until the hole
                    for c in range(r+1, cut_up+1):
                        system_dup[c]=system[c-1]
                        if c==cell:
                            mov=1

    return(system_dup, mov)

##### END event of Duplication
####################################################

