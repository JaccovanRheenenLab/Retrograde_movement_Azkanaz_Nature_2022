#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 12:02:49 2022

@author: bernat
"""

import matplotlib.pyplot as plt
import numpy as np
import random as ran
import Migration_to_kr_header as Mkr

### Location of the bottommost cell after ablation in time
### Using the retrograde movement dynamics and the SCB dynamics of Mkr header
### To reproduce figure 4 h,i

##Real data !!!

############### Histograms. Location of the bottommost cell

Hist_1_Real_SI=np.array([394, 48, 36, 18, 64])
Hist_1_Real_SI=Hist_1_Real_SI/sum(Hist_1_Real_SI)
Hist_1_Real_LI=np.array([353, 96, 62, 40, 60])
Hist_1_Real_LI=Hist_1_Real_LI/sum(Hist_1_Real_LI)

Hist_2_Real_SI=np.array([524,30,7,2, 6])
Hist_2_Real_SI=Hist_2_Real_SI/sum(Hist_2_Real_SI)
Hist_2_Real_LI=np.array([325, 89, 35, 20, 4])
Hist_2_Real_LI=Hist_2_Real_LI/sum(Hist_2_Real_LI)

Hist_3_Real_SI=np.array([378,3,1,0,0])
Hist_3_Real_SI=Hist_3_Real_SI/sum(Hist_3_Real_SI)
Hist_3_Real_LI=np.array([276,46,16,4,2])
Hist_3_Real_LI=Hist_3_Real_LI/sum(Hist_3_Real_LI)

Hist_4_Real_SI=np.array([287,4,0,0,0])
Hist_4_Real_SI=Hist_4_Real_SI/sum(Hist_4_Real_SI)
Hist_4_Real_LI=np.array([476,3,0,0,0])
Hist_4_Real_LI=Hist_4_Real_LI/sum(Hist_4_Real_LI)

###################

Pos_0_Real_SI=[]
Pos_0_Real_LI=[]
Pos_1_Real_SI=[]
Pos_1_Real_LI=[]
Pos_2_Real_SI=[]
Pos_2_Real_LI=[]
Pos_3_Real_SI=[]
Pos_3_Real_LI=[]
Pos_4_Real_SI=[]
Pos_4_Real_LI=[]

Pos_0_Real_SI=np.array([0, Hist_1_Real_SI[0], Hist_2_Real_SI[0], Hist_3_Real_SI[0], Hist_4_Real_SI[0]])
Pos_0_Real_LI=np.array([0, Hist_1_Real_LI[0], Hist_2_Real_LI[0], Hist_3_Real_LI[0], Hist_4_Real_SI[0]])

Pos_1_Real_SI=np.array([0, Hist_1_Real_SI[1], Hist_2_Real_SI[1], Hist_3_Real_SI[1], Hist_4_Real_SI[1]])
Pos_1_Real_LI=np.array([0, Hist_1_Real_LI[1], Hist_2_Real_LI[1], Hist_3_Real_LI[1], Hist_4_Real_SI[1]])

Pos_2_Real_SI=np.array([0, Hist_1_Real_SI[2], Hist_2_Real_SI[2], Hist_3_Real_SI[2], Hist_4_Real_SI[2]])
Pos_2_Real_LI=np.array([0, Hist_1_Real_LI[2], Hist_2_Real_LI[2], Hist_3_Real_LI[2], Hist_4_Real_SI[2]])

Pos_3_Real_SI=np.array([0, Hist_1_Real_SI[3], Hist_2_Real_SI[3], Hist_3_Real_SI[3], Hist_4_Real_SI[3]])
Pos_3_Real_LI=np.array([0, Hist_1_Real_LI[3], Hist_2_Real_LI[3], Hist_3_Real_LI[3], Hist_4_Real_SI[3]])

Pos_4_Real_SI=np.array([1, Hist_1_Real_SI[4], Hist_2_Real_SI[4], Hist_3_Real_SI[4], Hist_4_Real_SI[4]])
Pos_4_Real_LI=np.array([1, Hist_1_Real_LI[4], Hist_2_Real_LI[4], Hist_3_Real_LI[4], Hist_4_Real_SI[4]])

########################################################################
#############               SIMULATIONS         ########################
########################################################################

##################################################
######################## LI ######################
##################################################

############################################
##### rates and real time. Units days
kd=0.5
kr=kd
days=15
##### END rates and real time. Units days
############################################

N=5 ##### Size of the whole system
cells=4 #### Size of the crypt starting point
cell_0=4 #### starting point
steps=1000 #### 15 days

alpha=kr/kd #### prefactor of the v with respect to kd v=alpha*kd
qkd=days*kd*N/steps ### probability that a duplication event occurs per step
v=alpha*qkd ##### speed of the migration max=1/v timesteps
q=1 ### threshold for triggering swap

################################################
################# Numerical simulation starts
################################################
    
replicas=1000
plt.figure()

traj_glob=[]
density_glob=[]
Deepest_cell=[]

for u in range(0, replicas):
    cell=cell_0
    system=Mkr.create_system_1D_ablation(N,cell) #### Empty system

    traj=[]
    density=[]
    Deepest_cell.append([])
    
    for i in range(0, len(system)):
        if system[i]==0:
            cell=i
    traj.append(cell)
    
    count_dup=0
   
    for i in range(0, steps):
        mov=0
        system, mov=Mkr.migration_swap_trajectory(system, v,q,cell) ### Downwards migration
                                                                    ### Swap in case previous position is occupied
        cell=cell+mov
        traj.append(cell)
    
        mov=0
        system_dup=[] ### event of Duplication
        if ran.uniform(0,1)<qkd:
            system_dup, mov=Mkr.duplication_event_trajectory_T(system,cell)
            count_dup=count_dup+1
            system=[]### updating the system
            for j in range(len(system_dup)):
                system.append(system_dup[j])
                
        cell=cell+mov
        traj.append(cell)        
        Deepest_cell[u].append(Mkr.first_cell(system))
            
    traj_glob.append(traj)
 
###################################################
##### Statistics over trajectories
    
traj_av=[]
traj_std=[]
for k in range(0, len(traj_glob[0])):
    temp=[]
    for i in range (0, len(traj_glob)):
        temp.append(traj_glob[i][k]/(cells+1))
    traj_av.append(np.mean(temp))
    traj_std.append(np.std(temp))

traj_av=np.array(traj_av)
traj_std=np.array(traj_std)

##### END Statistics over trajectories
###################################################

###################################################
##### Plots

####### Plotting the distribution of the position of the last cell through time

Hist_first_cell_days=[]

for i in range(0, steps):
    if i%int(steps/14)==0:
        #print (i)
        Hist_temp=[]
        Hist_temp=np.zeros(N)#[0,0,0,0,0]
  
        for k in range(0, replicas):
            w=Deepest_cell[k][i]
            Hist_temp[w]=Hist_temp[w]+1
        if sum(Hist_temp[:N])>0:
            Hist_first_cell_days.append(np.array(Hist_temp[:N])/sum(Hist_temp[:N]))
        else:
            Hist_first_cell_days.append(np.array(Hist_temp[:N]))

####################################################

plt.figure(figsize=(9,7))
labels_N = ['Day 0', 'Day 02', 'Day 04', 'Day 07', 'Day 15']
labels= np.arange(len(labels_N))

Pos_0=np.array([0,Hist_first_cell_days[2][0], Hist_first_cell_days[4][0],Hist_first_cell_days[6][0],Hist_first_cell_days[14][0],])
Pos_1=np.array([0,Hist_first_cell_days[2][1], Hist_first_cell_days[4][1],Hist_first_cell_days[6][1],Hist_first_cell_days[14][1],])
Pos_2=np.array([0,Hist_first_cell_days[2][2], Hist_first_cell_days[4][2],Hist_first_cell_days[6][2],Hist_first_cell_days[14][2],])
Pos_3=np.array([0,Hist_first_cell_days[2][3], Hist_first_cell_days[4][3],Hist_first_cell_days[6][3],Hist_first_cell_days[14][3],])
Pos_4=np.array([1,Hist_first_cell_days[2][4], Hist_first_cell_days[4][4],Hist_first_cell_days[6][4],Hist_first_cell_days[14][4],])

width = 0.35       # the width of the bars: can also be len(x) sequence

plt.bar(labels-0.2, Pos_0, width, label='Pos_0', color='blue')
plt.bar(labels+0.2, Pos_0_Real_LI, width, color='blue')

plt.bar(labels-0.2, Pos_1, width, bottom=Pos_0, label='Pos_1', color='green')
plt.bar(labels+0.2, Pos_1_Real_LI, width,  bottom=Pos_0_Real_LI, color='green')

plt.bar(labels-0.2, Pos_2, width, bottom=Pos_0+Pos_1, label='Pos_2', color='orange')
plt.bar(labels+0.2, Pos_2_Real_LI, width, bottom=Pos_0_Real_LI+Pos_1_Real_LI, color='orange')

plt.bar(labels-0.2, Pos_3, width, bottom=Pos_0+Pos_1+Pos_2, label='Pos_3', color='gray')
plt.bar(labels+0.2, Pos_3_Real_LI, width, bottom=Pos_0_Real_LI+Pos_1_Real_LI+Pos_2_Real_LI, color='gray')

plt.bar(labels-0.2, Pos_4, width, bottom=Pos_0+Pos_1+Pos_2+Pos_3, color='k', label='Empty')
plt.bar(labels+0.2, Pos_4_Real_LI, width, bottom=Pos_0_Real_LI+Pos_1_Real_LI+Pos_2_Real_LI+Pos_3_Real_LI, color='k')

plt.ylim(0,1.35)
plt.ylabel('Fraction', fontsize=16)
plt.title('Bottom-most cell after ablation. Left: Model, Right Real. LI', fontsize=14)
plt.xticks(labels, labels_N)

plt.legend()

plt.show()

###### END Plots LI
####################################################

##################################################
######################## SI ######################
##################################################

############################################
##### rates and real time. Units days
kd=0.5
kr=2.5*kd
days=15
##### END rates and real time. Units days
############################################

N=5 ##### Size of the whole system
cells=4 #### Size of the crypt starting point
cell_0=4 #### starting point
steps=1000 #### 15 days

alpha=kr/kd #### prefactor of the v with respect to kd v=alpha*kd
qkd=days*kd*N/steps ### probability that a duplication event occurs per step
v=alpha*qkd ##### speed of the migration max=1/v timesteps
q=1 ### threshold for triggering swap

################################################
################# Numerical simulation starts
################################################
    
replicas=1000
plt.figure()

traj_glob=[]
density_glob=[]
Deepest_cell=[]

for u in range(0, replicas):
    cell=cell_0
    system=Mkr.create_system_1D_ablation(N,cell) #### Empty system

    traj=[]
    density=[]
    Deepest_cell.append([])
    
    for i in range(0, len(system)):
        if system[i]==0:
            cell=i
    traj.append(cell)
    
    count_dup=0
   
    for i in range(0, steps):
        mov=0
        system, mov=Mkr.migration_swap_trajectory(system, v,q,cell) ### Downwards migration
                                                                    ### Swap in case previous position is occupied
        cell=cell+mov
        traj.append(cell)
    
        mov=0
        system_dup=[] ### event of Duplication
        if ran.uniform(0,1)<qkd:
            system_dup, mov=Mkr.duplication_event_trajectory_T(system,cell)
            count_dup=count_dup+1
            system=[]### updating the system
            for j in range(len(system_dup)):
                system.append(system_dup[j])
                
        cell=cell+mov
        traj.append(cell)        
        Deepest_cell[u].append(Mkr.first_cell(system))
            
    traj_glob.append(traj)
 
###################################################
##### Statistics over trajectories
    
traj_av=[]
traj_std=[]
for k in range(0, len(traj_glob[0])):
    temp=[]
    for i in range (0, len(traj_glob)):
        temp.append(traj_glob[i][k]/(cells+1))
    traj_av.append(np.mean(temp))
    traj_std.append(np.std(temp))

traj_av=np.array(traj_av)
traj_std=np.array(traj_std)

##### END Statistics over trajectories
###################################################

###################################################
##### Plots

####### Plotting the distribution of the position of the last cell through time

Hist_first_cell_days=[]

for i in range(0, steps):
    if i%int(steps/14)==0:
        #print (i)
        Hist_temp=[]
        Hist_temp=np.zeros(N)#[0,0,0,0,0]
  
        for k in range(0, replicas):
            w=Deepest_cell[k][i]
            Hist_temp[w]=Hist_temp[w]+1
        if sum(Hist_temp[:N])>0:
            Hist_first_cell_days.append(np.array(Hist_temp[:N])/sum(Hist_temp[:N]))
        else:
            Hist_first_cell_days.append(np.array(Hist_temp[:N]))

####################################################

plt.figure(figsize=(9,7))
labels_N = ['Day 0', 'Day 02', 'Day 04', 'Day 07', 'Day 15']
labels= np.arange(len(labels_N))

Pos_0=np.array([0,Hist_first_cell_days[2][0], Hist_first_cell_days[4][0],Hist_first_cell_days[6][0],Hist_first_cell_days[14][0],])
Pos_1=np.array([0,Hist_first_cell_days[2][1], Hist_first_cell_days[4][1],Hist_first_cell_days[6][1],Hist_first_cell_days[14][1],])
Pos_2=np.array([0,Hist_first_cell_days[2][2], Hist_first_cell_days[4][2],Hist_first_cell_days[6][2],Hist_first_cell_days[14][2],])
Pos_3=np.array([0,Hist_first_cell_days[2][3], Hist_first_cell_days[4][3],Hist_first_cell_days[6][3],Hist_first_cell_days[14][3],])
Pos_4=np.array([1,Hist_first_cell_days[2][4], Hist_first_cell_days[4][4],Hist_first_cell_days[6][4],Hist_first_cell_days[14][4],])

width = 0.35       # the width of the bars: can also be len(x) sequence

plt.bar(labels-0.2, Pos_0, width, label='Pos_0', color='blue')
plt.bar(labels+0.2, Pos_0_Real_SI, width, color='blue')

plt.bar(labels-0.2, Pos_1, width, bottom=Pos_0, label='Pos_1', color='green')
plt.bar(labels+0.2, Pos_1_Real_SI, width,  bottom=Pos_0_Real_SI, color='green')

plt.bar(labels-0.2, Pos_2, width, bottom=Pos_0+Pos_1, label='Pos_2', color='orange')
plt.bar(labels+0.2, Pos_2_Real_SI, width, bottom=Pos_0_Real_SI+Pos_1_Real_SI, color='orange')

plt.bar(labels-0.2, Pos_3, width, bottom=Pos_0+Pos_1+Pos_2, label='Pos_3', color='gray')
plt.bar(labels+0.2, Pos_3_Real_SI, width, bottom=Pos_0_Real_SI+Pos_1_Real_SI+Pos_2_Real_SI, color='gray')

plt.bar(labels-0.2, Pos_4, width, bottom=Pos_0+Pos_1+Pos_2+Pos_3, color='k', label='Empty')
plt.bar(labels+0.2, Pos_4_Real_SI, width, bottom=Pos_0_Real_SI+Pos_1_Real_SI+Pos_2_Real_SI+Pos_3_Real_SI, color='k')

plt.ylim(0,1.35)
plt.ylabel('Fraction', fontsize=16)
plt.title('Bottom-most cell after ablation. Left: Model, Right Real. SI', fontsize=14)
plt.xticks(labels, labels_N)

plt.legend()

plt.show()

###### END Plots SI
####################################################
