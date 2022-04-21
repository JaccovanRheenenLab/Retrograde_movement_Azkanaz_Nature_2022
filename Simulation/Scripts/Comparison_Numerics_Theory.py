#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 10:50:30 2022

@author: bernat
"""

import numpy as np
import matplotlib.pyplot as plt
import Migration_to_kr_header as Mkr
    
###############################################################
### Test of the model predictions, contained in the file 
### "Model_Predictions_General.txt"
### Reproduce figure Ext. Data 9i
### Simulation parameters kr=2, kd=1
###############################################################

###Load simulation data

Simulation_data=np.loadtxt('Model_Predictions_General.txt')
ST=np.transpose(Simulation_data)

######## Size and parametrization
L=9
kr=2
kd=1
tday=np.log(2)/1.4
krkd=kr/kd

##### SDay: Integration results

SDay1_SI=[]
SDay2_SI=[]
SDay3_SI=[]
SDay4_SI=[]
SDay_inf=[]

choice=0

x=np.arange(0,9,1)

if choice==0:    
    for i in range (0,L):
        SDay1_SI.append(Mkr.Rel_psurvival_approx(i,0.01,krkd,L))
        SDay2_SI.append(Mkr.Rel_psurvival_approx(i,tday,krkd,L))
        SDay3_SI.append(Mkr.Rel_psurvival_approx(i,tday*2,krkd,L))
        SDay4_SI.append(Mkr.Rel_psurvival_approx(i,tday*3,krkd,L))
        SDay_inf.append(Mkr.p_gauss_infinity(i,krkd))

SDay1_SI=np.array(SDay1_SI)
SDay2_SI=np.array(SDay2_SI)
SDay3_SI=np.array(SDay3_SI)
SDay4_SI=np.array(SDay4_SI)
SDay_inf=np.array(SDay_inf)

################ Plotiing the figure with the comparison

plt.figure(figsize=(10,6))

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.title('Theory (solid lines) Vs simulation (circles)', fontsize=18)
plt.ylabel('Retention probability',fontsize=18)
plt.xlabel('Initial position',fontsize=18)
plt.plot(SDay1_SI*ST[0][0]/SDay1_SI[0], color='red')
plt.plot(ST[0], 'o',color='red', markersize=12, mfc='none')
plt.fill_between(x,ST[0]+ST[5],ST[0]-ST[5], color='red', alpha=0.2)
plt.plot(SDay2_SI*ST[1][0]/SDay2_SI[0], color='green')#, linestyle='--')
plt.plot(ST[1], 'o',color='green',markersize=12, mfc='none')
plt.fill_between(x,ST[1]+ST[6],ST[1]-ST[6], color='green', alpha=0.2)
plt.plot(SDay3_SI*ST[2][0]/SDay3_SI[0], color='blue')#,linestyle='--')
plt.plot(ST[2], 'o',color='blue',markersize=12, mfc='none')
plt.fill_between(x,ST[2]+ST[7],ST[2]-ST[7], color='blue', alpha=0.2)
plt.plot(SDay4_SI*ST[3][0]/SDay4_SI[0],color='orange')#,linestyle='--')
plt.plot(ST[3],'o', color='orange', markersize=12, mfc='none')
plt.fill_between(x,ST[3]+ST[8],ST[3]-ST[8], color='orange', alpha=0.2)
plt.plot(SDay_inf/np.sum(SDay_inf),color='red')#,linestyle='--')
plt.plot(ST[4]/np.sum(ST[4]),'o',color='red',markersize=12, mfc='none')
plt.fill_between(x,5*ST[4]+5*ST[9],5*ST[4]-5*ST[9], color='red', alpha=0.2)
plt.ylim(0,1.1)

#plt.savefig('Theory_Vs_Simulation.svg', dpi=600, format='svg')

plt.show()

