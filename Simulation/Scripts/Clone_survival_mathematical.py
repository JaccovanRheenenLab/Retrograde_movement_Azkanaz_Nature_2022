#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 14:28:57 2022

@author: bernat
"""

import numpy as np
import matplotlib.pyplot as plt
import Migration_to_kr_header as Mkr
from scipy.optimize import curve_fit
    
##################################################
##################################################
##################################################
###########   REAL DATA ANALYSIS   ###############
##################################################
##################################################
##################################################

days=4

##### SI					
#		day 1	day 2	day 3	day 4
C1=np.array([1,	0.987341772151899,	0.925373134328358,	0.897959183673469])
C2=np.array([1,	0.955056179775281,	0.89873417721519,	0.87037037037037])
B1=np.array([1,	0.9,	0.619047619047619,	0.527777777777778])
B2=np.array([1,	0.527777777777778,	0.217391304347826,	0.125])

SI_infinity=np.array([0.260,	0.146,0.104,	0.046875])
SI_infinity=SI_infinity/sum(SI_infinity)

day_1SI=[]
day_2SI=[]
day_3SI=[]
day_4SI=[]

day_1SI=np.array([C1[0], C2[0],B1[0],B2[0]])	
day_2SI=np.array([C1[1], C2[1],B1[1],B2[1]])	
day_3SI=np.array([C1[2], C2[2],B1[2],B2[2]])
day_4SI=np.array([C1[3], C2[3],B1[3],B2[3]])
 
##### LI					
#		day 1	day 2	day 3	day 4
C1C=np.array([1,	0.964912280701754,	0.974358974358974,	0.927272727272727])
C2C=np.array([1,	0.847826086956522,	0.716417910447761,	0.510204081632653])
B1C=np.array([1,	0.617021276595745,	0.395833333333333,	0.3125])
B2C=np.array([	1,	0.264150943396226,	0.180327868852459,	0.127659574468085])

cecum_infinity=np.array([0.2058, 0.045,	0,	0])
cecum_infinity=cecum_infinity/sum(cecum_infinity)

day_1C=[]
day_2C=[]
day_3C=[]
day_4C=[]

day_1C=np.array([C1C[0], C2C[0],B1C[0],B2C[0]])	
day_2C=np.array([C1C[1], C2C[1],B1C[1],B2C[1]])	
day_3C=np.array([C1C[2], C2C[2],B1C[2],B2C[2]])	
day_4C=np.array([C1C[3], C2C[3],B1C[3],B2C[3]])	

####################################################
######## SI Confrontation with the numerics.... ####
####################################################

##################################################
##################################################
####### Fiting parameters ###########
##################################################
##################################################

levels=np.array([0,1,2,3])    

SI_infinity_U=SI_infinity/SI_infinity[0]
cecum_infinity_U=cecum_infinity/cecum_infinity[0]

krkd_SI, rest1 = curve_fit(Mkr.p_gauss_infinity, levels, SI_infinity_U)
krkd_CC, rest2 = curve_fit(Mkr.p_gauss_infinity, levels, cecum_infinity_U)

##################################################
##################################################
####### End Fiting parameters ###########
##################################################
##################################################

L=4
kr=1
kd=1
tday=1/4
krkd=krkd_SI[0]

##### SDay: Integration results

SDay1_SI=[]
SDay2_SI=[]
SDay3_SI=[]
SDay4_SI=[]
SDay_inf=[]

choice=0

if choice==0:    
    for i in range (0,L):
        SDay1_SI.append(Mkr.Rel_psurvival_approx(i,0.01,krkd,L))
        SDay2_SI.append(Mkr.Rel_psurvival_approx(i,tday,krkd,L))
        SDay3_SI.append(Mkr.Rel_psurvival_approx(i,tday*2,krkd,L))
        SDay4_SI.append(Mkr.Rel_psurvival_approx(i,tday*3,krkd,L))
        SDay_inf.append(Mkr.p_gauss_infinity(i,krkd))

if choice==1:
    for i in range (0,L):
        SDay1_SI.append(Mkr.Rel_psurvival(i,0.01,krkd,L))
        SDay2_SI.append(Mkr.Rel_psurvival(i,tday,krkd,L))
        SDay3_SI.append(Mkr.Rel_psurvival(i,tday*2,krkd,L))
        SDay4_SI.append(Mkr.Rel_psurvival(i,tday*3,krkd,L))
        SDay_inf.append(Mkr.p_gauss_infinity(i,krkd))
        
if choice==2:
    for i in range (0,L):
        SDay1_SI.append(Mkr.Rel_psurvival_Norm(i,0.01,krkd,L))
        SDay2_SI.append(Mkr.Rel_psurvival_Norm(i,tday,krkd,L))
        SDay3_SI.append(Mkr.Rel_psurvival_Norm(i,tday*2,krkd,L))
        SDay4_SI.append(Mkr.Rel_psurvival_Norm(i,tday*3,krkd,L))

SDay1_SI=np.array(SDay1_SI)
SDay2_SI=np.array(SDay2_SI)
SDay3_SI=np.array(SDay3_SI)
SDay4_SI=np.array(SDay4_SI)
SDay_inf=np.array(SDay_inf)
SDay_inf=(SDay_inf/SDay_inf[0])*SI_infinity[0]

plt.figure(figsize=(20,6))

plt.subplot(1,2,1)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.title('SI', fontsize=18)
plt.ylabel('Retention probability',fontsize=18)
plt.xlabel('Position in crypt',fontsize=18)

plt.plot(day_2SI[0]*SDay2_SI/SDay2_SI[0], color='green')
plt.plot(day_2SI, 'o',markersize=12, color='green')

plt.plot(day_3SI[0]*SDay3_SI/SDay3_SI[0], color='blue')
plt.plot(day_3SI, 'o',markersize=10, color='blue')

plt.plot(day_4SI[0]*SDay4_SI/SDay4_SI[0],color='orange')
plt.plot(day_4SI, 'o',markersize=10,  color='orange')

plt.plot(SDay_inf,color='red')
plt.plot(SI_infinity, 'o',markersize=10,  color='red')


#####################################################
###### LI Confrontation with the numerics....    ####
#####################################################

tday=1/2
krkd=krkd_CC[0]

SDay1_C=[]
SDay2_C=[]
SDay3_C=[]
SDay4_C=[]
SDay_inf_C=[]

for i in range (0,L):
    SDay1_C.append(Mkr.Rel_psurvival_approx(i,0.01,krkd,L))
    SDay2_C.append(Mkr.Rel_psurvival_approx(i,tday,krkd,L))
    SDay3_C.append(Mkr.Rel_psurvival_approx(i,tday*2,krkd,L))
    SDay4_C.append(Mkr.Rel_psurvival_approx(i,tday*3,krkd,L))
    SDay_inf_C.append(Mkr.p_gauss_infinity(i,krkd))

SDay1_C=np.array(SDay1_C)
SDay2_C=np.array(SDay2_C)
SDay3_C=np.array(SDay3_C)
SDay4_C=np.array(SDay4_C)
SDay_inf_C=np.array(SDay_inf_C)
SDay_inf_C=(SDay_inf_C/SDay_inf_C[0])*cecum_infinity[0]

plt.subplot(1,2,2)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Position in crypt',fontsize=18)


plt.title('LI', fontsize=18)

plt.plot(day_2C[0]*SDay2_C/SDay2_C[0], color='green')
plt.plot(day_2C, 'o',markersize=10, color='green')

plt.plot(day_3C[0]*SDay3_C/SDay3_C[0], color='blue')
plt.plot(day_3C, 'o',markersize=10, color='blue')

plt.plot(day_4C[0]*SDay4_C/SDay4_C[0],color='orange')
plt.plot(day_4C, 'o',markersize=10, color='orange')

plt.plot(SDay_inf_C,color='red')
plt.plot(cecum_infinity, 'o',markersize=10, color='red')

plt.ylim(-0.1,1.1)
plt.show()
