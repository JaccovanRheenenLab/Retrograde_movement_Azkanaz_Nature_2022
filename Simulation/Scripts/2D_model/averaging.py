import sys
from numpy import *
import re,os, os.path
import numpy as np
import pickle
import glob
from scipy.stats import kendalltau
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


T=5
st=(1-0.95)/2.0

data=np.zeros((5,4,1000))
for file in glob.iglob("surv*"):
    file=open(file)
    lines=file.readlines()
    file.close
    c=-1
    for i in range(len(lines)):
         cols=lines[i].strip().split()
         if float(cols[0])==0:
            c+=1
         if float(cols[0])==0*T:
             for j in range(4):
                data[0][j][c]=1
         if float(cols[0])==1*T:
             for j in range(4):
                 data[1][j][c]=float(cols[1+j])
         if float(cols[0])==2*T:
             for j in range(4):
                 data[2][j][c]=float(cols[1+j])
         if float(cols[0])==3*T:
             for j in range(4):
                 data[3][j][c]=float(cols[1+j])
         if float(cols[0])==56*T:
            for j in range(4):
                data[4][j][c]=float(cols[1+j])


sorted=np.sort(data,2)
av=np.mean(data,2)
std=np.std(data,2)

max=np.zeros((5,4))
min=np.zeros((5,4))
for i in range(5):
    for j in range(4):
        max[i][j]=sorted[i][j][int(st*1000)]
        min[i][j]=sorted[i][j][int((1-st)*1000)]

# the output is for  each time point considered in the experiment (day 1, 2, 3, 4, 56) the average, standard deviation, upper 95% confidence bound and lower 95% confidence bound the clonal survival probability as a function of starting position (1, 2, 3, 4) 

f=open('model_pos.txt','w')
for j in range(4):
    for i in range(5):
        f.write(str(av[i][j]))
        f.write(' ')
    for i in range(5):
        f.write(str(std[i][j]))
        f.write(' ')
    for i in range(5):
        f.write(str(max[i][j]))
        f.write(' ')
    for i in range(5):
        f.write(str(min[i][j]))
        f.write(' ')
    f.write('\n')
f.close()





