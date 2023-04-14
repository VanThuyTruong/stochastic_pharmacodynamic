# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 13:39:02 2022

@author: kppb302
"""

# Nov 2022
# death and birth of cells
# drug dose and recovery
from random import random
from pylab import subplot,hist,plot,show,legend,xlabel,xticks,ylabel,xlim,ylim,savefig,figure,tight_layout,semilogy
from numpy import exp,random as nprandom
import matplotlib
from copy import deepcopy
import numpy as np

nreal = 1
tmax = 10
N = 100
mu = 1.0
lam = np.arange(0.2,1,0.01)#0.2#1.0 =[]

delta = 2.5
alpha = 1.5
T = 1#3.0/delta
num_cycles=10
dt = 0.01

Alow = delta/mu*(exp(mu/delta)-1)

x = mu/delta
Astar = exp(x)*(4**x-0.25)/(1+x) + (exp(x)-1)/(4*x)

ext_lamb=np.empty((0,2))
figure(figsize=(7,4))
params = {'font.family': 'serif'}
matplotlib.rcParams.update(params)

class Cell(object):
   ''' each cell is born with an initial k value '''
   def __init__(self):
      self.k = random()

def drate(k):
   ''' death rate = mu if k(t) < 0.25 '''
   rate = 0
   if k < 0.25:
#      rate = mu*(1-4*kt)
      rate = mu
   return rate

def brate(k):
   ''' division rate = lambda if k(t) > 0.5 '''
   rate = 0
   if k > 0.5:
      rate = lamb
   return rate

def onedose(celllist,T,a):
   ''' start with ncells, drug applied to time T 
   k(t) = k(0)exp(-delta*t)'''

   t = a
   print('dose t',t)

   ncells = len(celllist)

   while t < a+T:
      ranu = nprandom.random(ncells)
      drates = [drate(exp(-delta*t)*cell.k) for cell in celllist]
      celllist = [cell for i,cell in enumerate(celllist) if ranu[i] > drates[i]*dt] ##cells in the death pool die
      ranu = nprandom.random(ncells)
      brates = [brate(exp(-delta*t)*cell.k) for cell in celllist]    ##cells in the birth pool get a decrease in their k value
      newcellist = [deepcopy(cell) for i,cell in enumerate(celllist) if ranu[i] < brates[i]*dt]  ### cells from the birth pool giving birth despite k value decrease if k value is still high
      celllist += newcellist  ##add new daughter cells to celllist
      t += dt
      tlist.append(t)
      nlist.append(ncells)
      ncells = len(celllist)
      tlist.append(t)
      nlist.append(ncells)
      print('t',t)
   return celllist

def onerecovery(celllist,T,a):
   ''' start with ncells after drug applied to time T.
   k(t) = k(0) - (k(0)-k(T))exp(-delta*t)'''

   t = a
   print('recovery t',t)
   ncells = len(celllist)

   while t < a+2*T:  
      ranu = nprandom.random(ncells)
      kfac = 1 - exp(-alpha*(t-T))*(1-exp(-delta*T))   #cells recover and regain their initial k/perk value
      drates = [drate(kfac*cell.k) for cell in celllist]
      celllist = [cell for i,cell in enumerate(celllist) if ranu[i] > drates[i]*dt] #some cells still die
      ranu = nprandom.random(ncells)
      brates = [brate(kfac*cell.k) for cell in celllist]
      newcellist = [deepcopy(cell) for i,cell in enumerate(celllist) if ranu[i] < brates[i]*dt] #birth of daughter cells
      celllist += newcellist
      t += dt
      tlist.append(t)
      nlist.append(ncells)
      ncells = len(celllist)
      tlist.append(t)
      nlist.append(ncells)
   return celllist

def onecourse(ncells):
   ''' start with N cells, administer drug dose and relaxation '''
   print('start')
   celllist = [Cell() for i in range(ncells)]  #get all the initial cells into celllist
   for i in range(num_cycles):
       celllist = onedose(celllist,T,i*3)   #giving drug
       celllist = onerecovery(celllist,T,i*3+1)  #drug pause/recovery period


   return celllist

for i in range(nreal):
   lamb=lam[i]
   tlist,nlist = [0],[N]
   celllist = onecourse(N)

    #mybins = [t*0.1 for t in range(10*int(max(exttimelist)+2))]
   myt = [t*0.1 for t in range(1,10*tmax)]
   a=np.array([lamb,nlist[-1]])
   ext_lamb=np.vstack((ext_lamb, a)) 
#subplot(211)
   plot(tlist,nlist,label='{j}'.format(j=lamb))
   ylim([0,max(nlist)])
   xlim([0,max(tlist)])
   xlabel('$t$')
   ylabel('number of cells')
    
   #plot([T,T],[0,N],'--k',lw=0.5)
   legend()
   
#savefig('recovery.png')

#show()
np.save('ext_lamb',ext_lamb)

