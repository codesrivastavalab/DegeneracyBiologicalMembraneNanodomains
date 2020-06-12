import numpy as np
import matplotlib.pyplot as plt
import itertools 
from numpy.random import rand
import time
import random
from numba import jit

start = time.time()
data = np.loadtxt('kaisq-val-dppc-dopc.txt') #read the kai^2 values for the system 
ener_chi = np.loadtxt('energy-tabel-dppc-dopc.txt') #read the interaction energy table
chi_i = ener_chi[:,0] #column 1 of the energy table has kai^2_i
chi_j = ener_chi[:,1] #coulmn 2 of the energy table has kai^2_j
ener_ij_avg = ener_chi[:,3] #column 4 of the energy table has the interaction Hamiltonian
chi_ij_avg = ener_chi[:,2] #column 3 of the energy table has the kai^2_i and kai2_j written as a unique combination so as to correspond to the H_{ij} of the combination pair

chid = np.reshape(data, (13,14)) #convert the kai^2 values to a 13*14 lattice
random.shuffle(chid) #start from an initial random configuration
bins=10 #for histogram
original,chi_interval = np.histogram(chid,bins) #histogram of the original kai^2 values
config_ener = 0 #initializing energy
chi=chid 
#function to calculate energy
@jit
def calcEner(chi):
    energy = 0
    for i in range(13):
        for j in range(14):
	    #find the nearest nbrs kai^2_i and kai^2_j combination
            avg_cal_chi_i_plus = round((chi[i,j]/chi[(i+1)%13,j]) - (chi[(i+1)%13,j]/(chi[i,j])),13)
            avg_cal_chi_i_minus = round((chi[i,j]/chi[(i-1)%13,j]) - (chi[(i-1)%13,j]/(chi[i,j])),13)
            avg_cal_chi_j_plus = round((chi[i,j]/chi[i,(j+1)%14]) - (chi[i,(j+1)%14]/(chi[i,j])),13)
            avg_cal_chi_j_minus = round((chi[i,j]/chi[i,(j-1)%14]) - (chi[i,(j-1)%14]/(chi[i,j])),13)
            if avg_cal_chi_i_plus and avg_cal_chi_i_minus and avg_cal_chi_j_plus and avg_cal_chi_j_minus != 0:
                chi_ndx_1 = list(chi_ij_avg).index(avg_cal_chi_i_plus) #look for the interactions corresponding to the kai^2_i and kai^2_j combination in the enerngy table
                chi_ndx_2 = list(chi_ij_avg).index(avg_cal_chi_i_minus)
                chi_ndx_3 = list(chi_ij_avg).index(avg_cal_chi_j_plus)
                chi_ndx_4 = list(chi_ij_avg).index(avg_cal_chi_j_minus)
                ener_1 = ener_ij_avg[chi_ndx_1]
                ener_2 = ener_ij_avg[chi_ndx_2]
                ener_3 = ener_ij_avg[chi_ndx_3]
                ener_4 = ener_ij_avg[chi_ndx_4]
                ener = ener_1 + ener_2 + ener_3 + ener_4 #nearest neighbour interaction
                energy += ener #summing the total energy of the lattice configuration
    config_ener = energy #total energy
    return(config_ener)

beta= 1 #temperature to perform Simulated Annealing

#function for the Monte Carlo Runs
@jit
def mcmove(beta,chi,step,rejection):
    jj = np.random.randint(14)
    ii = np.random.randint(13)
    mm = np.random.randint(13)
    nn = np.random.randint(14)
    R = chi[ii,jj] #select random lattice site 1
    S = chi[mm,nn] #select another random lattice site 2
    E_old = calcEner(chi) #calulate energy for this configuration
    chi[mm,nn] = R #swap lattice site 1 with lattice site 2 and vice versa
    chi[ii,jj] = S
    E_new = calcEner(chi) #calculate energy of the configuration after the swapping
    cost = E_new - E_old #energy change
    step = step+1 #next step
    factor = np.exp(-cost*beta) #MC criteria
    randno = rand() #random number
    if randno < factor: #acceptance or rejection of the swap
        chi[mm,nn] = R
        chi[ii,jj] = S
    else:
        rejection = rejection + 1 #for acceptance ratio calculation
        chi[ii,jj] = R
        chi[mm,nn] = S
    return(chi,step,rejection,randno,factor,cost)

iteration = 250000 #number of steps
chi = chid

step = 0
rejection = 0
randno = 0
factor = 0

#perform the MC run and give out energy of the configuration, \delta E, random number, MC acceptance exp(-beta \delta E)
for u in range(iteration):
    chi,step,rejection,randno,factor,cost=mcmove(beta,chi,step,rejection)
    average_energy= calcEner(chi)
    print(average_energy,cost,randno,factor) #print energy for each iteration

print(step,acceptance)
plt.imshow(chi) #plot the lattice configuration at the energy of the run
plt.show()    
