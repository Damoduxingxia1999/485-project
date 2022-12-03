#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 14:37:07 2021

@author: admin
"""

import numpy as np
import random
import matplotlib.pyplot as plt
from Wisdom2 import chemical_potential,m,beta,kb

def minimum_image(r, lbox):
  return r - lbox*np.round(r / lbox)


def surface_lattice(nsites, lbox):
    '''
    create sites lattice
    
    Parameters
    ----------
    tiling : int
        the number of atoms
    L : int
        the length of the box

    Returns
    -------
    surface : 2D array(2,nsites)
        surface[0] is postion;surface[1] is state

    '''
    surface=np.zeros([nsites])
    for i in range(nsites):
                surface[i]=(i/nsites-0.5)*lbox
    return    surface


def check_full(surface):
    '''
    check which site is full

    Parameters
    ----------
    surface : 2D array(2,nsites)
        surface[0] is postion;surface[1] is state

    Returns
    -------
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    '''
    full=[]
    nonfull=[]
    nsites=surface.shape[0]

    for i in range(nsites):
        if surface[i]==1:
            full.append(i)
        else:
            nonfull.append(i)
    return full,nonfull
 
           
def my_neighbor_list(i, nsites):
  """ Find all neighbors of site (i, j).

  Args:
    i (int): site index along x
    j (int): site index along y
    N (int): number of sites along each dimension
  Return:
    list: a list of 2-tuples, [(i_left, j_left), (i_above, j_above),
     (i_right, j_right), (i_below, j_below)]
  """
  left   = i-1
  right  = i+1
  if (i==nsites-1):
      right  = 0
  if (i==0):
      left   = nsites-1
  return left, right

def E_neighbor(surface,i,w):
    '''
    calculate the decrease of energy caused by neigbor interaction when absorb or deabsorb the molecule)
    w is negative

    Parameters
    ----------
    surface : 2D array(2,nsites)
        surface[0] is postion;surface[1] is state
    i : int
        the postion of sites which is absorb or deaborb    
    w : float
        neighbor interaction parameter

    Returns
    -------
    energy : float
        energy change caused by neighbor effect after absorb or adabsorb

    '''

    left,right=my_neighbor_list(i, surface.shape[0])
    energy=0
    energy+=w*2*(surface[left]+surface[right])
    return energy


# def potential_system(surface,binding,w):
#     '''   
#     calculate the whole energy of system, including binding energy 
#     and neighbor interaction

#     Parameters 
#     ----------
#     surface : 2D array(2,nsites)
#         surface[0] is postion;surface[1] is state
#     binding : float
#         binding energy
#     w : float
#         neighbor interaction parameter

#     Returns
#     -------
#     energy : float

# '''
#     nsites=surface.shape[0]
#     full,nonfull=check_full(surface)
#     energy=len(full)*binding
#     for i in range(nsites):
#         neighbor=my_neighbor_list(i, nsites)       
#         for a in neighbor:
#             energy+=w*(surface[1,a[0]]+surface[1,a[1]])
#     return energy





# def check_potential(chemical_potential,beta,mass):
#     '''calculate the exponential part in the acceptence ratio'''
#     nabta=matter_wave_length(mass,beta)
#     return np.exp(beta*chemical_potential/nabta**3)

def acc_absorb(nsites,natom,binding,surface,beta,chemical_potential,i):
    En=E_neighbor(surface, i, w)
    print('1',np.exp(-beta*binding))
    print('2',np.exp(beta*chemical_potential))
    print('3',np.exp(-beta*En))
    print('4',(nsites-natom)*np.exp(-beta*binding)*np.exp(beta*chemical_potential)*np.exp(-beta*En))
    return (nsites-natom)*np.exp(-beta*binding)*np.exp(beta*chemical_potential)*np.exp(-beta*En)

def acc_desorb(nsites,natom,binding,surface,beta,chemical_potential,i):
    En=-E_neighbor(surface, i, w)
    return np.exp(-beta*En)/(nsites-natom+1)/np.exp(-beta*binding)/np.exp(beta*chemical_potential)

def grand_canonical(surface, lbox, beta, chemical_potential,mass,binding,w,full,nonfull):
    '''
    gand_canonical move
    

    Parameters
    ----------
    surface : 2D array(2,nsites)
        surface[0] is postion;surface[1] is state
    lbox : float
        length of  the whold sites
        cutt_off radius
    beta : float
        1/(k_b*T)
    chemical_potential : float
        chemical_potential get from Widom
    mass : int
        mass of Ne atom
    binding : float
        binding energy
    w : float
        neighbor interaction parametersurface : TYPE

    Returns
    -------
    s 1D array(nsites)
        the total number of biding atoms per sites


    '''
    nsites=surface.shape[0]
    natom=len(full)
    naccept_remove = 0
    nreject_remove = 0
    naccept_add = 0
    nreject_add = 0

    if (random.random()<0.5 and len(full)!=0):
        '''remove'''
        i=random.choice(full)
        surface_old=np.copy(surface)
        surface[i]=0
        acc=acc_absorb(nsites, natom, binding, surface, beta, chemical_potential, i)
        print('d acc',acc)
        if (random.random()>acc):
            surface[i]=surface_old[i]
            nreject_remove+=1
        else:
            naccept_remove+=1
    else:
        '''absorb'''
        surface_old=np.copy(surface)
        #U_old=potential_neighborsite_all(surface_old,b)
        i=random.choice(nonfull)
        surface[i]=1
        #U_new=potential_site_all(surface,lbox,rc)+b
        acc=acc_absorb(nsites, natom, binding, surface, beta, chemical_potential, i)
        print('absorb acc',acc)
        if (random.random()>acc):
            surface=surface_old
            nreject_add+=1
        else:
            naccept_add+=1
    return surface



nsites=10
lbox=1e-10

w=2.4*kb
binding=5*kb

surface=surface_lattice(nsites, lbox)
s=np.zeros([nsites])

ncycle=10
for _ in range(ncycle):
    full,nonfull=check_full(surface)
    surface=grand_canonical(surface, lbox, beta, chemical_potential,m,binding,w,full,nonfull)
    s+=surface
    sites=np.asarray(np.linspace(0,nsites-1,nsites))
#sites=np.stack(sites,surface,axis=0)
    fig = plt.figure()
    ax = plt.axes()
    ax.scatter(sites,surface)
    plt.show()
s=np.average(s)/ncycle
print (s)

