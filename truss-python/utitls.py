#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provides utilities used by FE analysis.
  1. assembly: Global stiffness matrix assembly.
  2. solvedr: Solving the stiffness equations by the reduction approach.

Created on Sat May 9 17:39:00 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""

import numpy as np
import FEData as model
import copy
import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt



def assembly(e, ke):
    """
    Assemble element stiffness matrix.
    
    Args:
        e   : (int) Element number
        ke  : (numpy(nen*ndof,nen*ndof)) element stiffness matrix
    """
    model.K[np.ix_(model.LM[:,e], model.LM[:,e])] += ke

def solvedr():
    """
    Partition and solve the system of equations
        
    Returns:
        f_E : (numpy.array(nd,1)) Reaction force vector
    """
    beta=10000000000000000000
    nd = model.nd
    neq = model.neq
    K_e = copy.deepcopy(model.K)
    K_e1 = copy.deepcopy(model.K)
    #Tr=K_e1.trace
    
    K_e1[0,0]=beta
    K_e1[1,1]=beta
    K_e1[2,2]=beta
    K_e1[3,3]=beta
    
    K_E = model.K[0:nd, 0:nd]
    K_F = model.K[nd:neq, nd:neq]
    K_EF =model. K[0:nd, nd:neq]
    f_F = model.f[nd:neq]
    d_E = model.d[0:nd]
    d_F = model.d[nd:neq]

    f_F1 = model.f
    f_F1[0]=0
    f_F1[1]=0
    f_F1[2]=0
    f_F1[3]=0

    # solve for d_F
    d = np.linalg.solve(K_e1, f_F1) 
    d_E=d[0:nd]
    d_F=d[nd:neq]
    
    # reconstruct the global displacement d
    model.d = np.append(d_E,d_F)
    
    # compute the reaction r    
    f_E =K_E@d_E + K_EF@d_F
   
    #K_E@d_E1  K_EF@d_F1
                      
    # write to the workspace
    print('\nsolution d');  print(model.d)
    print('\nreaction f =', f_E)
    mag=1e+4


    for i in range(model.nel):
                XX1 = np.array([model.x[model.IEN[i, 0]-1] + mag*model.d[2*model.IEN[i, 0]-1-1], 
                               model.x[model.IEN[i, 1]-1] + mag*model.d[2*model.IEN[i, 1]-1-1]])
                YY1= np.array([model.y[model.IEN[i, 0]-1] + mag*model.d[2*model.IEN[i, 0]-1], 
                               model.y[model.IEN[i, 1]-1] + mag*model.d[2*model.IEN[i, 1]-1]])
                plt.plot(XX1, YY1, "red")
    plt.show()

    
    return f_E