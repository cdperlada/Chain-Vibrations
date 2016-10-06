# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 12:43:38 2016

@author: cdperlada
"""

#!/usr/bin/env python

"""
Code for Chain Vibrations (W.Kinzel/G.Reents, Physics by Computer)

This code is based on chain.m listed in Appendix E of the book and will
replicate Fig. 2.12 of the book.
"""

__author__ = "Christian Alis"
__credits__ = "W.Kinzel/G.Reents"

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigvals, inv, eig

# what does ** in **kwargs do?
def main(f=1.0, m1=0.4, m2=1.0, **kwargs):
    
    '''
    main(f=1.0, m1=0.4, m2=1.0, **kwargs)
	Plot frequencies of four eigenmodes of a linear chain as a function of the wave number (q). 
	
	Parameters
	----------
	f : float
		Spring constant
	m1, m2: float
		masses of the atoms
		
	Returns
	-------
	eigensys: array
		The eigenvalues and eigenmodes
    '''
    
    print "\n Oscillator Chain \n"
    F = lambda q: np.array([[2*f            ,  -f,   0, -f*np.exp(-1j*q)],
                               [-f             , 2*f,  -f,                0],
                               [0              ,  -f, 2*f,               -f],
                               [-f*np.exp(1j*q),   0,  -f,              2*f]])
    # will cyclic rearrangements of m1 and m2 below change the results? Try it
    '''	No, the results were not changed since we have a periodic repetition 
    of the linear chains so a cylic rearrangement would not the dynamics'''
    
    # what will be the dimensions/shape of massmat?
    '''Massmat is a 4x4 diagonal matrix'''
    
    massmat = np.diag([m1, m1, m1, m2])
    # what happens if inv(massmat) * mat1(q) is used instead of 
    # inv(massmat).dot(mat1(q))?
    '''They will multiplied as arrays, element by element, instead of
		matrix multiplication.'''
    
    mat2 = lambda q: inv(massmat).dot(F(q))
    # what is the python type of kwargs?
    '''Kwargs is a dictionary'''
    
    plot_step = kwargs.get('plot_step', np.pi/50)
    # complete the following line. you should use plot_step
    x_axis = np.arange(-np.pi,np.pi+plot_step,plot_step)
    # what is the difference between eigvals() and eig()?
    '''eigvals() returns the eigenvalues while eig() returns
		both the eigenvalues and the correspnding eigenvectors'''
    
    # replace the following line to use eig() instead of eigvals()
    eigenlist = [eig(mat2(x))[0] for x in x_axis]
    # complete the following lines:
    plt.plot(x_axis, np.sqrt(eigenlist),'.')
    plt.xticks(np.arange(-np.pi,np.pi+np.pi/2,np.pi/2),[r"$\pi$", r"$\frac{\pi}{2}$", "0", r"$\frac{\pi}{2}$",
                r"$\pi$"])
    plt.xlim([-np.pi,np.pi])
    plt.xlabel('q',fontsize='x-large')
    plt.ylabel('$\omega$',fontsize='x-large')
    plt.tick_params('x', labelsize='x-large')
    plt.show()

    # why is the argument of mat2(), 0.0?
    '''We only need to consider the unit cell.'''
     
    # what does the parameter of mat2() mean?
    '''The parameter of mat2() is q.'''
     
    eigensys = eig(mat2(0.0))
    return eigensys
    
    # additional: 
    # * rename mat1 and mat2 to be more descriptive
    '''In the book, mat1 is named as F. mat2 is the product of the inverse of the matrix massmat
		and the matrix F)'''
    
    # * describe the result if the values of m1 and m2 are interchanged
    '''When the values of m1 and m2 are interchanged, the bands are lower.'''
    
    # * describe the effect of different values of f
    '''As f (spring constant) increases, the angular frequency omega increases.'''
    
    # optional challenging exercise: animate the eigenmodes as in Fig. 2.13

if __name__ == "__main__":
    print main()