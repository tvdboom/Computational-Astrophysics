# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 14:37:27 2017

@author: vansluijs
"""

import numpy as np
from matplotlib import pyplot as plt

# us tableau color scheme
# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
#Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)

def plot_aoeo():
    """
    Plot outer semi-major axis as a function of eccentricty for different
    mass loss rates
    """
    plt.figure()
    ax = plt.subplot(111)
    labels = ['A', 'C', 'B']
    xmin, xmax = (np.inf, -np.inf)
    for i in (0,2,1):
        
        # load the outer semi-major axis and eccentricity
        data = np.load('aoeo'+str(i)+'.npy')
        ao = data[:,0]
        eo = data[:,1]
        
        # plot the results
        plt.plot(ao, eo, color = tableau20[i*2], label = labels[i], lw = 1.5)
        plt.xlabel('Relative semi-major axis', size = 15)
        plt.ylabel('Relative eccentricity', size = 15)
        plt.legend(loc=3)
        plt.tight_layout()
        
        # right range
        if np.min(ao) < xmin: xmin = np.min(ao)
        if np.max(ao) > xmax: xmax = np.max(ao)
            
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.xlim(xmin, xmax)
    plt.savefig('plots/aeo_plot.png', dpi=300)
    plt.savefig('plots/aeo_plot.pdf')
    plt.show()

def plot_aiei():
    """
    Plot inner semi-major axis as a function of eccentricty for different
    mass loss rates
    """
    plt.figure()
    ax = plt.subplot(111)
    labels = ['A', 'C', 'B']
    xmin, xmax = (np.inf, -np.inf)
    for i in (0,2,1):
        
        # load the outer semi-major axis and eccentricity
        data = np.load('aiei'+str(i)+'.npy')
        ai = data[:,0]
        ei = data[:,1]
        
        # plot the results
        plt.plot(ai, ei, color = tableau20[i*2], label = labels[i], lw = 1.5)
        plt.xlabel('Relative semi-major axis', size = 15)
        plt.ylabel('Relative eccentricity', size = 15)
        plt.legend(loc=3)
        plt.tight_layout()
        
        # right range
        if np.min(ai) < xmin: xmin = np.min(ai)
        if np.max(ai) > xmax: xmax = np.max(ai)
            
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.xlim(xmin, xmax)
    plt.savefig('plots/aei_plot.png', dpi=300)
    plt.savefig('plots/aei_plot.pdf')
    plt.show()
    
def plot_massloss_timestep():
    """
    Plot massloss per timestep for the different simulations.
    """
    plt.figure()
    ax = plt.subplot(111)
    labels = ['A', 'C', 'B']
    xmin, xmax = (np.inf, -np.inf)
    for i in (0,2,1):
        
        # load the outer semi-major axis and eccentricity
        data = np.load('massloss'+str(i)+'.npy')
        timestep = np.arange(0, len(data[:,0]))+1
        massloss = data[:,1]
        
        # plot the results
        plt.plot(timestep, massloss, color = tableau20[i*2], label = labels[i], lw = 1.5)
        plt.xlabel('Time step', size = 15)
        plt.ylabel(r'Massloss [M$_{\odot}$ / Myr ]', size = 15)
        plt.legend(loc=3)
        plt.tight_layout()
        
        # right range
        if np.min(timestep) < xmin: xmin = np.min(timestep)
        if np.max(timestep) > xmax: xmax = np.max(timestep)
            
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.xlim(xmin, xmax)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.savefig('plots/massloss_plot.png', dpi=300)
    plt.savefig('plots/massloss_plot.pdf')
    plt.show()
    
def plot_massloss_time():
    """
    Plot massloss per timestep for the different simulations.
    """
    plt.figure()
    ax = plt.subplot(111)
    labels = ['A', 'C', 'B']
    xmin, xmax = (np.inf, -np.inf)
    for i in (0,2,1):
        
        # load the outer semi-major axis and eccentricity
        data = np.load('massloss'+str(i)+'.npy')
        time = data[:,0]
        dmdt = data[:,2]
        
        # plot the results
        plt.plot(time, dmdt, color = 'k', label = labels[i], lw = 1.5)
        plt.xlabel('Time [Myr]', size = 15)
        plt.ylabel(r'Massloss [M$_{\odot}$ / Myr ]', size = 15)
        plt.legend(loc=3)
        plt.tight_layout()
        
        # right range
        if np.min(time) < xmin: xmin = np.min(time)
        if np.max(time) > xmax: xmax = np.max(time)
            
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.xlim(xmin, xmax)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.savefig('plots/dmdt_plot.png', dpi=300)
    plt.savefig('plots/dmtdt_plot.pdf')
    plt.show()
    
def plot_time_timestep():
    """
    Plot time as a function of timestep.
    """
    plt.figure()
    ax = plt.subplot(111)
    labels = ['A', 'C', 'B']
    xmin, xmax = (np.inf, -np.inf)
    for i in (0,2,1):
        
        # load the outer semi-major axis and eccentricity
        data = np.load('massloss'+str(i)+'.npy')
        time = data[:,0]
        dt = np.diff(time)
        timestep = np.arange(1, len(data[:,0]))
        
        # plot the results
        plt.plot(timestep, dt, color = tableau20[i*2], label = labels[i], lw = 1.5)
        plt.xlabel('Timestep', size = 15)
        plt.ylabel('dt [Myr]', size = 15)
        plt.legend(loc=3)
        plt.tight_layout()
        
        # right range
        if np.min(timestep) < xmin: xmin = np.min(timestep)
        if np.max(timestep) > xmax: xmax = np.max(timestep)
            
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.xlim(xmin, xmax)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.savefig('plots/timetimestep_plot.png', dpi=300)
    plt.savefig('plots/timetimestep_plot.pdf')
    plt.show()
    
def plot_incoinci():
    """
    Plot time as a function of timestep.
    """
    plt.figure()
    ax = plt.subplot(111)
    labels = ['A', 'C', 'B']
    xmin, xmax = (np.inf, -np.inf)
    for i in (0,2,1):
        
        # load the outer semi-major axis and eccentricity
        data = np.load('incoinci'+str(i)+'.npy')
        inco = data[:,0]
        inci = data[:,1]
        
        # plot the results
        plt.plot(inco, inci, color = tableau20[i*2], label = labels[i], lw = 1.5)
        plt.xlabel('Relative outer inclination', size = 15)
        plt.ylabel('Relative innner inclination', size = 15)
        plt.legend(loc=3)
        plt.tight_layout()
        
        # right range
        if np.min(inco) < xmin: xmin = np.min(inco)
        if np.max(inco) > xmax: xmax = np.max(inco)
            
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.xlim(xmin, xmax)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.savefig('plots/incoinci_plot.png', dpi=300)
    plt.savefig('plots/incoicni_plot.pdf')
    plt.show()

# plot everything
#plot_aoeo()
#plot_aiei()
#plot_massloss_timestep()
#plot_massloss_time()
plot_time_timestep()
#plot_incoinci()

# print system times
labels = ['A','C','B']
for i in (0,2,1):
    data = np.load('time'+str(i)+'.npy')
    print ''
    print labels[i]
    print 'time: ',data[0],'s'
    print 'framework time: ',data[1],'s'
    print '% in amuse',data[2]
