"""
CAP assignment 5 - The lost Solar siblings

Lennart van Sluijs and Marco TvdBoom
"""

from amuse.lab import *
from amuse.couple import bridge
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.gridspec as gridspec
from amuse.community.fractalcluster.interface import new_fractal_cluster_model
import matplotlib.patches
import scipy.stats

"""
Tableau is a color scheme used for the plots here.
"""
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

# Class from AMUSE github example files
class MilkyWay_galaxy(object):
    """
    Create the Milky way Galactic potential.
    """
    def get_gravity_at_point(self, eps, x,y,z):
        phi_0 = self.get_potential_at_point(eps, x,y,z)
        grav = AdaptingVectorQuantity()
        dpos = 0.001*(x**2+y**2+z**2).sqrt()
        phi_dx = self.get_potential_at_point(0,x+dpos,y,z) - phi_0
        phi_dy = self.get_potential_at_point(0,x,y+dpos,z) - phi_0
        phi_dz = self.get_potential_at_point(0,x,y, z+dpos) - phi_0
        return phi_dx/dpos, phi_dy/dpos, phi_dz/dpos

    def disk_or_bulge_potentials(self, x, y, z, a, b, mass):
        r2 = x**2 + y**2
        b2 = (a + (z**2 + b**2).sqrt())**2
        return constants.G * mass / (r2 + b2).sqrt()

    def halo_potential(self, x,y,z, Mc=5.0E+10|units.MSun, Rc=1.0|units.kpc**2):
        r=(x**2 + y**2 + z**2).sqrt()
        rr = (r/Rc)
        return -constants.G * (Mc/Rc)*(0.5*np.log(1 +rr**2) + np.arctan(rr)/rr)

    def get_potential_at_point(self, eps, x, y, z):
        pot_disk = self.disk_or_bulge_potentials(x,y,z, 
            0.0|units.kpc, 0.277|units.kpc, 1.12E+10|units.MSun) 
        pot_bulge = self.disk_or_bulge_potentials(x,y,z, 
            3.7|units.kpc, 0.20|units.kpc, 8.07E+10|units.MSun) 
        pot_halo = self.halo_potential(x,y,z, 
            Mc=5.0E+10|units.MSun, Rc=6.0|units.kpc)
        return pot_disk + pot_bulge + pot_halo

def sim_Sun_and_siblings(N = 50, perc = 0.01, t_end = 4.56 | units.Gyr,
                         dt_save = 0.2 | units.Myr, dt = 0.01 | units.Myr,
                         pointparticles = True, dirname = 'None'):
    """
    Simulate Sun and siblings in the Galactic potential.
    """
    
    if pointparticles:
        """
        Create a set of point partciles to integrate around the galactic center.
        """
        
        # for the pointparticles
        if dirname == 'None': dirname = './pointparticles/'
        
        #create outputfolder for this star
        outputfolder = os.path.join(dirname)
        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)
        
        # initalize bodies
        bodies = Particles(N+1)
        
        # create the Sun
        Sun = bodies[0]
        Sun.mass = 1e-9 | units.kg # assume Sun is point mass moving in the
        Sun.radius = 0 | units.RSun # galactic potential
        Sun.velocity = (-1.35, -232.1, -7.41) | units.kms
        Sun.position = (-8400, 0.0, 17.0) | units.parsec
    
        # create a Gaussian cloud of stars as point particles
        siblings = bodies[1:]
        siblings.mass = 1e-9 | units.kg 
        siblings.radius = 0 | units.RSun
        
        # create a cloud of particles
        siblings.x = np.random.normal(Sun.x.value_in(units.parsec), abs(Sun.x.value_in(units.parsec)) * perc, N) | units.parsec
        siblings.y = 0 | units.parsec # no spread in y-direction since we start at y=0
        siblings.z = np.random.normal(Sun.z.value_in(units.parsec), abs(Sun.z.value_in(units.parsec)) * perc, N) | units.parsec
        
        siblings.vx = np.random.normal(Sun.vx.value_in(units.kms), abs(Sun.vx.value_in(units.kms)) * perc, N) | units.kms
        siblings.vy = np.random.normal(Sun.vy.value_in(units.kms), abs(Sun.vy.value_in(units.kms)) * perc, N) | units.kms
        siblings.vz = np.random.normal(Sun.vz.value_in(units.kms), abs(Sun.vz.value_in(units.kms)) * perc, N) | units.kms
    
        # Integrator
        converter = nbody_system.nbody_to_si(bodies.mass.sum(), Sun.x)
        
        bodies.velocity *= -1 # Integrate backwards in time

        # integrator to use
        system = ph4(converter)

    else:
        """
        Create an open star cluster to integrate around the galactic center.
        """
        
        # if dirname not specified
        if dirname == 'None': dirname = './opencluster/'
        
        #create outputfolder for this star
        outputfolder = os.path.join(dirname)
        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)
        
        # parameters to use
        Rvir = 1.0 | units.parsec
        Qvir = 0.5
        
        # create masses
        Mmax = 100 | units.MSun
        masses = new_kroupa_mass_distribution(N+1, Mmax) # cluster with N other stars
        Mtot_init = masses.sum()
        
        converter = nbody_system.nbody_to_si(Mtot_init, Rvir)

        # assign particles
        bodies = new_plummer_model(N+1, converter)
        bodies.scale_to_standard(converter, virial_ratio=Qvir)
        bodies.mass = masses
        Sun = bodies[0]
        Sun.mass = 1 | units.MSun # first body becomes the Sun
        bodies.radius = 0 | units.RSun
        bodies.velocity += (30.6471439995, -200.098514469, 2.52218463835) | units.kms # our birthplace of the Sun
        bodies.position += (9499.1059096, 1999.79260586, 80.1313697564) | units.parsec
	
        # integrator with softening length scale        
        system = BHTree(converter)
        system.parameters.epsilon_squared = ( 1 | units.parsec )**2
        
    system.particles.add_particle(bodies)
    channel_from_system = system.particles.new_channel_to(bodies)

    gravity = bridge.Bridge()
    gravity.add_system(system, (MilkyWay_galaxy(),) )
    gravity.timestep = dt
    
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    # duration of simulation
    time = 0 | units.yr
    time_save = 0 | units.yr
    
    # allocate memory for positions and velocities
    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []
    
    while time < t_end:
        
        print 'Integrating...', time
        time += dt
        time_save += dt 

        gravity.evolve_model(time)
        channel_from_system.copy()

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        print "T=", time, "M=", bodies.mass.sum(), 
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
        Etot_prev = Etot
        
        if time_save > dt_save:

            # save the postion and velocities of all objects
            x = np.append(x, bodies.x.value_in(units.parsec))
            y = np.append(y, bodies.y.value_in(units.parsec))
            z = np.append(z, bodies.z.value_in(units.parsec))
            vx = np.append(vx, bodies.vx.value_in(units.kms))
            vy = np.append(vy, bodies.vy.value_in(units.kms))
            vz = np.append(vz, bodies.vz.value_in(units.kms))
            print 'Appended velocities and positions.'

            # reset timer
            time_save = 0 | units.yr
    
    # explicitly save last position
    x = np.append(x, bodies.x.value_in(units.parsec))
    y = np.append(y, bodies.y.value_in(units.parsec))
    z = np.append(z, bodies.z.value_in(units.parsec))
    vx = np.append(vx, bodies.vx.value_in(units.kms))
    vy = np.append(vy, bodies.vy.value_in(units.kms))
    vz = np.append(vz, bodies.vz.value_in(units.kms))
    print 'Appended final postions.'
 
    gravity.stop()
    
    # Correct the inverse velocity
    if pointparticles: bodies.velocity *= -1

    # plot the positions of all particles
    x = x.reshape((N+1, x.size/(N+1)))
    y = y.reshape((N+1, y.size/(N+1)))
    z = z.reshape((N+1, z.size/(N+1)))
    vx = vx.reshape((N+1, vx.size/(N+1)))
    vy = vy.reshape((N+1, vy.size/(N+1)))
    vz = vz.reshape((N+1, vz.size/(N+1)))
    
    # save positions seperately
    np.save(os.path.join(dirname, 'x.npy'), x)
    np.save(os.path.join(dirname, 'y.npy'), y)
    np.save(os.path.join(dirname, 'z.npy'), z)
    np.save(os.path.join(dirname, 'vx.npy'), vx)
    np.save(os.path.join(dirname, 'vy.npy'), vy)
    np.save(os.path.join(dirname, 'vz.npy'), vz)
    print 'Saved velocities and positions.'
    
def plot_Sun_and_siblings(N = 50, dirname = './output/', Nbins = 25):
    """
    Plot orbits, final positions and corresponding distributions of the simulation
    of the Sun and siblings in the Galactic potential.
    """
    
    # load positions seperately
    x = np.load(os.path.join(dirname, 'x.npy'))
    y = np.load(os.path.join(dirname, 'y.npy'))
    z = np.load(os.path.join(dirname, 'z.npy'))
    vx = np.load(os.path.join(dirname, 'vx.npy'))
    vy = np.load(os.path.join(dirname, 'vy.npy'))
    vz = np.load(os.path.join(dirname, 'vz.npy'))
    
    # plot paths of all particles
    plt.scatter(x[:,0], y[:,0], color = tableau20[0], label = 'Sun', zorder=1)
    for i in range(100): # plot first 100 positions
        plt.scatter(x[:,i+1], y[:,i+1], color = tableau20[i%20], zorder=0)
    plt.legend(loc=2)
    plt.xlabel('X [pc]', size = 15)
    plt.ylabel('Y [pc]', size = 15)
    plt.savefig(os.path.join(dirname,'orbits.png'), dpi=300)
    plt.savefig(os.path.join(dirname,'orbits.pdf'))
    #plt.show()
    plt.close()
    
    # plot final postions of all particles
    fig, ax = plt.subplots()
    plt.scatter(x[-1,0], y[-1,0], facecolor = tableau20[0], marker = '*', label = 'Sun', s = 250, edgecolor = tableau20[0], zorder =1)
    for i in range(N+1):
        plt.scatter(x[-1,i+1], y[-1,i+1], color = tableau20[i%20], zorder =0)
    plt.legend(loc=2, scatterpoints = 1)
    plt.xlabel('X [pc]', size = 15)
    plt.ylabel('Y [pc]', size = 15)
    circ = matplotlib.patches.Circle(xy = (x[-1,0], y[-1,0]), radius=100, fill = False)
    ax.add_patch(circ)
    plt.savefig(os.path.join(dirname,'finalpositions.png'), dpi=300)
    plt.savefig(os.path.join(dirname,'finalpositions.pdf'))
    #plt.show()
    plt.close()

    # cumulative distribution of the x- y- z- vx- vy- vz
    fig = plt.figure(figsize = (15, 10))
    gs = gridspec.GridSpec(2, 3)
    gs.update(wspace = 0.20, hspace = 0.20, bottom = 0.20)
    plt.suptitle('Normalized distribution', size = 20)
    ax1 = plt.subplot(gs[0:1,0:1])
    ax2 = plt.subplot(gs[0:1,1:2])
    ax3 = plt.subplot(gs[0:1,2:3])
    ax4 = plt.subplot(gs[1:2,0:1])
    ax5 = plt.subplot(gs[1:2,1:2])
    ax6 = plt.subplot(gs[1:2,2:3])
    
    ax1.hist(x[-1,:], Nbins, normed=1, color = tableau20[0])
    ax2.hist(y[-1,:], Nbins, normed=1, color = tableau20[0])
    ax3.hist(z[-1,:], Nbins, normed=1, color = tableau20[0])
    ax4.hist(vx[-1,:], Nbins, normed=1, color = tableau20[2])
    ax5.hist(vy[-1,:], Nbins, normed=1, color = tableau20[2])
    ax6.hist(vz[-1,:], Nbins, normed=1, color = tableau20[2])
    
    ax1.set_title('x [pc]', size = 15)
    ax2.set_title('y [pc]', size = 15)
    ax3.set_title('z [pc]', size = 15)
    ax4.set_title(r'v$_x$ [km/s]', size = 15)
    ax5.set_title(r'v$_y$ [km/s]', size = 15)
    ax6.set_title(r'v$_z$ [km/s]', size = 15)
    plt.savefig(os.path.join(dirname,'distfinalpositions.png'), dpi=300)
    plt.savefig(os.path.join(dirname,'distfinalpositions.pdf'))
    #plt.show()
    plt.close()

def get_cumulative_dist2Sun(N = 50, dirname = './output/'):
    """
    Plot the cumulative distribution of the final distances to the Sun for 
    all star in the cluster.
    """
    
    # Nbins is number of stars
    Nbins = N
    
    # load positions seperately
    x = np.load(os.path.join(dirname, 'x.npy'))
    y = np.load(os.path.join(dirname, 'y.npy'))
    z = np.load(os.path.join(dirname, 'z.npy'))
    vx = np.load(os.path.join(dirname, 'vx.npy'))
    vy = np.load(os.path.join(dirname, 'vy.npy'))
    vz = np.load(os.path.join(dirname, 'vz.npy'))
    
    # get the position/velocity of the Sun at the end of the simulation
    x_Sun = x[-1,0]
    y_Sun = y[-1,0]
    z_Sun = z[-1,0]
    vx_Sun = vx[-1,0]
    vy_Sun = vy[-1,0]
    vz_Sun = vz[-1,0]
    print 'Birthplace of the Sun:(',x_Sun,',',y_Sun,',',z_Sun,')'
    print 'Birthplace velocity of the Sun:(',vx_Sun,',',vy_Sun,',',vz_Sun,')'
    
    # caculate the distance/rel velocity of the final positions with that of the Sun
    distance = np.zeros(N)
    relvelocity = np.zeros(N)
    for i in range(N):
        distance[i] = np.sqrt( (x[-1,i+1]-x_Sun)**2 + (y[-1,i+1]-y_Sun)**2 + (z[-1,i+1]-z_Sun)**2)
        relvelocity[i] = np.sqrt( (vx[-1,i+1]-vx_Sun)**2 + (vy[-1,i+1]-vy_Sun)**2 + (vz[-1,i+1]-vz_Sun)**2)

    # calculate histogram and from histogram the ECDF
    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(distance, Nbins, normed=1, histtype='step', cumulative=True, lw = 1.5, color = tableau20[0])
    plt.xlim(bins[0],bins[-2])
    plt.ylim(0,1)
    plt.axvline(x = 100, lw = 1.5 , ls = '--', color = 'k')
    plt.xlabel('Final distance to Sun [parsec]', size = 15)
    plt.ylabel('Cumulative distribution', size = 15)
    plt.savefig(os.path.join(dirname,'finalcumdist.png'), dpi=300)
    plt.savefig(os.path.join(dirname,'finalcumdist.pdf'))
    #plt.show()
    plt.close()

    # average velocity of objects < 100 parsec
    ind = np.where(distance < 100)[0]
    avg_relvelocity = np.mean(relvelocity[ind])
    print '|Objects < 100 parsec:',len(ind)
    print 'Average velocity of object < 100 parsec:',avg_relvelocity,'km/s'
    
def sim_multiple_cluster_sizes(dirname = './multipleclustersizes'):
    """
    Simulate Sun around Galactic potential in a open cluster for multiple cluster
    sizes.
    """
    
    # run all simulations
    Nstars_values = np.array([50,100,200,500,1000,1500,2000])
    for Nstars in Nstars_values:
        
        #create outputfolder for this star
        outputfolder = os.path.join(dirname, str(Nstars))
        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)
            
        # run simulation
        sim_Sun_and_siblings(N = Nstars, dirname = outputfolder, pointparticles = False)
        
        # plot the results of this simulation
        plot_Sun_and_siblings(dirname = outputfolder, N = Nstars)
        get_cumulative_dist2Sun(dirname = outputfolder, N = Nstars)
        
def plot_multiple_cluster_sizes(dirname = './multipleclustersizes'):
    """
    Plot results of simulations of multiple cluster sizes.
    """
    
    # run all simulations
    Nstars_values = np.array([50,100,200,500,1000,1500,2000])
    N100pc = np.zeros(len(Nstars_values))
    for n, Nstars in enumerate(Nstars_values):

        # get the outputfolder
        outputfolder = os.path.join(dirname, str(Nstars))      
        
        # load positions seperately
        x = np.load(os.path.join(outputfolder, 'x.npy'))
        y = np.load(os.path.join(outputfolder, 'y.npy'))
        z = np.load(os.path.join(outputfolder, 'z.npy'))
        
        # get the position/velocity of the Sun at the end of the simulation
        x_Sun = x[-1,0]
        y_Sun = y[-1,0]
        z_Sun = z[-1,0]

        # caculate the distance/rel velocity of the final positions with that of the Sun
        distance = np.zeros(Nstars)
        for i in range(Nstars):
            distance[i] = np.sqrt( (x[-1,i+1]-x_Sun)**2 + (y[-1,i+1]-y_Sun)**2 + (z[-1,i+1]-z_Sun)**2)
            
        # objects < 100 parsec
        ind = np.where(distance < 100)[0]
        N100pc[n] = len(ind)
    
    # linear fit to the line
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Nstars_values, N100pc)
    yfit = slope * Nstars_values + intercept    
    
    # create a plot
    plt.annotate('a = '+str(np.round(slope,2)), xy = (1650, 4.5), size = 12)
    plt.annotate('b = '+str(np.round(intercept,2)), xy = (1650, 3), size = 12)
    plt.scatter(Nstars_values, N100pc, color = tableau20[0], zorder = 1, label = 'Simulations')
    plt.plot(Nstars_values, yfit, color = tableau20[2], zorder = 0, label = 'Best fit')
    plt.xlabel('Number of other stars in cluster', size = 15)
    plt.ylabel('Number of stars < 100 parsec', size = 15)
    plt.xlim(0,2050)
    plt.legend(loc=2, scatterpoints = 1)
    plt.savefig(os.path.join(dirname,'multicluster.png'), dpi=300)
    plt.savefig(os.path.join(dirname,'multicluster.pdf'))
    #plt.show()
    plt.close()
  
def main():
    """
    CAP - Assignment 5.
    """

    # part 1
    #sim_Sun_and_siblings(pointparticles = True)
    #plot_Sun_and_siblings(dirname = './pointparticles')

    # part 2
    #sim_Sun_and_siblings(pointparticles = False, dirname = './opencluster')
    #plot_Sun_and_siblings(dirname = './opencluster_oursolarbirthplace', N = 50)
    get_cumulative_dist2Sun(dirname = './opencluster_oursolarbirthplace', N = 50)

    # part 3
    #sim_multiple_cluster_sizes()
    #plot_multiple_cluster_sizes(dirname = './multipleclustersizes/')
    #get_cumulative_dist2Sun(dirname = './multipleclustersizes/50', N = 50)
    #get_cumulative_dist2Sun(dirname = './multipleclustersizes/100', N = 100)
    #plot_Sun_and_siblings(dirname = './multipleclustersizes/500', N = 500)
    #get_cumulative_dist2Sun(dirname = './multipleclustersizes/500', N = 500)
    #plot_Sun_and_siblings(dirname = './multipleclustersizes/1000', N = 1000)
    #get_cumulative_dist2Sun(dirname = './multipleclustersizes/1000', N = 1000)
    #plot_Sun_and_siblings(dirname = './multipleclustersizes/1500', N = 1500)
    #get_cumulative_dist2Sun(dirname = './multipleclustersizes/1500', N = 1500)
    #plot_multiple_cluster_sizes(dirname = './multipleclustersizes/')

main()
