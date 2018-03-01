import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches
import matplotlib.animation as animation


def update(t):
    # Optionally. Clear axes and reset limits
    plt.gca().cla() 
    ax.set_xlim(-10000,10000)
    ax.set_ylim(-10000,10000)
    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    
    if t <101:
	
	for i in range(1, n):
	    ax.scatter(x[t, i], y[t, i], marker="o", c='k', s=10)
	ax.scatter(x[t, i], y[t, i], marker="*", c='y', s=40, label='Sun')
	title = str(round(4.56/frames * t, 2)) + ' Gyr'
    
    else:
	for i in range(1, n):
	    ax.scatter(x[-1, i], y[-1, i], marker="o", c='k', s=10)
	ax.scatter(x[-1, i], y[-1, i], marker="*", c='y', s=40, label='Sun')
	title = str(round(4.56/frames * 121, 2)) + ' Gyr'
   
    ax.set_title(title)
    plt.legend(loc=2, scatterpoints = 1)
    plt.tight_layout()





if __name__ in ('__main__', '__plot__'):
    
    x = np.load('./multipleclustersizes/100/x.npy')
    y = np.load('./multipleclustersizes/100/y.npy')
    
    frames = 121
    n, steps = x.shape
       
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)


    fig, ax = plt.subplots()
    ani = animation.FuncAnimation(fig, update, frames=frames, interval=200, repeat=False)
    ani.save('movie.mp4', writer=writer)
    plt.show()

    print '\nDone!!'
    

