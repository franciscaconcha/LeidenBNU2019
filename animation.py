import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani


def init():

    # Initializes plot
    # Viewing window is a box centered on the origin, such that all points are
    # within it
    ax_ani.set_xlim(-ext_ani, ext_ani)
    ax_ani.set_ylim(-ext_ani, ext_ani)
    return l1,

def update(frame):

    # Updates plot
    xdata = frame[0]
    ydata = frame[1]
    l1.set_data(xdata, ydata)
    return l1,


global ax_ani
global l1
global ext_ani

moviename = 'test'

animation_time = 10. # duration of movie in seconds
N = 1000 # number of points in time

# positions of 2 particles in 3 dimensions, at N moments in time
R = np.zeros((N, 3, 2))

# Particle 1 moves around; particle 2 isn't assigned any position, so remains
# at origin
t = np.linspace(0., 100., N)
R[:,0,0] = t
R[:,1,0] = np.cos(t)
R[:,2,0] = np.sin(t)


plot_dimensions = [0,1] # Dimensions to plot; 0 means X, 1 means Y, 2 means Z

SKIP = 1    # How many data points to skip between frames


fig = plt.figure(1)
ax_ani = fig.add_subplot(111)

l1, = ax_ani.plot([], [], 'ro', animated=True)

ext_ani = 1.2*np.max(np.abs(R[:,plot_dimensions])) # Defines viewing window

ANI = ani.FuncAnimation(fig, update, interval=1000.*animation_time/float(N),
                        frames=R[::SKIP,plot_dimensions], init_func=init,
                        blit=True)

writer = ani.FFMpegWriter(bitrate=1800, fps=15)
ANI.save(moviename+'.mp4', writer=writer)
