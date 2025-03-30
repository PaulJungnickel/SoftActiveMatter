import numpy as np
from pylab import *
from matplotlib import pyplot as plt
from matplotlib import animation

trajectory = loadtxt('rousecoords2d_sample.txt')

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, 20), ylim=(0,15))
line, = ax.plot([], [], "-o", color="b", markevery=1)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line, 

# animation function.  This is called sequentially
def animate(i):
    x = (trajectory[i*6,:]+2.5)
    y = (trajectory[i*6+1,:]+7.5)
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=1500, interval=10, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html

anim.save('polymer.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

# plt.show()
