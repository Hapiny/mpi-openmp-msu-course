import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


proc_list = [1,2,4,6,8,10,12,14,16]
size_list = [100,500,1000,2000,4000,6000,8000,10000,12000,14000,16000]
x = []
y = []
z = []
dx = []
dy = []
dz = []
colours = []

for proc in proc_list:
	for size in size_list:
		time = 0
		with open('omp-%d-%d.out'%(size, proc), 'r') as file:
			time = file.read().strip().split(' ')[-1]
			time = float(time)
			if proc == 1:
				colours.append('red')
			else:	
				colours.append('green')
				
		y.append(proc-0.5)
		x.append(size-0.5)
		z.append(0)

		dx.append(800)
		dy.append(1)
		dz.append(time)

ax3d = plt.figure(figsize=(8,8)).gca(projection='3d')
ax3d.bar3d(x, y, z, dx, dy, dz, color=colours)

ax3d.set_ylabel('Threads number')
ax3d.set_xlabel('Matrix size')
ax3d.set_zlabel('Time (in seconds)')

ax3d.invert_xaxis()
# ax3d.invert_yaxis()

plt.show()









