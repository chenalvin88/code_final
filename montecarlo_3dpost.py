import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm as tqdm

def random_three_vector():
    """
    Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    :return:
    """
    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = np.sin( theta) * np.cos( phi )
    y = np.sin( theta) * np.sin( phi )
    z = np.cos( theta )
    return [x,y,z]

# threetups1,threetups2 = [],[]
angles = []
for _ in tqdm(range(int(1e6))):
    v1 = random_three_vector()
    v2 = random_three_vector()
    # threetups1.append(v1)
    # threetups2.append(v2)
    angle = np.degrees(np.arccos(np.dot(v1,v2)))
    if angle<=90:
      angles.append(angle)

# threetups1 = list(map(list, zip(*threetups1)))
# plt.gca(projection='3d').scatter(threetups1[0],threetups1[1],threetups1[2])
# print(threetups1[0])
# print(angles)
# plt.hist(angles,bins=90)
# plt.xlim(0,90)

hist,deges,plot = plt.hist(angles,bins=90)
hist = hist/hist[-1]/(180/np.pi)
print(hist)
plt.close()
x = np.arange(0,90)
plt.bar(x,hist,width = 1)
plt.plot(x,np.sin(np.radians(x))/(180/np.pi),c='black')
plt.show()