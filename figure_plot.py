import numpy as np
from matplotlib import pyplot as plt
from pylab import *
img_array = np.load('mss.npy')
plt.imshow(img_array,'gray')
plt.show()