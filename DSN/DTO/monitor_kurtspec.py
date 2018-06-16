import cPickle
import numpy

from pylab import *

filename = "/data/HDF5/dss14/2018/167/mon-2018-167-183301.bin"

dfile = open(filename, "rb+")

data = cPickle.load(dfile)
data_array = data.reshape(data.shape + (1,))
count = 0
reading = True
while reading:
  try:
    data = cPickle.load(dfile)
    data_array = numpy.append(data_array, data.reshape(data.shape + (1,)), axis=2)
    count += 1
    print "\r"+str(count),
  except EOFError:
    reading = False
dfile.close()

pwrIimg = numpy.log10(data_array[:,0,:])
pwrQimg = numpy.log10(data_array[:,2,:])
pwrIspec = pwrIimg.mean(axis=1)
pwrQspec = pwrQimg.mean(axis=1)

krtIimg = data_array[:,1,:]
krtQimg = data_array[:,3,:]

figure()
plot(pwrIspec, label="I")
plot(pwrQspec, label="Q")
legend()
grid()
title("Average Passband")

figure(figsize=(30,6))
imshow(pwrIimg)
title("power I")
colorbar()

figure(figsize=(30,6))
imshow(pwrQimg)
title("power Q")
colorbar()

figure(figsize=(30,6))
imshow(krtIimg)
title("kurtosis I")
colorbar()

figure(figsize=(30,6))
imshow(krtQimg)
title("kurtosis Q")
colorbar()

