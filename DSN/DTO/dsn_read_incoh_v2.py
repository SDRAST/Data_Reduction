"""
Process packets from kurtosis spectrometer firmware
"""
import numpy
import struct

from pylab import *

filepath = '/home/kuiper/mnt/dto/home/kuiper/mnt/gpu1/var/tmp/'
filename = '20160705194100_1'
fd = open(filepath+filename,'r')

end_packet = 2
pow_rs = {}
kurt_rs = {}
for j in range(end_packet):
    # discard first two frames
    nbytes = struct.calcsize("q") # uint64
    discard = fd.read(nbytes) # read one uint64 to discard header frame
    nbytes = struct.calcsize("q") # ubit64 is same as uint64
    discard = fd.read(nbytes) # read one ubit64 to discard counter frame
    # get the remaining 1024 frames
    pktfmt = 2*4*'B'*512*2 # two pols x 4-byte ints x 512 words x (pwr,kurt)
    nbytes = struct.calcsize(pktfmt)
    data_all = numpy.array(struct.unpack_from(pktfmt, fd.read(nbytes)))
    # note:  this agrees with what the Matlab program reads
    # In [1]: data_all[:10]
    # Out[1]: array([  0, 217,   0,  56,  10, 120,  41, 122,   0,  16])
    # In [2]: data_all[-10:]
    # Out[2]: array([ 44,  15,   0,   8,   0,   6,  42,  18,  39, 107])

    # The data are unpacked into a multi-dimensional array
    data_all_rs = data_all.reshape((2,512,4,2))
    # In [12]: for i in range(5):
    #              print i,data_all_rs[0,i,:,:].transpose()
    #    ....: 
    #   0 [[  0   0  10  41] [217  56 120 122]]
    #   1 [[  0   0  61  41] [ 16  23 188 117]]
    #   2 [[  0   0  71  51] [ 13  32  50 202]]
    #   3 [[  0   0  67  42] [ 16  30 235 243]]
    #   4 [[  0   0  60  52] [ 15  23 227 239]]
    #   5 [[  0   0  68  49] [ 19  26 160 211]]
    #   6 [[  0   0  72  57] [ 19  23 252 135]]
    # ...
    # 510 [[  0   0  42  36] [ 41  14  37 207]]
    # 511 [[  0   0  46  42] [  9   6   4 137]]
    #   0 [[  0   0   9  45] [221  56 172 153]]
    #   1 [[  0   0  60  38] [ 14  27 116 108]]
    # ...
    # 508 [[  0   0  43  52] [ 45   7 121 237]]
    # 509 [[  0   0  50  43] [ 48   8  25  33]]
    # 510 [[  0   0  47  44] [ 36  15 248  15]]
    # 511 [[  0   0  42  39] [  8   6  18 107]]

    # In python the following collapses the fourth dimension (axis 3)
    data_all_packed = data_all_rs[:,:,:,1] + data_all_rs[:,:,:,0]*2**8
    # In [3]: data_all_packed.shape
    # Out[3]: (2, 512, 4)
    # In [21]: for i in range(5):
    #              print i,data_all_packed[0,i,:]
    #    ....: 
    #   0 [  217    56  2680 10618]
    #   1 [   16    23 15804 10613]
    #   2 [   13    32 18226 13258]
    #   3 [   16    30 17387 10995]
    #   4 [   15    23 15587 13551]
    # ...
    # 508 [   50    10  9930 10906]
    # 509 [   47    10 12108 12022]
    # 510 [   41    14 10789  9423]
    # 511 [    9     6 11780 10889]
    # In [4]: for i in range(5):
    #             print i,data_all_packed[1,i,:]
    #    ...: 
    #   0 [  221    56  2476 11673]
    #   1 [   14    27 15476  9836]
    #   2 [   12    26 12248 10065]
    #   3 [   15    21 14595 10727]
    #   4 [   14    23 12256 10888]
    # ...
    # 508 [   45     7 11129 13549]
    # 509 [   48     8 12825 11041]
    # 510 [   36    15 12280 11279]
    # 511 [    8     6 10770 10091]

    # The power and kurtosis are extracted from the array
    pow_all = data_all_packed[:,:,0:2]
    kurt_all = data_all_packed[:,:,2:4]
    # In [2]: pow_all.shape
    # Out[2]: (2, 512, 2)
    
    # In [6]: pow_all[0,:,:].transpose()
    # Out[6]: 
    # array([[217,  16,  13, ...,  47,  41,   9],
    #        [ 56,  23,  32, ...,  10,  14,   6]])
    # In [5]: pow_all[1,:,:].transpose()
    # Out[5]: 
    # array([[221,  14,  12, ...,  48,  36,   8],
    #        [ 56,  27,  26, ...,   8,  15,   6]])

    # kurt_all = data_all_permute1[2:4,:,:]
    # In [7]: kurt_all[0,:,:].transpose()
    # Out[7]: 
    # array([[ 2680, 15804, 18226, ..., 12108, 10789, 11780],
    #        [10618, 10613, 13258, ..., 12022,  9423, 10889]])
    # In [8]: kurt_all[1,:,:].transpose()
    # Out[8]: 
    # array([[ 2476, 15476, 12248, ..., 12825, 12280, 10770],
    #        [11673,  9836, 10065, ..., 11041, 11279, 10091]])
    
    # Now we need to combine the upper and lower halves of the spectra
    # The first index selects the upper and lower block and the last index
    # the left and right columns
    pow_rs[j] = numpy.append(pow_all[j,:,0], pow_all[j,:,1])
    kurt_rs[j] = numpy.append(kurt_all[j,:,0], kurt_all[j,:,1])
    
figure()
for pol in [0,1]:
    subplot(2,1,2)
    plot(abs(kurt_rs[pol])/2.**12, label="Pol "+str(pol))
    subplot(2,1,1)
    semilogy(abs(pow_rs[pol]), label="Pol "+str(pol)) 
grid()
legend()          
fd.close()
show()


