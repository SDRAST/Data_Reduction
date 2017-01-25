"""
Unpacks packets sent by kurt_spec firmware

Packet Structure::
  __________________________________________________________________
  |                         64-bit Frames                          |
  |________________________________________________________________|
  | Frame |          uint(3)          | uint(2) | uint(1)| uint(0) |
  |_______|___________________________|_________|__________________|
  |   0   |                user defined header                     |
  |   1   | (pol << 63) + pkt_cnt_sec | sec_cnt |   raw_pkt_cnt    |
  |----------------------------------------------------------------|
  |   2   |            F512           |    F0   |  P512  |   P0    |
  |   3   |            F513           |    F1   |  P513  |   P1    |
  |  ...  |             ...           |   ...   |   ...  |  ...    |
  |  512  |           F1022           |  F510   | P1022  | P510    |
  |  513  |           F1023           |  F511   | P1023  | P511    |
  |----------------------------------------------------------------|
  |  514  |            F512           |    F0   |  P512  |   P0    |
  |  515  |            F513           |    F1   |  P513  |   P1    |
  |  ...  |             ...           |   ...   |   ...  |  ...    |
  | 1024  |           F1022           |  F510   | P1022  | P510    |
  | 1025  |           F1023           |  F511   | P1023  | P511    |
  |----------------------------------------------------------------|
where P means 'power' and F means 'fourth moment'.

So the packets are 1026*64= 65664 bits or 65664/8 = 8028 bytes long.
The header packets have format "8cHHi" which is 16 bytes.  The data
comprise calcsize(1024*4*"H") = 8192 bytes.  The unpacking should be
"8cHHi"+1024*4*"H"
"""
import logging
import numpy

from os.path import getsize
from struct import calcsize, unpack_from

pkt_size = 1026*64/8
pkt_fmt = "!8cHHi"+1024*4*"H"

logging.basicConfig(level=logging.DEBUG)
mylogger = logging.getLogger()
mylogger.setLevel(logging.DEBUG)

def get_data(filename):
  """
  """
  def get_packet(pkt_buf):
    """
    """
    def unscramble(data):
      """
      """
      D = numpy.array(data).reshape((1024,4))
      power = {}
      power['I'] = numpy.append(D[:512,3],D[:512,2])
      power['Q'] = numpy.append(D[512:,3],D[512:,2])
      kurt = {}
      kurt['I'] = numpy.append(D[:512,1],D[:512,0])
      kurt['Q'] = numpy.append(D[512:,1],D[512:,0])
      return power, kurt

    result = unpack_from(pkt_fmt, pkt_buf)
    hdr = result[:8]
    logging.debug("get_data.get_packet: header is %s", hdr)
    pkt_cnt_sec = result[8]
    sec_cnt = result[9]
    raw_pkt_cnt = result[10]
    logging.debug(
             "get_data.get_packet: pkt_cnt_sec=%u, sec_cnt=%u, raw_pkt_cnt=%u",
             pkt_cnt_sec, pkt_cnt_sec, raw_pkt_cnt)
    data = result[11:]
    pwr, krt = unscramble(data)
    return (hdr, pkt_cnt_sec, sec_cnt, raw_pkt_cnt, pwr, krt)
    
  datafile = open(filename,'rb')
  databuf = datafile.read()
  datafile.close()
  filesize = getsize(filename)
  num_pkts = filesize/pkt_size
  mylogger.debug("get_data: got data for %d packets", num_pkts)
  pkt_data = {}
  for pkt in range(num_pkts):
    mylogger.debug("get_data: getting packet %d", pkt)
    pkt_data[pkt] = get_packet(databuf[pkt_size*pkt:pkt_size*(pkt+1)])
  return pkt_data
  
if __name__ == "__main__":
  filename = "/var/tmp/20160705194100_1"
  data = get_data(filename)

