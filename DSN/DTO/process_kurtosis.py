"""
This grabs data from the roach ten times.  Each time it reads 10,000 packets
and stores them in an array called 'packet'. It saves pkt_cnt_sec from
the first and the last packet in a list 'pkt_cnt_sec'.  Likewise it stores
raw_pkt_cnt from the first and last record in 'raw_pkt_cnt'.  It then prints
the differences for each of the ten times around the loop.


Packet Structure::
  __________________________________________________________________
  |                         64-bit Frames                          |
  |________________________________________________________________|
  | Frame |          uint(3)          | uint(2) | uint(1)| uint(0) |
  |_______|___________________________|_________|__________________|
  |   0   |                user defined header                     |
  |   1   | (pol << 63) + pkt_cnt_sec | sec_cnt |   raw_pkt_cnt    |
  |   2   |            F512           |    F0   |  P512  |   P0    |
  |   3   |            F513           |    F1   |  P513  |   P1    |
  |  ...  |             ...           |   ...   |   ...  |  ...    |
  |  512  |           F1022           |  F510   | P1022  | P510    |
  |  513  |           F1023           |  F511   | P1023  | P511    |
  |----------------------------------------------------------------||
  |  514  |            F512           |    F0   |  P512  |   P0    |
  |  515  |            F513           |    F1   |  P513  |   P1    |
  |  ...  |             ...           |   ...   |   ...  |  ...    |
  | 1024  |           F1022           |  F510   | P1022  | P510    |
  | 1025  |           F1023           |  F511   | P1023  | P511    |
  |----------------------------------------------------------------|
where P means 'power' and F means 'fourth moment'.

So the packets are 1026*64/8 = 8028 bytes long. 
The unpacking should be "8cHHi"+1024*4*"H"
"""
import logging
import numpy
import Queue
import signal
import socket
import sys
import time

from struct import calcsize, unpack_from

pkt_fmt = "!8cHHi"+1024*4*"H"

from MonitorControl import DeviceReadThread
from support import NamedClass

logger = logging.getLogger(__name__)

def signal_handler(signal, frame):
  """
  This does not end the thread
  """
  print 'You pressed Ctrl+C!'
  sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

class KurtosisProcessor(NamedClass):
  """
  Captures packets from ROACH over 10Gbe port and operates on them.
  
  public attributes::
    data_queue - Queue.Queue object where the incoming data are put
    logger     - logging.Logger object
    pkt_size   - (should probably be a class variable)
    reader     - DeviceReadThread obj. reads the data and puts in on the queue
    socket     - connection to ROACH 10Gbe
  """
  def __init__(self):
    """
    initializes KurtosisProcessor
    
    creates the reader and the queue
    """
    self.logger = logging.getLogger(logger.name+".KurtosisProcessor")
    if socket.gethostname() == 'gpu2':
      self.socket = self._open_socket(host='10.0.0.2')
    else:
      self.socket = self._open_socket()
    self.pkt_size = 1026*64/8
    self.data_queue = Queue.Queue()
    self.reader = DeviceReadThread(self, self.socket)
    self.reader.daemon = True # these don't seem to do any good for catching
    if self.reader.isAlive(): # keyboard interrupts
      self.reader.join(1)     #
    
  def _open_socket(self, host='10.0.0.12', port=60000):
    """
    opens socket to ROACH 10Gbe
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.bind((host,port))
    self.logger.info("_open_socket: socket bound")
    return s
  
  def action(self, s):
    """
    used by DeviceReadThread object to handle the data
    """
    try:
      self.logger.debug("Waiting...")
      data, addr = s.recvfrom(self.pkt_size)
      self.logger.debug('Packet len: %d', len(data))
      self.data_queue.put(data)
    except (socket.timeout):
      self.logger.debug('Socket timeout...')
    except (KeyboardInterrupt):
      self.logger.info("action: Ctrl-C noted")
      self.quit()
      
  def get_packet(self):
    """
    gets unscrambled packet
    """
    def unscramble(data):
      """
      unscrambles a packet
      """
      D = numpy.array(data).reshape((1024,4))
      power = {}
      power['I'] = numpy.append(D[:512,3],D[:512,2])
      power['Q'] = numpy.append(D[512:,3],D[512:,2])
      kurt = {}
      kurt['I'] = numpy.append(D[:512,1],D[:512,0])
      kurt['Q'] = numpy.append(D[512:,1],D[512:,0])
      return power, kurt

    pkt_buf = self.data_queue.get()
    result = unpack_from(pkt_fmt, pkt_buf)
    pkt = {}
    pkt['hdr'] = result[:8]
    logging.debug("get_data.get_packet: header is %s", pkt['hdr'])
    pkt['pkt cnt sec'] = result[8]
    pkt['sec cnt'] = result[9]
    pkt['raw pkt cnt'] = result[10]
    logging.debug(
             "get_data.get_packet: pkt_cnt_sec=%u, sec_cnt=%u, raw_pkt_cnt=%u",
              pkt['pkt cnt sec'], pkt['sec cnt'], pkt['raw pkt cnt'])
    data = result[11:]
    pkt['pwr'], pkt['krt'] = unscramble(data)
    return pkt
  
  def quit(self):
    """
    terminates the KurtosisProcessor object
    """
    self.logger.warning("Quitting")
    self.reader.terminate()
       
if __name__ == "__main__":
  kp = KurtosisProcessor()
  kp.reader.start()
  durations = []
  pkt_cnt_sec = []
  raw_pkt_cnt = []
  for count in range(10):
    packet = {}
    npkts=10000
    startime = time.time()
    for index in range(npkts):
      packet[index] = kp.get_packet()
    duration = time.time()-startime
    durations.append(duration)
    pkt_cnt_sec.append([packet[0]['pkt cnt sec'], packet[npkts-1]['pkt cnt sec']])
    raw_pkt_cnt.append([packet[0]['raw pkt cnt'], packet[npkts-1]['raw pkt cnt']])
  kp.quit()
  print "Times to acquire 10,000 packets:"
  print durations
  print "   dif pkt_cnt_sec     dif raw_pkt_cnt"
  for count in range(10):
    if pkt_cnt_sec[count][1] < pkt_cnt_sec[count][0]:
      pkt_cnt_sec[count][1] += 2**15
    if raw_pkt_cnt[count][1] < raw_pkt_cnt[count][0]:
      raw_pkt_cnt[count][1] += 2**15
    print count, pkt_cnt_sec[count][1]-pkt_cnt_sec[count][0], raw_pkt_cnt[count][1]-raw_pkt_cnt[count][0]
   
