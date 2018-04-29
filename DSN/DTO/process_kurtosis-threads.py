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
#import Queue
import signal
import socket
import sys
#import threading
import time

from multiprocessing import Process, Queue
from struct import calcsize, unpack_from

pkt_fmt = "!8cHHi"+1024*4*"H"

from MonitorControl import DeviceReadThread
from support import NamedClass

logger = logging.getLogger(__name__)

npkts=10000

class Worker(Process):
  """
  """
  def __init__(self, queue):
    """
    """
    Process.__init__(self, target=self.run)
    self.logger = logging.getLogger(logger.name+"."+self.name)
    self.durations = []
    self.pkt_cnt_sec = []
    self.raw_pkt_cnt = []
    self.packet = {}
    self.data_queue = queue
    print self.name, "started"
    
  def run(self):
    """
    """
    self.startime = time.time()
    print "at", self.startime
    for index in range(npkts):
      self.packet[index] = self.get_packet()
    duration = time.time()-self.startime
    self.durations.append(duration)
    self.pkt_cnt_sec.append([self.packet[0]['pkt cnt sec'],
                             self.packet[npkts-1]['pkt cnt sec']])
    self.raw_pkt_cnt.append([self.packet[0]['raw pkt cnt'],
                             self.packet[npkts-1]['raw pkt cnt']])
      
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
      power['I'] = numpy.append(D[:512,0],D[:512,1])
      power['Q'] = numpy.append(D[512:,0],D[512:,1])
      kurt = {}
      kurt['I'] = numpy.append(D[:512,2],D[:512,3])/4096.
      kurt['Q'] = numpy.append(D[512:,2],D[512:,3])/4096.
      return power, kurt

    pkt_buf = self.data_queue.get()
    result = unpack_from(pkt_fmt, pkt_buf)
    pkt = {}
    pkt['hdr'] = result[:8]
    self.logger.debug("get_packet: header is %s", pkt['hdr'])
    pkt['pkt cnt sec'] = result[8]
    pkt['sec cnt'] = result[9]
    pkt['raw pkt cnt'] = result[10]
    data = result[11:]
    pkt['pwr'], pkt['krt'] = unscramble(data)
    self.logger.debug("get_packet: unscrambled")
    return pkt

  def terminate(self):
    """
    Thread termination routine
    """
    self.end_flag = True
    
      
def signal_handler(signal, frame):
  """
  This does not end the thread
  """
  global kp
  print 'You pressed Ctrl+C!'
  kp.quit()
  sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)


class Receiver(NamedClass):
  """
  Captures packets from ROACH over 10Gbe port and operates on them.
  
  public attributes::
    data_queue - Queue.Queue object where the incoming data are put
    logger     - logging.Logger object
    pkt_size   - (should probably be a class variable)
    reader     - DeviceReadThread obj. reads the data and puts in on the queue
    socket     - connection to ROACH 10Gbe
  """
  def get_data(self):
    """
    used by DeviceReadThread object to handle the data
    """
    self.logger.debug("get_data: started")
    while self.running:
      try:
        self.logger.debug("Waiting...")
        data, addr = self.socket.recvfrom(self.pkt_size)
        self.logger.debug('Packet len: %d', len(data))
        self.data_queue.put(data)
      except (socket.timeout):
        self.logger.debug('Socket timeout...')
      except (KeyboardInterrupt):
        self.logger.info("action: Ctrl-C noted")
        self.quit()
        
  def __init__(self):
    """
    initializes Receiver
    
    creates the reader and the queue
    """
    self.logger = logging.getLogger(logger.name+".Receiver")
    if socket.gethostname() == 'gpu2':
      self.socket = self._open_socket(host='10.0.0.2')
    else:
      self.socket = self._open_socket()
    self.pkt_size = 1026*64/8
    self.data_queue = Queue()
    #self.reader = DeviceReadThread(self, self.socket)
    self.running = True
    self.reader = Process(target=self.run)
    self.logger.debug("__init__: reader defined")
    self.reader.start()
    self.reader.join()
    #self.reader.daemon = True # these don't seem to do any good for catching
    print "Receiver started"
    
  def _open_socket(self, host='10.0.0.12', port=60000):
    """
    opens socket to ROACH 10Gbe
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.bind((host,port))
    self.logger.info("_open_socket: socket bound")
    return s
  
  def run(self):
    """
    used by DeviceReadThread object to handle the data
    """
    self.logger.debug("get_data: started")
    while self.running:
      try:
        self.logger.debug("Waiting...")
        data, addr = self.socket.recvfrom(self.pkt_size)
        self.logger.debug('Packet len: %d', len(data))
        self.data_queue.put(data)
      except (socket.timeout):
        self.logger.debug('Socket timeout...')
      except (KeyboardInterrupt):
        self.logger.info("action: Ctrl-C noted")
        self.quit()
  
  def quit(self):
    """
    terminates the Receiver object
    """
    self.logger.warning("Quitting")
    self.reader.terminate()
       
if __name__ == "__main__":
  logging.basicConfig()
  mylogger = logging.getLogger()
  mylogger.setLevel(logging.INFO)
  
  kp = Receiver()
  print "kp started"
  begin = time.time()
  mylogger.debug("Starting at %s", begin)
  wk = {}
  done = {}
  for count in range(10):
    wk[count] = Worker(queue = kp.data_queue)
    wk[count].start()
    wk[count].join()
  print 'pckts packets duration'
  for count in range(10):
    print "checking", count
    done[count] = False
    while done[count] == False:
      if len(wk[count].packet) == npkts:
        done[count] = True
        print count, "is done"
    if wk[count].pkt_cnt_sec[1] < wk[count].pkt_cnt_sec[0]:
      wk[count].pkt_cnt_sec[1] += 2**15
    if wk[count].raw_pkt_cnt[1] < wk[count].raw_pkt_cnt[0]:
      wk[count].raw_pkt_cnt[1] += 2**15
    print("%5d %5d %8.5f" % (wk[count].pkt_cnt_sec[1]-wk[count].pkt_cnt_sec[0],
                             wk[count].raw_pkt_cnt[1]-wk[count].raw_pkt_cnt[0],
                             wk[count].durations))
  print "Total",time.time()-begin,"s"
  kp.quit()
  for count in range(10):
    wk[count].terminate()

   
