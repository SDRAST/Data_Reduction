import logging
import numpy
import signal
import socket
import sys
import time

from multiprocessing import Process, Queue
from struct import unpack_from

logger = logging.getLogger(__name__)
pkt_size = 1026*64/8
pkt_fmt = "!8cHHI"+1024*4*"H"
      
def signal_handler(signal, frame):
  """
  This does not end the thread
  """  
  print 'You pressed Ctrl+C!'
  reader.terminate()
  for count in range(10):
    worker[count].terminate()
  processor.terminate()
  print "Total",time.time()-begin,"s"
  sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

def _open_socket(host='10.0.0.2', port=60000):
    """
    opens socket to ROACH 10Gbe
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.bind((host,port))
    logger.info("_open_socket: socket bound")
    return s

def capture(socket):
    while True:
      try:
        data, addr = socket.recvfrom(pkt_size)
        data_queue.put(data)
      except (KeyboardInterrupt):
        print "Got <CTRL-C>"
        reader.terminate()
      
def get_packet(data_queue, packet_queue):
    """
    gets unscrambled packet
    """
    global packets
    
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

    while True:
      try:
        pkt_buf = data_queue.get()
        result = unpack_from(pkt_fmt, pkt_buf)
        pkt = {}
        pkt['hdr'] = result[:8]
        pkt['pkt cnt sec'] = result[8]
        pkt['sec cnt'] = result[9]
        pkt['raw pkt cnt'] = result[10]
        data = result[11:]
        pkt['pwr'], pkt['krt'] = unscramble(data)
        logger.debug("get_packet: unscrambled packet %d", pkt['raw pkt cnt'])
        packet_queue.put(pkt)
      except (KeyboardInterrupt):
        # wait for reader to finish
        pass

def resequence(packet_queue, ordered_queue):
  """
  processes packets
  """
  packets = []
  sort_on = 'raw pkt cnt'
  while True:
    try:
      while len(packets) < 10000 + 100:
        packets.append(packet_queue.get())
      logger.info("resequence: got 100100 packets")
      decorated = [(dict_[sort_on], dict_) for dict_ in packets]
      decorated.sort()
      result = [dict_ for (key, dict_) in decorated]
      ordered_queue.put(result[:10000])
      packets = result[10000:]
    except KeyboardInterrupt:
      pass
  
logging.basicConfig()
mylogger = logging.getLogger()
mylogger.setLevel(logging.INFO)

socket = _open_socket()
data_queue = Queue()
packet_queue = Queue()
ordered_queue = Queue()

reader = Process(target=capture, args=(socket,))
processor = Process(target=resequence, args=(packet_queue, ordered_queue))
worker = {}
for count in range(10):
    worker[count] = Process(target=get_packet, args=(data_queue, packet_queue))
    
begin = time.time()
mylogger.debug("Starting at %s", begin)
processor.start()
for count in range(10):
    worker[count].start()
    mylogger.debug("started %d", count)
    time.sleep(0.1)
reader.start() # get all the others going before starting the reader

print "Total",time.time()-begin,"s"

one_second = []
try:
  one_second.append(ordered_queue.get())
  print "got one_second"
except KeyboardInterrupt:
  pass
  
reader.join()
for count in range(10):
    worker[count].join()
processor.join()
  
  
####################################### this does not work
      
class Worker(Process):
  """
  """
  def __init__(self, queue):
    """
    """
    Process.__init__(self, target=self.get_packet)
    self.logger = logging.getLogger(logger.name+"."+self.name)
    self.packet = {}
    self.data_queue = queue
    print self.name, "started"
      
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
    packets.append(pkt)
  


