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
if socket.gethostname() == 'gpu1':
  IP = '10.0.0.12'
else:
  IP = '10.0.0.2'
  
def signal_handler(signal, frame):
  """
  This does not end the thread
  """  
  print 'You pressed Ctrl+C!'
  socket.close()
  reader.terminate()
  for count in range(10):
    worker[count].terminate()
  print "Total",time.time()-begin,"s"
  sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

def _open_socket(host=IP, port=60000):
    """
    opens socket to ROACH 10Gbe
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.bind((host,port))
    logger.info("_open_socket: socket bound")
    return s

def get_packet(socket, unpacker_queue):
  """
  gets packets and assigns them to unscramblers
  """
  while True:
    try:
      for count in range(10000):
        data, addr = socket.recvfrom(pkt_size)
        unpacker_queue.put(data)
      self.logger.debug("get_packet: for %s", unpacker_queue)
    except KeyboardInterrupt:
      # nothing to do; signal_handler takes care of it
      pass
      
def unscramble_packet(unpacker_queue, ordered_queue):
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

    while True:
      try:
        pkt_buf = unpacker_queue.get()
        result = unpack_from(pkt_fmt, pkt_buf)
        pkt = {}
        pkt['hdr'] = result[:8]
        pkt['pkt cnt sec'] = result[8]
        pkt['sec cnt'] = result[9]
        pkt['raw pkt cnt'] = result[10]
        data = result[11:]
        pkt['pwr'], pkt['krt'] = unscramble(data)
        logger.debug("get_packet: unscrambled packet %d", pkt['raw pkt cnt'])
        ordered_queue.put(pkt)
      except (KeyboardInterrupt):
        # wait for reader to finish
        pass

if __name__ == "__main__":
  logging.basicConfig()
  mylogger = logging.getLogger()
  mylogger.setLevel(logging.DEBUG)

  socket = _open_socket()
  in_queue = Queue()
  ordered_queue = Queue()

  unpacker_queues = {} # one unpacker_queue for each worker
  worker = {}
  for count in range(10):
    unpacker_queues[count] = Queue()
    worker[count] = Process(target=unscramble_packet,
                            args=(in_queue, unpacker_queues[count]))
  reader = Process(target=get_packet, args=(socket, unpacker_queues))

  begin = time.time()
  mylogger.debug("Starting at %s", begin)
  for count in range(10):
    worker[count].start()
    mylogger.debug("started %d", count)
    time.sleep(0.1)
  reader.start() # get all the others going before starting the reader
  self.logger.debug("reader started")

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
  
  

