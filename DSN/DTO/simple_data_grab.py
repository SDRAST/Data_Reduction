"""
read packets from 10Gbe port, reformat, and store to disk


"""
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
  try:
    reader.terminate()
  except AttributeError:
    logger.error("signal_handler: no reader to terminate")
  for count in range(10):
    try:
      worker[count].terminate()
    except AttributeError:
      logger.error("signal_handler: no %s to terminate", worker[count])
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

def get_packet(socket, unpacker_queues):
  """
  gets packets and assigns them to unscramblers
  """
  while True:
    try:
      for worker in range(10):
        logger.debug("get_packet: getting data for worker %d", worker)
        for count in range(10000):
          data, addr = socket.recvfrom(pkt_size)
          unpacker_queues[worker].put(data)
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
      count = 10000
      one_second = []
      while count:
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
          one_second.append(pkt)
          count -= 1
        except (KeyboardInterrupt):
          # wait for reader to finish
          pass
      ordered_queue.put(one_second)
      logger.debug("get_packet: unscrambled 10000 packets")


if __name__ == "__main__":
  logging.basicConfig()
  mylogger = logging.getLogger()
  mylogger.setLevel(logging.INFO)

  socket = _open_socket()
  ordered_queue = Queue()

  unpacker_queues = {} # one unpacker_queue for each worker
  worker = {}
  for count in range(10):
    unpacker_queues[count] = Queue()
    worker[count] = Process(target=unscramble_packet,
                            args=(unpacker_queues[count], ordered_queue))
  reader = Process(target=get_packet, args=(socket, unpacker_queues))

  begin = time.time()
  mylogger.debug("Starting at %s", begin)
  for count in range(10):
    worker[count].start()
    mylogger.debug("started %d", count)
    time.sleep(0.1)
  reader.start() # get all the others going before starting the reader
  mylogger.debug("reader started")

  one_second = []
  not_done = True
  while not_done:
    try:
      one_second.append(ordered_queue.get())
      print "got one_second"
    except KeyboardInterrupt:
      not_done = False
  
  reader.join()
  for count in range(10):
    worker[count].join()
  processor.join()
  
  

