"""
Capture packets from 10Gbe port, unscramble, assemble into datasets and save

The 'reader' process captures 10000 packets and puts them round-robin into the
'unpacker' queues.  Each group of four unpackers has a common 'team' queue into
which it places the results.  Each 'team' just moves incoming packets to the
final 'ordered_queue'.
"""
import h5py
import logging
import numpy
import signal
import socket
import sys
import time

from multiprocessing import Lock, Process, Queue
from struct import unpack_from

from support import sync_second

logger = logging.getLogger(__name__)

pkt_size = 1026*64/8
pkt_fmt = "!8cHHI"+1024*4*"H"
if socket.gethostname() == 'gpu1':
  IP = '10.0.0.12'
else:
  IP = '10.0.0.2'
num_workers = 16
num_teams = 4

def signal_handler(signal, frame):
  """
  This does not end the thread
  """  
  print 'You pressed Ctrl+C!'
  socket.close()
  f.flush()
  try:
    reader.terminate()
  except AttributeError:
    logger.error("signal_handler: no reader to terminate")
  for count in range(num_workers):
    try:
      unpacker[count].terminate()
    except AttributeError:
      logger.error("signal_handler: no %s to terminate", unpacker[count])
  for count in range(num_teams):
    try:
      team[count].terminate()
    except AttributeError:
      logger.error("signal_handler: no %s to terminate", team[count])
  f.close()
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
      for unpacker in range(num_workers):
        logger.debug("get_packet: getting data for worker %d", unpacker)
        logger.debug("get_packet: putting data on %s at %f",
                     unpacker_queues[unpacker], time.time())
        for count in range(10000):
          data, addr = socket.recvfrom(pkt_size)
          unpacker_queues[unpacker].put(data)
        logger.debug("get_packet: finished 10000 packets at %f",
                     time.time())
    except KeyboardInterrupt:
      # nothing to do; signal_handler takes care of it
      pass
      
def unscramble_packet(input_queue, output_queue):
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
      one_second = {}
      one_second['hdr'] = []
      one_second['pkt cnt sec'] = []
      one_second['sec cnt'] = []
      one_second['raw pkt cnt'] = []
      one_second['pwr-I'] = []
      one_second['krt-I'] = []
      one_second['pwr-Q'] = []
      one_second['krt-Q'] = []
      count = 10000
      while count:
        try:
          pkt_buf = input_queue.get()
          if count == 10000:
            one_second['time'] = time.time()
          result = unpack_from(pkt_fmt, pkt_buf)
          one_second['hdr'].append(result[:8])
          one_second['pkt cnt sec'].append(result[8])
          one_second['sec cnt'].append(result[9])
          one_second['raw pkt cnt'].append(result[10])
          data = result[11:]
          power, kurtosis = unscramble(data)
          one_second['pwr-I'].append(power['I'])
          one_second['krt-I'].append(kurtosis['I'])
          one_second['pwr-Q'].append(power['Q'])
          one_second['krt-Q'].append(kurtosis['Q'])
          count -= 1
        except (KeyboardInterrupt):
          # wait for reader to finish
          pass
      output_queue.put(one_second)
      logger.debug("unscramble_packet: unscrambled 10000 packets from %s at %f",
                    input_queue, one_second['time'])
      logger.debug("unscramble_packet: unscrambling ended at %f", time.time())

def move_data(inqueue, outqueue):
  while True:
    try:
      data = inqueue.get()
      logger.debug("move_data: got data from %s at %s", inqueue, data['time'])
      logger.debug("move_data: sent data to %s at %s", outqueue, data['time'])
      outqueue.put(data)
    except KeyboardInterrupt:
      # nothing to do; signal_handler takes care of it
      pass
    
if __name__ == "__main__":
  logging.basicConfig()
  mylogger = logging.getLogger()
  mylogger.setLevel(logging.DEBUG)

  socket = _open_socket()
  ordered_queue = Queue()
  lock = Lock()
  
  unpacker_queue = {} # one unpacker_queue for each worker
  unpacker = {}
  team_queue = {}
  team = {}
  for count in range(num_workers):
    unpacker_queue[count] = Queue()
    if count % num_teams == 0:
      # define the team and team queue for this unpacker
      teamID = count/num_teams
      team_queue[teamID] = Queue()
      team[teamID] = Process(target=move_data,
                             args=(team_queue[teamID], ordered_queue))
      mylogger.debug("team %d takes data from %s", teamID, team_queue[teamID])
      mylogger.debug("team %d puts data on %s", teamID, ordered_queue)
    unpacker[count] = Process(target=unscramble_packet,
                              args=(unpacker_queue[count],
                                    team_queue[teamID]))
    mylogger.debug("unpacker %d takes data from %s",
                   count, unpacker_queue[count])
    mylogger.debug("unpacker %d puts data on %s", count, team_queue[teamID])
  reader = Process(target=get_packet, args=(socket, unpacker_queue))

  f = h5py.File('one_second.hdf5')
  
  for count in range(num_teams):
    team[count].start()
    mylogger.debug("started team %d", count)
  for count in range(num_workers):
    unpacker[count].start()
    mylogger.debug("started unpacker %d", count)
  sync_second()
  reader.start() # get all the others going before starting the reader

  not_done = True
  grpnum = 0
  while not_done:
    try:
      one_second = ordered_queue.get()
      grpname = "one_second %5d" % grpnum
      grp = f.create_group(grpname)
      for d in one_second.keys():
        ds = grp.create_dataset(d, data=numpy.array(one_second[d]))
      print("got %s" % grpnum)
      grpnum += 1
    except KeyboardInterrupt:
      not_done = False
  
  reader.join()
  for count in range(num_workers):
    unpacker[count].join()
  for count in range(num_teams):
    team[count].join()

  
  

