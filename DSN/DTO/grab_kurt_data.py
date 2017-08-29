"""
Monitors the 10 GBe interface and grabs packets.

The packets are 1026*64/8 = 8028 bytes long.
"""
import socket
import time

pkt_size = 1026*64/8

if __name__ == "__main__":
  port = 60000
  host = '10.0.0.12'

  s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
  s.bind((host,port))
  print "socket bound"
  
  start=time.time() + 1
  tmp = time.gmtime(start)
  tmp2 = time.strftime("%Y%m%d%H%M%S",tmp)

  i=1;
  running = True
  while (running):
    try: 
      filename = "/var/tmp/%s_%d" % (tmp2, i)
      outfile = open (filename, "wb")
      while (time.time()-start < 60*i and running):
        try:
          print "Waiting..."
          data, addr = s.recvfrom(pkt_size)
          print 'Packet len: ', len(data)
          if outfile:
            outfile.write (data)
        except (socket.timeout):
          print 'Socket timeout...'
          break
        except (KeyboardInterrupt):
          running = False
        i=i+1
      outfile.close()
    except (KeyboardInterrupt):
      running = False  
  
  if outfile:
    outfile.close()
  s.close()

