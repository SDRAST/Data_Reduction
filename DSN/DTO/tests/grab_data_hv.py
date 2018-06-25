"""
Monitors the 10 GBe interface(s) and grabs packets.
"""
import socket
#from time import sleep, time
#import numpy as np
#from optparse import OptionParser
import time
import telnetlib

roach1 = '192.168.100.41'
roach2 = '192.168.100.42'
katcp_port = 7147

if __name__ == "__main__":
  roach = roach1
  tn = telnetlib.Telnet(roach1,katcp_port)
  buf =  9000 # maximum packet size
  port = 4001
  port1 = 4000
  thishost = socket.gethostname()
  if thishost == 'gpu1':
    host = '10.0.0.12' # roach1
    host1 = '10.0.0.13' #'193.166.42.23'
  elif thishost == 'gpu2':
    host = '10.0.0.2' # roach2
    host1 = '10.0.0.3'

  s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
  s1 = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
  s.bind((host,port))
  print "socket bound"
  s1.bind((host1,port1))
  print "socket 1 bound"
  start=time.time() + 1
  tmp = time.gmtime(start)
  tmp2 = time.strftime("%Y%m%d%H%M%S",tmp)
  tn.write("?wordwrite adc_reg_arm 0 1\n")
  tn.write("?wordwrite adc_reg_arm 0 0\n")
  tn.write("?wordwrite adc_reg_arm 0 1\n")

  i=1;
  starttime = time.time()
  while (1):
    try: 
      filename = "/var/tmp/%s_%d" % (tmp2, i)
      filename1 = "/var/tmp/%s_%d" % (tmp2, i)
      outfile = open (filename, "wb")
      outfile1 = open (filename1, "wb")
      while (time.time()-start < 60*i):
        try:
          print "Waiting..."
          data, addr = s.recvfrom(8192)
          data1, addr1 = s1.recvfrom(8192)
          print "Got it!"
          #seq = int(data[1],2)
          #print "seq: %u" % seq 
          #print "%u" int(data[0],2)
          #fdata=np.append(fdata,data)
   
          if outfile: outfile.write (data)
          if outfile1: outfile1.write (data1)
        except (socket.timeout):
          print 'Socket timeout...'
          break
        except (KeyboardInterrupt):
          print "%d packets in %d seconds" % (i,time.time()-start)
          exit()
      i=i+1
      outfile.close()
    except (KeyboardInterrupt):
      print "%d packets in %d seconds" % (i,time.time()-start)
      exit()  
   
    if outfile: outfile.close()
    if outfile1: outfile1.close()

    print 'Packet len: ', len(data)