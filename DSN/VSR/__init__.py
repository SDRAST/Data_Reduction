"""
VSR Data Reduction

Background
==========

Each "side" (A and B) of a VSR is a virtual VSR. It produces a log called
vrt[0|1]_log.{DOY}{HHMM[SS]} which records the script configuration commands
issued and the responses to them.

Each virtual VSR side has two 16-MHz wide hardware channels. Each channel
can have up to four subchannels, two "narrow" (xN1 and xN2) and two "wide"
(xW1 and xW2) where x is the channel number.  In the vrt logs the subchannels
are numbered 0-7 for 1N1, 1N2, 1W1, 1W2, 2N1, 2N2, 2W1, 2W2.
There is a limit of four subchannels total for both sides.  A typical
configuration for astronomy with well separated frequencies might be::
  1W1, 1W2, 2W1, 2W2.

So a typical data hierarchy is::
  files:           vrt0       vrt1
  virtual VSRs:   VSRxA      VSRxB     where x = 1 only for DSS13
                 --------   --------
  vrt channels    2    6     2    6
  subchannels    1W1  1W2   2W1  2W2

We manage this with nested dictionaries::
  {"1A":{"1W1":{data},"1W2":{data}}, "1B":{"1W1":{data},"2W1":{data}}}

Raw VSR data files in UNIX format
=================================

A subchannel record structure consists of::

    typedef struct
    {
      schan_rcd_hdr_t hdr;         /* record header */
      U32 data[1];                 /* first entry in variably sized array */
    }
    schan_rcd_t;

This also has functions for handling UNIX format VSR raw data files.  These
each have data for one subchannel. Thesubchannel record header structure
consists of::

    typedef struct
    {
       U32 label;                       /* four byte label always = 'VSRD' */
       U32 record_length;               /* the length of the record in bytes */
       U16 record_version;              /* version number (increment with changes
                                           to this file) */
       U16 software_version;            /* version number (increment with changes
                                           to this file */
       U16 spc_id;                      /* station id - 10, 40, 60, 21 */
       U16 vsr_id;                      /* vsr1a, vsr1b ... from enum */
       U16 schan_id;                    /* subchannel id 0,1,2,3 */
       U16 bits_per_sample;             /* number of bits per sample - 1,2,4,8,
                                           or 16 */
       U16 sample_rate;                 /* number of samples per second in
                                       kilo-samples per second */
       U16 sdplr;                       /* algorithm used to apply doppler
                                       correction to sfro */
       U16 prdx_dss_id;                 /* DSS number */
       U16 prdx_sc_id;                  /* spacecraft number listed in predicts
                                           file */
       U16 prdx_pass_num;               /* pass number listed in predicts file */
       char prdx_uplink_band[2];        /* uplink RF band */
       char prdx_downlink_band[2];      /* downlink RF band */
       U16 trk_mode;                    /* tracking mode */
       U16 uplink_dss_id;              /* uplink dss id when in 3 way */
       U16 ddc_lo;                      /* DDC LO frequency in MHz */
       U16 rf_to_if_lo;                 /* RF to IF LO frequency in MHz */
       U16 data_error;                  /* hw error flag, dma error or
                                           num_samples error, 0 ==> no errors */
       U16 year;                        /* time tag - year */
       U16 doy;                         /* time tag - day of year */
       I32 sec;                         /* time tag - second of day */
       I32 data_time_offset;            /* in nano-seconds */
       double frov;                     /* in Hz */
       double fro;                      /* in Hz */
       double frr;                      /* in Hz/sec */
       double sfro;                     /* in Hz */
       double rf_freq;                  /* in Hz */
       double schan_accum_phase;        /* number of accumulated full turns */
       double schan_phase_poly[4];      /* coefficients for baseband phase
                                           polynomial */
       char schan_label[16];            /* ascii label - indicates source of data
                                           - carrier, subcarrier, etc */
    }
    schan_rcd_hdr_t;
    
This unpacks with ">4s l 11h 2s 2s 7h 2l 10d 16s"
"""

import struct
import DatesTimes as DT
import time
from glob import glob
from Data_Reduction.DSN import obs_dir, ravi_data_dir

subchannel = ["1N1", "1N2", "1W1", "1W2", "2N1", "2N2", "2W1", "2W2"]

BITS_PER_32BIT_WORD = 32
VSR_HEADER_SIZE = 152
diag = False

def find_VSR_log_block_times(ses_date,year,month,day,DOY):
  """
  Get the recording start and stop times
  
  This examines the VSR log files to get the recording start and
  stop times.  The log files have the start and stop times for each
  subchannel of the VSR.  The output is written to files.

  @param ses_date : string -
    Date of observations as 'YYYY-MM-DD'

  @param year : int -
    Year of observation

  @param month : int -
    Month of observation

  @param day : int -
    Day of the month

  @param DOY : int -
    Day of year

  @return: True|False -
    True if succesful.  Otherwise False
  """
  TT = VSR_to_timetuple((year,DOY,0))
  datestr = time.strftime("%y-%j",TT)
  pattern = obs_dir+ses_date+"/vrt*log."+("%03d" % DOY)+"*"
  datafiles = glob(pattern)
  if datafiles != []:
    for datafile in datafiles:
      if diag:
        print "find_VSR_log_block_times: processing",basename(datafile)
      # First we need to match sub-channel numbers with IDs we understand
      channels = {}
      response = os.popen("grep 'fft_file' "+datafile,'r').readlines()
      for line in response:
        parts = line.strip().split()
        chan = int(parts[4].strip(':'))
        chanID = parts[-1][-6:]
        channels[chan] = chanID
      if diag:
        print channels
      rec_times = os.popen("grep 'rec set to' "+datafile,'r').readlines()
      start_times = {}
      stop_times = {}
      subchans = []
      for rec_time in rec_times:
        timetuple = strptime(rec_time[:16],"%y/%j %H:%M:%S>")
        TS = time.mktime(timetuple)
        parts = rec_time[16:].strip().split()
        chan  = int(parts[1])
        state = int(parts[5])
        try:
          pos = subchans.index(chan)
        except ValueError:
          pos = -1
          subchans.append(chan)
          start_times[chan] = []
          stop_times[chan] = []
          if diag:
            print "find_VSR_log_block_times: created lists for channel",chan
        if state == 0:
          stop_times[chan].append(TS)
        elif state == 1:
          start_times[chan].append(TS)
        else:
          print "find_VSR_log_block_times: unrecognized state",state
          sys.exit(0)
        if diag:
          print timetuple,chan,state,"index",pos
      if diag:
        print "find_VSR_log_block_times - start times:",start_times
        print "find_VSR_log_block_times - stop times:",stop_times
      for chan in start_times.keys():
        outfile = obs_dir+ses_date+"/recording_times."+channels[chan]
        if diag: print "find_VSR_log_block_times: writing output to",\
          basename(outfile)
        outfd = open(outfile,'w')
        outfd.write("#Processing " + basename(datafile)
                                  + " subchannel " + str(chan)
                                  + " for block times\n")
        if diag:
          print "find_VSR_log_block_times: recording",basename(datafile),\
            "subchan",chan
        for index in range(len(start_times[chan])):
          if diag:
            print "From",ctime(start_times[chan][index]),"to",
            print        ctime(stop_times[chan][index])
          outfd.write(str(start_times[chan][index])
                      + ' '
                      + str(stop_times[chan][index]))
          outfd.write(" # From "+ctime(start_times[chan][index])+" to ")
          outfd.write(ctime(stop_times[chan][index])+'\n')
        outfd.close()
    return True
  else:
    return False

def get_RAVI_spectra(year,DOY,dest_path):
  """
  Copies spectrum files from the RAVI to the local host.

  Spectra produced with the program 'fft' on the RAVI are
  named 'FFT_vsr*', at least for the PESD project.

  @param year : int -
    Year of observation.

  @param DOY : int -
    Day of observation.

  @param dest_path : string -
    Where the files will go on the local host.
  """
  spectra = glob(ravi_data_dir
                 +"FFT_vsr*"+str(year)[2:]+("-%03d" % DOY)+"*")
  for spectrum in spectra:
    shutil.copy(spectrum,dest_path)
    if diag:
      print basename(spectrum),"copied"
  return spectra

def get_raw_record_header(fd):
  """Reads the header of a UNIX format raw VSR data file.

  Notes
  =====
  The file pointer
  mist be positioned at the front of the header.

  @param fd : file descriptoy
  """
  buf = fd.read(VSR_HEADER_SIZE)
  label,reclen,recver,softver,spcid,vsrid,schanid,bps,srate,sdplr, \
    prdx_dss_id,prdx_sc_id,prdx_pass_num,prdx_uplink_band, \
    prdx_downlink_band,trk_mode,uplink_dss_id,ddc_lo,rf_to_if_lo, \
    data_error,year,doy,sec,data_time_offset,frov,fro,frr,sfro, \
    rf_freq,schan_accum_phase,scpp0,scpp1,scpp2,scpp3,schan_label \
    = struct.unpack_from(">4sl11h2s2s7h2l10d16s",buf)
  return label,reclen,recver,softver,spcid,vsrid,schanid,bps,srate,sdplr, \
    prdx_dss_id,prdx_sc_id,prdx_pass_num,prdx_uplink_band, \
    prdx_downlink_band,trk_mode,uplink_dss_id,ddc_lo,rf_to_if_lo, \
    data_error,year,doy,sec,data_time_offset,frov,fro,frr,sfro, \
    rf_freq,schan_accum_phase,(scpp0,scpp1,scpp2,scpp3),schan_label

def get_raw_record(fd, filesize):
  """
  Obtains a record from the raw data file specified.
  
  This obtains one record from the current file position. A record consists
  of the header
  and one second of raw complex data samples.  The number of samples is
  'ksamps_per_sec*1000.  This leaves the file position pointer at the beginning
  of the next record.
    
  Notes
  =====
  To get a specific record from an arbitrary position use::
    get_binary_record(fd,file_header_size,record_size,index)

  @param fd : file descriptor

  @param filesize : long -
    To prevent attempts to read past the end of the file.
  """
  # need to make sure that there is enough data left in file
  pos = fd.tell()
  if pos + VSR_HEADER_SIZE > filesize:
    return [],[]
  header = get_raw_record_header(fd)
  bits_per_sample = header[7]
  ksamps_per_sec = header[8]
  reclen = header[1]
  if diag:
    print "get_raw_record:",ksamps_per_sec,"complex ksamples/sec"
  cmplx_samp_per_word = (BITS_PER_32BIT_WORD/(2*bits_per_sample))
  if diag:
    print cmplx_samp_per_word,"complex",bits_per_sample,\
    "bit samples per",str(BITS_PER_32BIT_WORD)+"-bit word"
  num_words = ksamps_per_sec*1000/cmplx_samp_per_word
  if diag:
    print "Reading",num_words,str(BITS_PER_32BIT_WORD)+"-bit words"
  try:
    data = fd.read(num_words*4) # bytes
  except MemoryError, details:
    print "MemoryError:",details
    print "get_raw_record: failed to get record a file position",fd.tell()
    print "get_raw_record:",ksamps_per_sec,"complex ksamples/sec"
    print cmplx_samp_per_word,"complex",bits_per_sample,\
    "bit samples per",str(BITS_PER_32BIT_WORD)+"-bit word"
    print "Reading",num_words,str(BITS_PER_32BIT_WORD)+"-bit words"
    return [],[]
  if diag:
    print len(data),"get_raw_record: data bytes read"
  return header, data

def time_to_index(rec_time,blk_times,file_start,file_end):
  """
  Convert a UNIX timestamp to a record index.
  
  Given a record time as a UNIX timestamp, it returns the index
  of the record in the file.  It's a bit complicated because the
  recording block times are obtained from the log file but the
  actual recorded data may not conform to that due to disk drive
  errors.

  @param rec_time : float -
    UNIX timestamp

  @param blk_times : list -
    Tuples of recording start/stop time pairs.

  @param file_start : float -
    UNIX timestamp for first record in file

  @param file_end : float -
    UNIX timestamp for last record in file
  """
  if rec_time < file_start:
    print "Before start of recording"
    return -1
  elif rec_time > file_end:
    print "After end of recording"
    return -1
  else:
    # initialize the index to the start of the file
    index = 0
    for pair in blk_times:
      if rec_time < pair[0]:
        # Time was during a period with no recording
        return -1
      elif rec_time < pair[1]:
        # Time is in this recording interval
        index += (rec_time - pair[0])
        break
      else:
        # Time is after the current interval; move to the next one
        index += (pair[1] - pair[0])
  return index

def recording_times(datadir,vsr,channel):
  """
  Get recording times from a recording times file.

  This reads a recording times and extracts the blocks of time for
  which means data exist.  It returns a list of [start,stop] pairs.
    
  Notes
  =====
  A recording times file consists of a series of blocks like this::
    Processing vsr1b.2w1.10-113-2200-2210.log for block times
    From Fri Apr 23 22:10:00 2010 to Fri Apr 23 23:59:59 2010
    From Sat Apr 24 00:15:00 2010 to Sat Apr 24 01:59:59 2010
    From Sat Apr 24 02:15:00 2010 to Sat Apr 24 04:59:59 2010
    From Sat Apr 24 05:15:00 2010 to Sat Apr 24 07:14:21 2010
    From Sat Apr 24 07:14:23 2010 to Sat Apr 24 07:14:59 2010

  @param datadir : string -
    Directory for an obsering day

  @param vsr: string -
    A VSR identifier, like "1A".

  @param channel : string -
    A channel identifier, like "2"
  """
  fd = open(datadir+"/recording_times."+vsr.upper()+"."+channel.upper(),"r")
  finished = False
  times = []
  total_time = 0.
  while not finished:
    stripped_line = fd.readline().strip()
    if diag:
      print "Line:",stripped_line
      print "length=",len(stripped_line)
    if len(stripped_line) > 0:
      if stripped_line[0] != "#":
        items = stripped_line.split()
        start = float(items[0])
        end   = float(items[1])
        times.append([start,end])
        block_time = end - start
        if diag:
          print block_time,"samples in this block"
        total_time += block_time
    else:
      finished = True
  if diag:
    print total_time,"samples total"
  return times

subchannel = ["1N1", "1N2", "1W1", "1W2", "2N1", "2N2", "2W1", "2W2"]

def parse_VSR_log(ses_date,year,month,day,DOY):
  """
  Parse VSR 'vrt*' log files for observing parameters.
  
  This examines the VSR log files to get observing paramaters such as
  frequency and sampling rate.  It also gets the recording start and
  stop times.  The log files have the start and stop times for each
  subchannel of the VSR.

  Notes
  =====
  This supercedes find_VSR_log_block_times(), which only writes the
  start and stop times to a disk file.

  @param ses_date : string -
    Date of observations as 'YYYY-MM-DD'

  @param year : int -
    Year of observation

  @param month : int -
    Month of observation

  @param day : int -
    Day of the month

  @param DOY : int -
    Day of year

  @return: dictionary -
    The dictionary is keyed to the VSR id.  The values of the VSR
    dictionaries are also dictionaries, keyed to the VSR channel number.
    The values of these VSR dictionaries are also dictionaries with
    keys like:
      - 'dss': station identifier
      - 'vsr': VSR identifier
      - 'frov':
      - 'LO':
      - 'sfro':
      - 'restfreq':
      - 'start':
      - 'stop':
  """
  TT = DT.VSR_to_timetuple((year,DOY,0))
  datestr = time.strftime("%y-%j",TT)
  pattern = obs_dir+ses_date+"/vrt*log."+("%03d" % DOY)+"*"
  if diag:
    print "Looking for",pattern
  datafiles = glob(pattern)
  if diag:
    print "Found:",datafiles
  vsr_data = {}
  if datafiles != []:
    for datafile in datafiles:
      if diag:
        print "\nparse_VSR_log: processing",basename(datafile)
      df = open(datafile,'r')
      lines = df.readlines()
      df.close()
      # The line with bandwidth does not have a subchannel identifier.
      # The previous line does.
      expecting_bandwidth = 0
      for line in lines:
        parts = line.strip().split()
        date_str = parts[0]
        time_str = parts[1].rstrip('>')
        # The following apply to all channels
        if len(parts) > 7 and parts[2] == "spc_id":
          # This is encountered very early in the file
          dss = parts[4].rstrip(',')
          vsr = parts[7]
          # create a dictionary for this VSR
          vsr_data[vsr] = {}
          if diag:
            print "DSS",dss,"VSR",vsr
        elif len(parts) > 7 and parts[2] == "year":
          # The year is needed since it's not in the time stamp
          year = 2000+int(parts[4].rstrip(','))
          DOY = int(parts[7])
          if diag:
            print year,DOY
        elif len(parts) > 4 and parts[2] == "rf_to_if_lo":
          # This is assumed to apply to all subchannels
          LO = int(parts[4])
        elif len(parts) > 5 and parts[2] == 'frov':
          # This applies to all subchannels
          frov = parts[-1]
        elif parts[5:8] == ['configuring', 'ASD', 'slot']:
          # the rest of the 'elif's below apply on a
          # channel-by-channel basis
          chan_num = int(parts[3])
          chan = subchannel[chan_num]
          expecting_bandwidth = chan
          if not vsr_data[vsr].has_key(chan):
            vsr_data[vsr][chan] = {}
            if diag:
              print "Channel",chan,"has been added"
            vsr_data[vsr][chan]["start"] = []
            vsr_data[vsr][chan]["stop"] = []
        elif len(parts) > 8 and parts[2] == "bandwidth" \
              and expecting_bandwidth != 0:
          vsr_data[vsr][expecting_bandwidth]['bw'] = int(parts[4])
          vsr_data[vsr][expecting_bandwidth]['bpsamp'] = \
              int(parts[8].strip(','))
          expecting_bandwidth = 0
        elif len(parts) > 7 and parts[4] == 'sfro':
          chan = subchannel[int(parts[3])]
          vsr_data[vsr][chan]['sfro'] = float(parts[-1])
        elif len(parts) > 7 and parts[4] == 'label':
          chan = subchannel[int(parts[3])]
          vsr_data[vsr][chan]['restfreq'] = float(parts[-1])
        elif len(parts) > 6 and parts[6] == 'recfn':
          chan = subchannel[int(parts[7])]
          vsr_data[vsr][chan]["datafile"] = parts[8].strip('"')
        elif parts[4:7] == ['rec','set','to']:
          chan = subchannel[int(parts[3])]
          state = int(parts[7])
          timetuple = time.strptime(parts[0]+parts[1],"%y/%j%H:%M:%S>")
          TS = time.mktime(timetuple)
          if state == 1:
            vsr_data[vsr][chan]["start"].append(TS)
            vsr_data[vsr][chan]["LO"] = LO
            vsr_data[vsr][chan]["frov"] = frov
            vsr_data[vsr][chan]["dss"] = dss
            vsr_data[vsr][chan]["vsr"] = vsr
            if diag:
              print "Chan",chan,"recording started at",time.ctime(TS)
              print "LO =", vsr_data[vsr][chan]["LO"]
          elif state == 0:
            vsr_data[vsr][chan]["stop"].append(TS)
            if diag:
              print "Chan",chan,"recording stopped at",time.ctime(TS)
              print "Channel",chan,"recorded to",basename(datafile)
          else:
              print "Invalid recording state",state
  else:
    vsr_data = {}
  return vsr_data
