"""
Support for STATS data files.

STATS binary file structure
===========================

A stats binary output files begins with a stats_hdt_t structure::
  
  typedef struct
  {
    (unsigned short header_size   /* bytes, may or may not be there */
    unsigned short spcid;         /* station id - 10, 40, 60, 21 */
    unsigned short vsrid;         /* vsr1a, vsr1b ... from enum */      
    unsigned short chanid;        /* subchannel id 0,1,2,3 */
    unsigned short bps;           /* number of bits per sample - 1, 2, 4, 8,
                                     or 16 */
    unsigned long  srate;         /* number of samples per second in kilo-
                                     samples per second */
    unsigned short error;         /* hw err flag, dma error or num_samples
                                     error,
                                     0 ==> no errors */
    unsigned short year;          /* time tag - year */
    unsigned short doy;           /* time tag - day of year */
    unsigned long  sec;           /* time tag - second of day */
    double         freq;          /* in Hz */
    unsigned long  orate;         /* number of statistics samples per
                                     second */
    unsigned short nsubchan;      /* number of output sub chans */
  }
  stats_hdr_t;

This unpacks with "=4H Q HHH Q d Q H"

A data record looks like this::
        fwrite(&doy,   sizeof(int),    1, ofp);
        fwrite(&sec,   sizeof(int),    1, ofp);
        fwrite(&i,     sizeof(int),    1, ofp);
        fwrite(&mean,  sizeof(double), 1, ofp);
        fwrite(&var,   sizeof(double), 1, ofp);
        fwrite(&skew,  sizeof(double), 1, ofp);
        fwrite(&kurt,  sizeof(double), 1, ofp);
        fwrite(&mean,  sizeof(double), 1, ofp);
        fwrite(&var,   sizeof(double), 1, ofp);
        fwrite(&skew,  sizeof(double), 1, ofp);
        fwrite(&kurt,  sizeof(double), 1, ofp);

which unpacks with "=LLL dddd dddd"


STATS ASCII file structure
==========================

A STATS file begins with a header.  The data lines consist of
  - column 0: second, starting with 1
  - column 1: sample number within the second (typically 0-999)
  - column (subchannel + 2): mean
  - column (subchannel + 3): r.m.s.
  - column (subchannel + 4): kurtosis
  - column (subchannel + 5): skewness
where subchannel = 0, 1.
"""
from DatesTimes import *
from time import *
from glob import glob
from Data_Reduction.DSN import obs_dir, \
                               get_binary_record, \
                               STATS_binary_record_size
from os.path import basename
from numpy import array, append

diag = True
diag_read = False

def process_STATS_ASCII_header(fd):
  """Process the header in a STATS file
   
  STATS files are created from VSR data.  The earliest
  versions of a STATS file did not preface header data
  with #. This was added later for use with some plotting
  programs.

  @param fd : file descriptor

  @return: dictionary
    The keys are the same ones as in the header.
  """
  header = {}
  doing_header = True
  while(doing_header):
    line = fd.readline().strip()
    if re.search('::',line):
      [k,v] = line.split('::')
      key = k.strip().lstrip('#')
      if re.search('\.',v) == None:
        header[key] = int(v.strip())
      else:
        header[key] = float(v.strip())
    elif re.search('HEADER_END',line):
      doing_header = False
  return header

def process_STATS_ASCII_data_line(line,nchan):
  """This processes one line of a STATS data file

  @param line : string
    One line from an ASCII STATS data file

  @param nchan : number of signal channels in the file
  
  @return: list of lists
    Lists of the means, rms, kurtoses and skews for the subchannels
    at one sampling time.
  """
  data = line.split()
  mean = []
  rms = []
  skew = []
  kurtosis = []
  sec = int(data[0])
  ms = int(data[1])
  for i in range(2,2+4*nchan,4):
    mean.append(float(data[i]))
    rms.append(float(data[i+1]))
    kurtosis.append(float(data[i+2]))
    skew.append(float(data[i+3]))
  return (sec,ms,mean,rms,kurtosis,skew)

def get_STATS_ASCII_data_block(fd,nchan,nsamps):
  """This reads data for one block from the data file.
    
  This should sit and wait until data are available, line by line.
  When the required number of lines have been read it
  returns the data as a tuple of arrays.

  @param fd : file descriptor

  @param nchan : int
    Number of data channels processed

  @param nsamps : int
    Number of samples in a block

  @return: tuple
    The tuple consists of five arrays
    (means,rootmeansquare,kurtosis,skewness,sec)
    Each array shape is (nsamps,nchan).
  """
  if diag_read:
    print("Reading",fd)
  counter = 0
  while(counter < nsamps):
    fd.flush()
    line = fd.readline()
    if line == '':
      # end-of-file
      return zeros(1),zeros(1),zeros(1),zeros(1),0
    # Handle incomplete lines
    while line[-1] != '\n':
      fd.flush()
      line += fd.readline()
    # process the line
    line = line.strip()
    if line != '':
      if diag_read:
        print("Read:",line)
      # mean, rms, kurt and skew are lists whose length is the number of
      # channels
      sec,ms,mean,rms,kurt,skew = process_STATS_ASCII_data_line(line,nchan)
      if counter == 0:
        # initialize the arrays
        # ndmin forces the arrays to have 2 dimensions, the first
        # dimension being 1 and the second num_subch.
        means = array(mean,ndmin=2)
        rootmeansquare = array(rms,ndmin=2)
        kurtosis = array(kurt,ndmin=2)
        skewness = array(skew,ndmin=2)
      else:
        # append to the moment arrays
        means = append(means,array(mean,ndmin=2),axis=0)
        rootmeansquare = append(rootmeansquare,array(rms,ndmin=2),axis=0)
        kurtosis = append(kurtosis,array(kurt,ndmin=2),axis=0)
        skewness = append(skewness,array(skew,ndmin=2),axis=0)
      counter += 1
  return means,rootmeansquare,kurtosis,skewness,sec

 
def get_data_block(fd,nchan,nsamps):
  """
  Alias for get_STATS_ASCII_data_block
   
  For backward compatibility

  @param fd : file descriptor

  @param nchan : int
    Number of data channels processed

  @param nsamps : int
    Number of samples in a block

  @return: tuple
    The tuple consists of five arrays
    (means,rootmeansquare,kurtosis,skewness,sec)
    Each array shape is (nsamps,nchan).
  """
  return get_STATS_ASCII_data_block(fd,nchan,nsamps)

def parse_STATS_ASCII_header(header):
  """
  Parses the header of a STATS
   
  @param header : dictionary
    Header dictionary of a STATS file.

  @return: tuple
    (year,doy,start_sec,freq,spc,vsr,nchan,bw,bps,nsamps)
  """
  year      = header['YEAR']
  doy       = header['DOY']
  start_sec = header['START_SEC']
  freq      = header['RF_FREQ[HZ]']/1.e6     # MHz
  spc       = header['SPC_ID']
  vsr       = header['VSR_ID']
  nchan     = header['NOCHAN']               # number of sub-channels
  bw        = header['SAMPLE_RATE[HZ]']/2.e6 # MHz
  bps       = header['BITS_PER_SAMPLE']
  nsamps    = header['OUTPUT_RATE[HZ]']      # output samples/sec
  return year,doy,start_sec,freq,spc,vsr,nchan,bw,bps,nsamps

def parse_STATS_header(header):
  """
  Extract the header from a binary STATS data file.

  This extracts the STATS binary file header information
  into variables with meaningful names.  It also converts year
  and day_of_year to seconds at midnight since the epoch used by
  UNIX systems.

  @param header : string of binary data

  @return: tuple
    (spcid, vsrid, chanid, bps, srate, errflg, year, doy, sec, freq,
    orate,nsubchan)
  """
  (spcid,             #  1) station id - 10, 40, 60, 21
   vsrid,             #  2) vsr1a, vsr1b ...
   chanid,            #  3) subchannel id 0,1,2,3
   bps,               #  4) number of bits per sample - 1, 2, 4, 8, or 16
   srate,             #  5) number of samples per second in samples per second
   errflg,            #  6) hardware error flag, dma error or num_samples
                      #     error, 0 ==> no errors
   year,              #  7) time tag - year
   doy,               #  8) time tag - day of year

   sec,               #  9) time tag - second of day

   freq,              # 10) frequency in Hz

   orate,             # 11) number of statistics samples per second

   nsubchan           # 12) number of output sub chans
   ) = header
  return spcid, vsrid, chanid, bps, srate, errflg, year, doy, sec, freq, \
         orate,nsubchan

def get_binary_stats_header(fd):
  """Get the header from a binary stats file.

  There is an old format in which the first datum is a short with the
  station ID of 13.  Otherwise the first long has a header size.  This
  handles either case.

  Notes
  =====
  Function parse_STATS_header() is one-liner that
  translates the tuple members to variables with meaningful names.

  @param fd : file descriptor.

  @return: tuple.
    A string with binary data followed by the size of the header (int).
  """
  first_word = fd.read(2)
  first = struct.unpack_from('=H',first_word)
  if diag_read:
    print("Header size =",first)
  # Unpack returns a tuple
  if first[0] != 13:
    # header_size = first[0]
    header_size = 52 # header length prepends header
    fd.seek(2,1) # advance past the long which alledgedly has the
                 # header size
    buf = fd.read(header_size-4)
  else:
    header_size = 48
    # read the remaining header
    buf = first_word + fd.read(header_size-2)
  # This will change if header_size changes
  header = struct.unpack_from('=4H Q HHH Q d Q H',buf)
  return header,header_size

def write_binary_stats_header(header):
  """
  Write a header in binary format.
  
  This packs a header into a buffer for creating a binary file.

  @param header : tuple

  @return: binary string
  """
  buf = struct.pack('=4H Q HHH Q d Q H',*header)
  return buf
  
def get_binary_stats_record(fd,header_size,index):
  """
  Extracts a binary record at the specified record index.

  If a particular time is wanted, then it is necessary to read the seconds
  since midnight and the index (usually milliseconds) to verify and
  possibly adjust the position.

  Notes
  =====
  Two data channels are assumed.

  @param fd : file descriptor.
  
  @param header_size : int.
    Header size in bytes.

  @param index : long.
    Index of record to be retrieved.

  @return: tuple.
    (DOY, start_sec, record_index , mean0, variance0, skewness0, kurtosis0,
    mean1, variance1, skewness1, kurtosis1)
  """
  buf = get_binary_record(fd,
                          header_size,
                          STATS_binary_record_size,
                          index)
  data = struct.unpack_from("=LLL dddd dddd", buf)
  return data

def get_STATS_block(fd,year,orate,blk_index):
  """
  Gets the signal statistics data for a one second block.

  Read Returns
  the sample times in UNIX time and arrays with the statistics for
  each channel.  The array shape is (n_samples,n_channels).

  Notes
  =====
  This can't work because 'header_size' is needed by
  get_binary_stats_record() but is not defined either locally
  or globally.
  """
  # position to the first record
  first_rec_index = blk_index*orate
  times = []
  means = []
  variances = []
  skews = []
  kurts = []
  # get all the records in this block.  There are 'orate' records.
  for i in range(first_rec_index,first_rec_index+orate):
    TS,mean,variance,skewness,kurtosis \
      = parse_record(year,orate,get_binary_stats_record(fd,header_size,i))
    times.append(TS)
    means.append(mean)
    variances.append(variance)
    skews.append(skewness)
    kurts.append(kurtosis)
  return array(times),array(means),array(variances),array(skews),array(kurts)
 
def parse_record(year,orate,data):
  """
  This parses a data record that is in list format.

  It converts the time data into a UNIX-like
  timestamp such as used in module 'time'.  Unlike the the UNIX timestamp
  it has resolution up to a microsec.

  @param year : int
    Year of observation

  @param orate : int
    Number of averages per second

  @param data : list

  @return:
    (UNIX time stamp, (mean0, mean1), (variance0,  variance1),
    (skewness0,  skewness1), (kurtosis0, kurtosis1))
  """
  # the doy may have advanced through midnight
  doy = data[0]
  TS0 = VSR_tuple_to_timestamp(year,doy,0)
  sec = data[1]
  rec_index = data[2]
  TS = TS0 + sec + float(rec_index)/orate
  ave1 = data[3]; ave2 = data[7]
  var1 = data[4]; var2 = data[8]
  skw1 = data[5]; skw2 = data[9]
  krt1 = data[6]; krt2 = data[10]
  return TS,(ave1,ave2),(var1,var2),(skw1,skw2),(krt1,krt2)

def find_STATS_bin_block_times(ses_date,year,month,day,DOY):
  """
  Gets the recording start and stop times from a binary STATS file.
  
  This looks for recording gaps in the binary data files, that is, for
  discontinuities in record 'start_sec'.  It writes the start and
  stop time pairs to a file called 'recording_times.STATS-bin' in the
  current data directory.  This one is very slow.  It's faster to examine
  the STATS log files - see 'find_STATS_log_block_times'.

  @param ses_date : string
    Date of observations as 'YYYY-MM-DD'

  @param year : int
    Year of observation

  @param month : int
    Month of observation

  @param day : int
    Day of the month

  @param DOY : int
    Day of year

  @return: None
  """
  TT = VSR_to_timetuple((year,DOY,0))
  datestr = strftime("%y-%j",TT)
  datafiles = glob(ravi_data_dir+"STATS*"+datestr+"*-bin")
  if diag:
    print("STAT bin files:\n",datafiles)
  for datafile in datafiles:
    fd = open(datafile,'r')
    st_mode, st_ino, st_dev, st_nlink, st_uid, st_gid, st_size, st_atime, \
	   st_mtime, st_ctime = os.stat(datafile)
    if diag:
      print("\nProcessing ",basename(datafile)," for block times")
      print("File size =",st_size)
    chanID = basename(datafile[15:22]).upper()
    outfile = obs_dir+ses_date+"/recording_times."+chanID
    if diag:
      print("Writing output to",basename(outfile))
      print("Be patient.  Have lunch.  Play solitaire.")
    outfd = open(outfile,'w')
    outfd.write("Processing "+basename(datafile)+" for block times\n")
    header,header_size = get_binary_stats_header(fd)
    num_recs = (st_size-header_size)/STATS_binary_record_size
    if diag:
      print("File data size =",st_size-header_size)
      print(STATS_binary_record_size,"bytes per record")
      print(num_recs,"samples of data")
    spcid, vsrid, chanid, bps, srate, errflg, year, doy, sec, freq, \
         orate,nsubchan = parse_STATS_header(header)
    if diag:
      print("Number of records per second =",orate)
      print((st_size-header_size)/STATS_binary_record_size/orate/60., \
            "minutes of data")
    # First recording start time
    TSprev = VSR_tuple_to_timestamp(year,doy,sec)
    if diag:
      print("Header time:",ctime(TSprev))
    sleep(5)
    outfd.write("From "+ctime(TSprev))
    for rec_id in range(1,num_recs,orate):
      data = get_binary_stats_record(fd,
                                     header_size,
                                     rec_id)
      doy = data[0]
      sec = data[1]
      TS = VSR_tuple_to_timestamp(year,doy,sec)
      if diag:
        print("Examining  record", rec_id,"\r", end=' ')
      if TS - TSprev > 1:
        if diag:
          print("Break at    record",rec_id, \
                "time =",doy,sec, \
                "\nfrom",ctime(TSprev),'to',ctime(TS))
        outfd.write(" to "+ctime(TSprev)+'\n')
        outfd.write("From "+ctime(TS))
      TSprev = TS
      # final end
    if diag:
      print("Finished at record",rec_id, \
            "time =",doy,sec, \
            'at',ctime(TS))
    outfd.write(" to "+ctime(TS)+'\n')
    fd.close()
    outfd.close()

def find_STATS_log_block_times(ses_date,year,month,day,DOY):
  """
  Gets the recording times from the STATS log files.
   
  This looks for gaps in the times that STATS data were recorded.
  This gives the blocks
  of data that have been processed by STATS.  It may differ from actual
  recording times if the program 'stats' started at a time later than the
  start of recording or terminated early.  The results are written to a
  file 'recording_times.STATS-log' in the current data directory.

  @param ses_date : string
    Date of observations as 'YYYY-MM-DD'

  @param year : int
    Year of observation

  @param month : int
    Month of observation

  @param day : int
    Day of the month

  @param DOY : int
    Day of year

  @return: None  
  """
  TT = VSR_to_timetuple((year,DOY,0))
  datestr = strftime("%y-%j",TT)
  datafiles = glob(obs_dir+ses_date+"/vsr1*"+datestr+"*log")
  if datafiles != []:
    if diag:
      print("find_STATS_log_block_times: STAT log files:\n",datafiles)
    for datafile in datafiles:
      if diag:
        print("\nfind_STATS_log_block_times: processing ",\
          basename(datafile)," for block times")
      chanID = basename(datafile)[4:9]
      outfile = obs_dir+ses_date+"/recording_times.1"+chanID.upper()
      if diag:
        print("find_STATS_log_block_times: writing output to",basename(outfile))
      outfd = open(outfile,'w')
      outfd.write("#find_STATS_log_block_times: processing "+\
        basename(datafile)+" for block times\n")
      fd = open(datafile,'r')
      not_EOF = True
      first_record = True
      while not_EOF:
        line = fd.readline()
        if line[:12] == "year:doy:sec" and len(line) < 30:
          # There could be a bad last line like:
          # year:doy:sec 2010:7/media/disk-5/PESD_data/STATS_NP1000_vsr1a.2w1.10-079-mars complete: Illegal seek
          datecode = line[12:]
          date_parts = datecode.split(":")
          year = int(date_parts[0])
          DOY = int(date_parts[1])
          sec_str = date_parts[2].split('#') # strips off any comment
          sec = int(sec_str[0])
          TS = VSR_tuple_to_timestamp(year,DOY,sec)
          if first_record:
            start = TS
            TSprev = TS
            first_record = False
            if diag:
              print("From ",ctime(start), end=' ')
          if TS - TSprev > 1:
            stop = TSprev
            outfd.write(str(start)+" "+str(stop)
                        +" # From "+ctime(start)+" to "+ctime(stop)+'\n')
            start = TS
            if diag:
              print("Break at",year,DOY,sec)
              print("to",ctime(TSprev))
              print() 
          TSprev = TS
        elif line == "":
          not_EOF = False
      outfd.write(str(start)+" "+str(stop)
                        +" # From "+ctime(start)+" to "+ctime(stop)+"\n")
      if diag:
        print("to",ctime(TSprev))
      fd.close()
      outfd.close()
    return True
  else:
    return False

def print_record(data):
  """
  Pretty-prints a STATS record that is in list format.

  @param data : tuple
    See parse_record() for the format

  @return: None
  """
  TS,means,variances,skews,kurts = data
  print(timestamp_to_str_with_ms(TS))
  print("Mean:     %8.5f, %8.5f (should be 0)" % means)
  print("Variance: %8.3f, %8.3f (power)" % variances)
  print("Skewness: %8.5f, %8.5f (should be zero)" % skews)
  print("Kurtosis: %8.5f, %8.5f" % kurts)

def get_block_time(fd,TS0,orate,blk_index):
  """
  This gets the time of a particular STATS binary file record.
    
  Notes
  =====
  Not currently used and probably doesn't work anymore.

  @param fd : file descriptor.

  @param TS0 : float.
    UNIX timestamp for the record

  @param orate : int.
    Averages per second

  @param blk_index : int.
    Index to the record at the start of a block
  """
  first_rec_index = blk_index*orate
  TS,mean,variance,skewness,kurtosis \
      = parse_record(TS0,orate,get_record(fd,first_rec_index))
  return TS
