# -*- coding: utf-8 -*-
"""
Handles old Tidbinbilla data files

NEEDS MAJOR UPDATES FROM MODULE Observatory TO MonitorControl
"""
import Astronomy as A
from time import strftime
import Physics as P
#import Observatory as O
import Data_Reduction.FITS.SDFITS as SDF

diag = False

class Tid_scan():
  """
  Make a Tidbinbilla spectral scan object

  This covers spectra taken with the HP9825 controller.  The header
  variables are all flots and the assigments are::

    01   file number of next scan to be stored
    02   integration time (secs)
    03   year (i.e. 1776)
    04   day of year
    05   UT (hrs) for the start of the scan
    06   current R.A. (hrs)
    07   current Decl. (deg)
    08   R.A. 1950 (hrs)
    09   Decl. 1950 (deg)
    10   Vlsr (km/s)
    11   Longitude (deg)
    12   Latitude (deg)
    13   doppler factor
    14   L.S.T. (hrs) for the end of the scan
    15   line rest freq (MHz)
    16   synthesizer (MHz)
    17   signal I.F. (MHz)
    18   phase-lock I.F. (MHz)
    19   sideband number
    20   calibration temp (K)
    21   time per on or off source integration (sec)
    22   observing mode number
      =1   TEST mode
                data                      -> C[] -> tape
      =2   POSITION SWITCHING (front end not switching: (-) beam)
                off source noise diode cal data -> K[]
                off source                      -> C[]
                on source                       -> D[]
                     Tcal*(D[]-C[])/(K[]-C[]) -> tape
      =3   LOAD SWITCHING (switching between (-) beam and load)
           TOTAL POWER (not switching: (-) beam (+) ref beam)
                off source                 -> C[]
                on source                 -> D[]
                     Tload*(D[]/C[]-1)    -> tape
                     Tsys*(D[]/C[]-1)     -> tape
      =4   CHOPPING (symmetrical switching between two horns)
                noise diode cal in (-) beam   -> K[]
                source in (+) beam            -> C[]
                source in (-) beam            -> D[]
                     0.5*Tcal*(D[]-C[])/K[]   -> tape
      =5   BEAM SWITCHING (switching between two horns)
                cal in (-) beam           -> K[]
                source in (-) beam        -> C[]
                     Tcal*C[]/K[]         -> tape
      =6   FREQUENCY SWITCHING (source in (-) beam)
                cal  with r35 ch. offset        -> K[]
                data with r35 ch. offset        -> C[]
                data with 0   ch. offset        -> D[]
                     Tcal*(D[]-C[])/(K[]-C[])   -> tape
      =21  POSITION SWITCHING (front end not switching: (-) beam)
                off source noise diode switched -> K[]
                off source                    -> C[]
                on source                   -> D[]
                     Tcal*(D[]-C[])/K[]      -> tape
      =99  Y-Factor measured with power meter
    23   number of channels of center frequency offset
    24   bandwidth (MHz)
    25   SDS total counts/sec
    26   SDS d.c. counts/sec
    27   system temperature (K)
    28   if 1, use comet velocity formula
    29   azimuth pointing offset
    30   elevation pointing offset
    31   H.A. reference offset
    32   dec. reference offset
    33   H.A. source offset
    34   dec. source offset
    35   channel offset for frequency switching
    36   synthesizer harmonic
    38-49 source name in ASCII codes
    50   number of channels (64: Y-factor, 256: SDS or 16384: DAVOS)
    51   DSS
    52   calibration time
    53   test tone frequency offset (relative to receiver freq)
    54   test tone level (dBm)

  Notes
  =====

  Receiver configurations
  -----------------------
  1. A K-band maser was installed in 1981. It had one linear
     polarization. The receiver was based on a phase-lock Gunn of an OVRO
     down-converting to S-band design. The back-end was a 256 channel
     SpectraData Fourier Transform Spectrometer built by George Morris.

  2. In early 1988 a Ku-band HEMT LNA front end was installed.  It
     used an HP synthesizer as a first LO, down-converting to S-band.

  3. Around 1995 an 18-26 GHz HEMT was installed.  Also, at some point,
     SpectraData was replaced with the DAVOS spectrometer.

  Frequency calculations
  ----------------------
  1. RESTFREQ is the frequency of interest, usually a line frequency,
     in the rest frame of the source, in Hz

  2. obs_freq is the frequency in the rest frame of the observer where
     this feature of interest appears.  It is computed from the Doppler
     factor.

  3. exp_freq is the frequency in the rest frame of the observer where
     this feature of interest should appear, calculated as an offset
     from the receiver frequency.  We expect exp_freq = obs_freq.

  Source name encoding
  --------------------
  In 1989, the convention was to put the source name on the first
  line of the scan block, after the scan number.

  1991/123 the header numbers[37:49] were reserved for source name
  encoding but filled with the negative of the index number.

  1991/174-194,1992/144-1145 zero-filled again.

  1992/165 source name encoding started but it was wrong, mostly::
    ['11 09 45.', '11 09 45.', '17 42 29.3', '17 42 29.3',....
     '17 42 29.3', 'G1.6 meth peak', 'G1.6 meth peak', 'G1.6 meth peak',
     '17 42 29.3',.... '17 42 29.3', '1', '1', '1']
  That was sorted out on 1992/236.

  1995/066 the source name encoding in numbers[37:49] was erroneously
  broken up by the number of spectrometer channels in numbers[45].
  Also, the number encoded was '256' whereas the DAVOS spectrometer
  (16384 channels) was already being used.

  On 1995/090, number of scans moved to numbers[49]

  If a source name cannot be extracted, then the coordinates are used
  to try to find the source name.
  """
  def __init__(self,dss,numbers):
    """
    Create an instance of the Tid_scan() class
    """
    self.header = {}
    self.get_scan_data(dss,numbers)

  def get_scan_data(self,dss,numbers):
    """
    Get data for a scan.

    This sets values in self.header and self.data
    
    @param dss : DSN station number
    @type  dss : int

    @param numbers : scan data
    @type  numbers : float

    @return: None
    """
    if diag:
      print "================================================"
    self.header['SCAN'] = int(numbers[0])        # header variable 1
    self.header['EXPOSURE'] = float(numbers[1])  # header variable 2
    year = int(numbers[2])
    doy = int(numbers[3])
    ut = float(numbers[4])
    yr,month,day = A.calendar_date(year,doy)
    hh = int(ut)
    mm = int((ut-hh)*60)
    ss = ((ut-hh)*60-mm)*60
    self.header['DATE-OBS'] = strftime('%Y-%m-%dT%H:%M:%S',
                                       (yr,month,day,hh,mm,ss,0,0,0))
    self.header['ra'] = float(numbers[5])        # decimal hours
    self.header['dec'] = float(numbers[6])       # decimal degrees
    self.header['CRVAL2'] = float(numbers[7])*15 # decimal deg B1950
    self.header['CRVAL3'] = float(numbers[8])    # decimal degrees B1950
    self.header['EQUINOX'] = 1950.0
    self.header['VELOCITY'] = float(numbers[9])  # source in the LSR
    # These station coordinate numbers are probably not as accurate as
    # what is in the database
    self.header['TELESCOP'] = "DSS "+str(dss)
    station_data = A.get_geodetic_coords()
    dss_long,dss_lat,dss_elev,tz,name, diam = station_data[dss]
    self.header['SITELONG'] = dss_long
    self.header['SITELAT'] = dss_lat
    self.header['SITEELEV'] = dss_elev
    # This is the Doppler factor computed by the HP9825
    # which is not as accurate as the full Manchester-Gordon code
    dopfac = float(numbers[12])                  # header variable 13
    if diag:
      print "HP9825 dopfac =",dopfac
      print "HP 9825 computed Doppler shift =",(1-dopfac)*P.c/1000,"km/s"
    # The Doppler calculation is probably not as good as this one
    # as it was Manchester Gordon FORTRAN converted to HP9825 language
    lst,vdop = \
      A.ManGord.dop(self.header['ra']*P.pi/12,
                    self.header['dec']*P.pi/180,
                    year, doy+ut/24.,
                    self.header['SITELAT'],
                    self.header['SITELONG'],
                    self.header['SITEELEV'])
    if diag:
      print "Manchester-Gordon computed Doppler shift =", \
            vdop+self.header['VELOCITY'],"km/s"
    MG_dopfac = (1-(vdop+self.header['VELOCITY'])/(P.c/1000))
    if diag:
      print "Manchester-Gordon computed dopfac =",MG_dopfac
    # It appears that the recorded LST is for the start of the scan
    # and the recorded UT for the end.
    if diag:
      print "Manchester-Gordon LST =",lst*24
    self.header['LST'] = lst*24 # float(numbers[13])
    HA = self.header['LST'] - self.header['ra']
    self.header['AZIMUTH'],self.header['ELEVATIO'] = \
      A.HaDec_to_AzEl(HA, self.header['dec'], self.header['SITELAT'])
    # Frequency of interest in the rest frame of the source
    # There is probably at least one digit missing off the end
    # of this number when the data are in ASCII:
    # 13778.80 instead of  Lovas 13778.804
    self.header['RESTFREQ'] = float(numbers[14])*1e6
    if diag:
      print "Frequency of interest in the source frame =", \
        self.header['RESTFREQ'],"Hz"
    # The observing frequency is what was calculated with the HP9825
    # Doppler factor calculation
    obs_freq = dopfac*self.header['RESTFREQ'] # Hz
    if diag:
      print "Frequency of interest in the observer frame =",obs_freq,"Hz"
      print "corresponding velocity =",(1-MG_dopfac)*P.c/1000,"km/s"
    # Generally, the frequency can be calculated from these numbers
    self.header['synth'] = float(numbers[15])
    self.header['IFsig'] = float(numbers[16]) # MHz
    self.header['IFpll'] = float(numbers[17]) # MHz
    if int(numbers[18]) == 0:
      self.header['SIDEBAND'] = 'U'
    else:
      self.header['SIDEBAND'] = 'L'
    self.header['TCAL'] = float(numbers[19])
    self.header['t_onoff'] = float(numbers[20])
    obs_mode = int(numbers[21])
    if obs_mode == 1:
      self.header['OBSMODE'] = 'LINETLPW'
    elif obs_mode == 2 or obs_mode == 21:
      self.header['OBSMODE'] = 'LINEPSSW'
    elif obs_mode == 3:
      self.header['OBSMODE'] = 'LINELDSW'
    elif obs_mode == 4:
      self.header['OBSMODE'] = 'LINECHOP'
    elif obs_mode == 5:
      self.header['OBSMODE'] = 'LINEBMSW'
    elif obs_mode == 6:
      self.header['OBSMODE'] = 'LINEFQSW'
    elif obs_mode == 99:
      self.header['OBSMODE'] = 'CONTLDSW'
    else:
      self.header['OBSMODE'] = 'UNKNOWN'
    chan_ofst = int(numbers[22])
    self.header['BANDWIDT'] = float(numbers[23])*1.e6 # Hz
    SDStotal = int(numbers[24])
    SDSdc = int(numbers[25])
    self.header['TSYS'] = float(numbers[26])
    comet = int(numbers[27])
    self.header['PTG_AOFF'] = float(numbers[28])
    self.header['PTG_EOFF'] = float(numbers[29])
    self.header['REF_HOFF'] = float(numbers[30])
    self.header['REF_DOFF'] = float(numbers[31])
    self.header['BEAMHOFF'] = float(numbers[32])
    self.header['BEAMDOFF'] = float(numbers[33])
    harmonic = int(numbers[35])
    self.header['MAXIS1'] = number_of_channels(numbers)
    unknown = float(numbers[36])
    # Ignore if name header data has negative integers
    if float(numbers[37]) > 31.0:
      # OK; looks like an ASCII code
      # Source name is in the header?
      name = ""
      for number in numbers[37:49]:
        if number == numbers[45] and float(numbers[45]) > 256.0:
          # number of channels in the wrong place
          name += ''
        else:
          name += chr(int(number))
      self.header['OBJECT'] = name.strip()
      # Now the source names may have some strange errors which
      # we'll clean up here.  Leading colon?
      self.header['OBJECT'] = self.header['OBJECT'].lstrip(":")
      # LaTeX formatting doesb't like this
      self.header['OBJECT'] = self.header['OBJECT'].replace('_',' ')
    else:
      # Once the coordinates are known we can try to find the source name
      self.header['OBJECT'] = A.find_source_name(self.header['CRVAL2'],
                                                  self.header['CRVAL3'])
    self.header['CDELT1'] = self.header['BANDWIDT']/self.header['MAXIS1'] # Hz
    # This assumes no tapering. This is the full width at half power
    # for a sinc function.
    self.header['FREQRES'] = 2*0.443*self.header['CDELT1'] # Hz
    # Calculate the receiver center frequency in MHz.  It could be zero
    # because some factor was zero in the header
    recvfreq = self.header['synth'] * harmonic \
             + (1-int(numbers[18])) * self.header['IFsig'] # MHz
    freq_offset = chan_ofst * self.header['CDELT1']/1.e6   # MHz
    if diag:
      print "Frequency offset =",freq_offset,"MHz"
    # Center frequency during reference phase of frequency switching
    self.header['FOFFREF1'] = float(numbers[34])*self.header['CDELT1'] # Hz
    if recvfreq != 0.0:
      # Receiver frequency calculation succeeded
      self.header['OBSFREQ'] = recvfreq*1.e6 # Hz
      exp_freq = (recvfreq + freq_offset)*1.e6      # Hz
      if diag:
        print "Actual receiver center frequency =",recvfreq,"MHz"
        print "Frequency where line is expected =",exp_freq,"MHz"
    else:
      # Receiver frequency calculation failed
      # Assume that the line is where it is expected to be
      # that is, the frequency of interest - offset
      # (Remember that obs_freq is in Hz)
      recvfreq = obs_freq - freq_offset*1.e6  # Hz
      self.header['OBSFREQ'] = recvfreq # Hz
      exp_freq = obs_freq # Hz
      if diag:
        print "Putative receiver center frequency =",recvfreq,"MHz"
    self.header['CRVAL1'] = exp_freq
    # Frequency offset in MHz from spectrum center
    self.header['CRPIX1'] = self.header['MAXIS1']/2 - chan_ofst
    # Doppler shifted frequency:
    #    f          v
    # ------ =  1 - -
    # f_rest        c
    # Doppler velocity:
    # v         f
    # - = 1 - ------
    # c       f_rest
    # Total Doppler shift between source and observatory
    v_doppler = (1 - (obs_freq/self.header['RESTFREQ'])) * (P.c/1000)
    self.header['RVSYS'] = v_doppler
    v_obs = v_doppler - self.header['VELOCITY']
    self.header['VFRAME'] = v_obs
    if diag:
      print "Velocity offset =",v_doppler,"km/s"
      print "Velocity of source w.r.t. LSR =",self.header['VELOCITY'],"km/s"
      print "Velocity of the observatory w.r.t. LSR =",v_obs,"km/s"
    # Make these up, sort of
    self.header['TEMPSCAL'] = 'TA'
    if self.header['OBSFREQ']/1e6 > 10000. and \
       self.header['OBSFREQ']/1e6 < 17000.:
      self.header['FRONTEND'] = 'Ku HEMT'
      self.header['CRVAL4'] = -2 # LCP ?
    elif self.header['OBSFREQ']/1e6 > 18000. and year < 1993:
      self.header['FRONTEND'] = 'K MASER'
      self.header['CRVAL4'] = -2 # LCP ?
    elif self.header['OBSFREQ']/1e6 < 10000.:
      self.header['FRONTEND'] = 'DSN X'
      self.header['CRVAL4'] = -1 # RCP ?
    else:
      self.header['FRONTEND'] = 'K HEMT'
      self.header['CRVAL4'] = -2 # LCP ?
    if self.header['MAXIS1'] == 256:
      self.header['BACKEND'] = 'SpectraData'
      self.data = numbers[64:320]
    elif self.header['MAXIS1'] == 16384:
      self.header['BACKEND'] = 'DAVOS'
      self.data = numbers[64:16448]
    else:
      self.header['BACKEND'] = 'unknown'
 

  def print_header(self):
    """
    Print the header of a Tid_scan instance

    @return: None
    """
    print "--------------------------------------------"
    for k in self.header.keys():
      print k,self.header[k]

def number_of_channels(numbers):
  """
  Deduce the number of channels from Tid ASCII file data

  This also handles the case when the number of channels was written
  in the middle of the source name field.

  Notes
  =====
  I wasn't thinking when I wrote this.  One only needs to
  check len(numbers)

  @param numbers : list of floats
    The data from the scan

  @return: int
  """
  # Number of channels in header?
  if int(numbers[45]) > 256:
    # number of channels in the wrong place
    return int(numbers[45])
  elif int(numbers[49]) > 256:
    # There is a number of channels count
    return int(numbers[49])
  else:
    # number of channels not in header
    return 256

def configure_station(dss):
  """
  Describe the station and report on the elements and their
  relationships.

  @param dss : int::
    station number

  @return: tuple::
    tid, tel, beam_routers, FEs, FE_channels, rf_switches, \
    downconverters, if_switches, backends
  """
  # Describe the antenna configuration
  tid, tel, beam_routers, FEs, FE_channels, rf_switches, downconverters, \
    if_switches, backends \
      = O.describe_antenna("Tidbinbilla", dss)

  # Report the on antenna configuration elements
  if diag:
    O.report_dictionaries(tid, tel, beam_routers, FEs, FE_channels,
                        rf_switches, downconverters, if_switches,
                        backends)
  return tid, tel, beam_routers, FEs, FE_channels, rf_switches, \
    downconverters, if_switches, backends

def update_backend(numbers,backends):
  """
  Sets the attributes 'polarizations' and 'num_chan'

  Sets these attrubutes for the backend used to obtain the data
  in 'numbers'
  """
  num_chans = number_of_channels(numbers)
  if num_chans == 256:
    backends['SpectraData'].polarizations = [SDF.pol_code('LCP')]
    backends['SpectraData'].num_chan = num_chans
  else:
    backends['DAVOS'].polarizations = [SDF.pol_code('LCP')]
    backends['DAVOS'].num_chan = num_chans
  return backends

def fill_bintable_header(tid_scan,num_backends,tbhead):
  """
  Initializes the scantable header

  Sets values for specific keys in the SDFITS scantable header.

  @param tid_scan : a Tid scan object
  @type  tid_scan : instance of Tid_scan()

  @param tbhead : scantable header
  @type  tbhead : dictionary

  @return: revised tbhead
  """
  tbhead.update('TELESCOP', tid_scan.header['TELESCOP'])
  tbhead.update('SITELONG', tid_scan.header['SITELONG'])
  tbhead.update('SITELAT',  tid_scan.header['SITELAT'])
  tbhead.update('SITEELEV', tid_scan.header['SITEELEV'])
  if num_backends == 1:
    tbhead.update('BANDWIDT', tid_scan.header['BANDWIDT'])
    tbhead.update('FREQRES',  tid_scan.header['FREQRES'])
  tbhead.update('FRONTEND', tid_scan.header['FRONTEND'])
  tbhead.update('TEMPSCAL', tid_scan.header['TEMPSCAL'])

def fill_bintable_row(tid_scan,num_backends,tbdata,table_row):
  """
  Writes a Tid scan to a scantable row

  @param tid_scan : a Tid scan
  @type  tid_scan : instance of Tid_scan()

  @param num_backends : 
  """
  scannum = tid_scan.header['SCAN']
  tbdata.field('SCAN')[table_row]       = scannum
  tbdata.field('CYCLE')[table_row]      = 0
  tbdata.field('DATE-OBS')[table_row]   = tid_scan.header['DATE-OBS']
  tbdata.field('OBJECT')[table_row]     = tid_scan.header['OBJECT']
  tbdata.field('OBSMODE')[table_row]    = tid_scan.header['OBSMODE']
  tbdata.field('RESTFREQ')[table_row]   = tid_scan.header['RESTFREQ']
  tbdata.field('MAXIS1')[table_row]     = tid_scan.header['MAXIS1']
  if num_backends > 1:
    # Otherwise it goes in the table header
    tbdata.field('BANDWIDT')[table_row] = tid_scan.header['BANDWIDT']
    tbdata.field('FREQRES')[table_row]  = tid_scan.header['FREQRES']
  tbdata.field('RVSYS')[table_row]      = tid_scan.header['RVSYS']
  tbdata.field('VFRAME')[table_row]     = tid_scan.header['VFRAME']
  tbdata.field('VELOCITY')[table_row]   = tid_scan.header['VELOCITY']
  tbdata.field('EXPOSURE')[table_row]   = tid_scan.header['EXPOSURE']
  tbdata.field('OBSFREQ')[table_row]    = tid_scan.header['OBSFREQ']
  tbdata.field('CRVAL1')[table_row]     = tid_scan.header['CRVAL1']
  tbdata.field('CDELT1')[table_row]     = tid_scan.header['CDELT1']
  tbdata.field('CRPIX1')[table_row]     = tid_scan.header['CRPIX1']
  tbdata.field('CRVAL2')[table_row]     = tid_scan.header['CRVAL2']
  tbdata.field('CRVAL3')[table_row]     = tid_scan.header['CRVAL3']
  tbdata.field('CRVAL4')[table_row]     = tid_scan.header['CRVAL4']
  tbdata.field('EQUINOX')[table_row]    = tid_scan.header['EQUINOX']
  tbdata.field('AZIMUTH')[table_row]    = tid_scan.header['AZIMUTH']
  tbdata.field('ELEVATIO')[table_row]   = tid_scan.header['ELEVATIO']
  tbdata.field('TSYS')[table_row]       = tid_scan.header['TSYS']
  tbdata.field('DATA')[table_row]       = tid_scan.data
  tbdata.field('REF_HOFF')[table_row]   = tid_scan.header['REF_HOFF']
  tbdata.field('REF_DOFF')[table_row]   = tid_scan.header['REF_DOFF']
  tbdata.field('BEAMHOFF')[table_row]   = tid_scan.header['BEAMHOFF']
  tbdata.field('BEAMDOFF')[table_row]   = tid_scan.header['BEAMDOFF']
  tbdata.field('FOFFREF1')[table_row]   = tid_scan.header['FOFFREF1']
  table_row += 1
  return tbdata, table_row
