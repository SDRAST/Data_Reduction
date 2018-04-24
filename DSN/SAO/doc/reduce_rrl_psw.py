# Python packages
import warnings
import re
import logging
import sys

from pylab import *
from os import chdir, mkdir, makedirs
from os.path import basename, curdir, dirname, exists
from scipy import interpolate
from scipy import signal
from optparse import OptionParser

# Third part packages
from novas import compat as novas
from novas.compat import eph_manager
jd_start, jd_end, number = eph_manager.ephem_open()
import ephem
import dill as cPickle


# Local packages
from Astronomy import B_epoch_to_J, MJD, v_sun
from Astronomy import c # m/s
from Astronomy.redshift import doppler_radio, doppler_optical, doppler_relat
from DatesTimes import UnixTime_to_datetime, UnixTime_to_MPL
from MonitorControl.BackEnds.ROACH1.SAOspec import SAOhdf5
from MonitorControl.Configurations.coordinates import DSS
from support import nearest_index
from support.logs import initiate_option_parser, init_logging #, loglevel
from support.text import clean_TeX, select_files 

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)
 
def compute_X_axis(dataset, frame, dss, ref_freq=None, vspline=None, 
                    time=None):
   """
   Computes the appropriate X-axis for the averaged difference spectrum
    
   'vobj' is velocity of the object in the rest velocity of frame specified.
   Acceptable frames are defined in the SAOhdf5.rel_freq_units() docstring.
   In addition we allow here RADIO-OBJ which defines the rest frame of the
   object.
   
   @param header : information about the target source
   @type  header : dict
   
   @param frame : the rest frame and the X-axis type (freq or vel)
   @type  frame : str
   
   @param dss : DSS station
   @type  dss : int
   
   @param ref_freq : frequency in MHz for computing relative frequencies
   @type  ref_freq : float
   
   @param vspline : radial velocity of a moving body as a function of time
   @type  vspline : function
   
   @param time : the time at which the spline is to ve evaluated
   @type  time : UNIX seconds
   """
   n_chans = dataset.data[1].shape[0]
   if ref_freq:
     f_ref = ref_freq
   else:
     f_ref = dataset.header['linefreq']/1e6 # MHz
   v_ref = dataset.header['VELOCITY']
   logger.debug(" compute_X-axis: requested frame is %s", frame)
   logger.debug(" compute_X_axis: reference frequency is %10.3f", f_ref)
   if frame == "CHAN-OBS" or frame == "FREQ-OBS" or frame == "RELA-OBS":
     x = dataset.rel_freq_units(frame=frame, ref_freq=f_ref)
     vobj = None
   elif frame == "RADI-OBS":
     vobj = V_LSR(dataset.header, time, dss)
     x = dataset.rel_freq_units(frame=frame, ref_freq=f_ref)
   elif frame == "RADI-LSR":
     x = dataset.rel_freq_units(frame=frame, ref_freq=f_ref, 
                             v_frame=V_LSR(dataset.header, time, dss))
     vobj = v_ref
     logger.debug("compute_X_axis: vobj = %.2f", vobj)
   elif frame == "RADI-OBJ":
     # This is the object's rest frame
     if vspline and time:
       vobj = vspline(time)
       x = -(c/1000)*dataset.rel_freq_units(frame="DELF-OBS")/f_ref - vobj
     else:
       vobj = dataset.header[0]['VELOCITY']
       x = dataset.rel_freq_units(frame=frame, ref_freq=f_ref,
                               v_frame=V_LSR(dataset.header, time, dss) + vobj)
   else:
     self.logger.warning(" frame %s is not valid", frame)
     return
   return x, frame, vobj
  
def V_LSR(header, dt, dss):
  """
  Computes the velocity of the local standard of rest w.r.t. the observer
  
  @param header : source data
  @type  header : dict
  
  @param dt : date/time of the observation
  @type  dt : datetime object
  
  @param dss : DSN station
  @type  dss : int
  """
  ra2000,dec2000 = B_epoch_to_J(header['R.A.'][0],
                                header['declination'][0],
                                'decimal')
  cat_entry = novas.make_cat_entry(header["OBJECT"][0],
                                   "",0,ra2000, dec2000, 0, 0, 0, 0)
  source = novas.make_object(2, 0, header["OBJECT"][0], cat_entry)
  dss43 = DSS(dss)
  DSS43 = novas.make_observer_on_surface(dss43.lat*180/pi,
                                         360+dss43.lon*180/pi,
                                         dss43.elev, 0, 0)
  jd = novas.julian_date(dt.year,dt.month,dt.day,dt.hour+dt.minute/60.)
  mjd = MJD(dt.year,dt.month,dt.day)
  earth = novas.make_object(0, 3, 'Earth', None)
  urthpos,urthvel = novas.ephemeris((jd,0), earth, origin=0)
  (obspos,obsvel) = novas.geo_posvel(jd,0,DSS43,0)
  totvel = tuple(array(urthvel)+array(obsvel))
  (srcpos,srcvel) = novas.starvectors(cat_entry)
  V = novas.rad_vel(source, srcpos, srcvel, totvel,0,0,0)
  logger.debug("V_LSR: velocity of LSR = %.2f", V)
  return V+v_sun(mjd,ra2000,dec2000)

def average_calibrated_spectrum(self):
    """
    Averages all the difference spectra in a scan.
    
    Because this works on scan pairs, the number of pairs is half the number of
    scans in the dataset.
    """
    scans = array(self.data.keys())
    self.logger.debug("average_calibrated_spectrum: for scans %s", scans)
    n_chans = self.data[1].shape[0]
    average_spectra = {0: zeros((n_chans)), 1: zeros((n_chans))}
    total_records = 0
    integr = {0: 0.0, 1: 0.0}
    if len(scans) % 2:
      # truncate odd number of scans
      scans = scans[:-1]
      self.logger.warning("average_calibrated_spectrum: odd scan removed")
      self.logger.debug("average_calibrated_spectrum: doing %s", scans)
    for scan in scans[::2]:
      self.logger.debug("average_calibrated_spectrum: processing scan pair %d",
                        scan)
      spectrum, n_recs, tau = calibrated_scan_pair_difference(self,scan)
      for pol in [0,1]:
        average_spectra[pol] += n_recs*spectrum[pol]
        total_records += n_recs
        integr[pol] += tau[pol]
    for pol in [0,1]:
      average_spectra[pol] /= (total_records/2)
    return average_spectra, integr

def calibrated_scan_pair_difference(self,scan):
    """
    Computes a calibrated spectrum from a pair of scans.
    
    The first scan number must be odd with source in beam 1.  The associated
    scan with the source in baem 2 will be 'scan'+1.
    
    Returns a dict with calibrated spectrum for each pol and the number of
    records processed. Note that the number of differences in the average
    is n_recs/2
    
    @param scan : scan number, starting with 1, not 0
    @type  scan : int
    
    @return: tuple of dict,int,float
    """
    self.logger.debug("calibrated_scan_pair_difference: for scan %d", scan)
    scan_index = scan-self.first_scan
    n_chans = self.data[scan].shape[0]
    n_recs = min(array(self.header['num_records'][scan].values()).min(),
                 array(self.header['num_records'][scan+1].values()).min())
    
    ratio_spectrum = {}
    tau = {}
    for pol in [0,1]:
      tau[pol] = 0.0
      ratio_spectrum[pol] = zeros((n_chans))
      self.logger.debug(
        "calibrated_scan_pair_difference: processing %d records for pol %d",
        pol, n_recs)
      for record in range(n_recs):
        on1data  = self.data[scan][:,0,0,pol,record,0]
        off1data = self.data[scan][:,0,0,pol,record,1]
        on2data  = self.data[scan+1][:,0,0,pol,record,1]
        off2data = self.data[scan+1][:,0,0,pol,record,0]
        try:
          ratio1 = on1data/off2data
        except RuntimeWarning,details:
          self.logger.warning(
            "calibrated_scan_pair_difference: "+
            "ratio1 divide by zero in record %d; %s",
            record, str(details))
        try:
          ratio2 = on2data/off1data
        except RuntimeWarning,details:
          self.logger.warning(
            "calibrated_scan_pair_difference: "+
            "ratio2 divide by zero in record %d: %s",
            record, str(details))
        on1Tsys  = self.header['TSYS'][scan][pol,record,0]
        off1Tsys = self.header['TSYS'][scan][pol,record,1]
        on2Tsys  = self.header['TSYS'][scan][pol,record,1]
        off2Tsys = self.header['TSYS'][scan][pol,record,0]
#        on1Tsys  = self.header['TSYS'][pol,record,0]
#        off1Tsys = self.header['TSYS'][pol,record,1]
#        on2Tsys  = self.header['TSYS'][pol,record,1]
#        off2Tsys = self.header['TSYS'][pol,record,0]
        try:
          ratio_spectrum[pol] += (ratio1 - 1)*off2Tsys + (ratio2 - 1)*off1Tsys
        except UnboundLocalError, details:
          self.logger.error("calibrated_scan_pair_difference: %s", str(details))
          continue
      tau[pol] += self.header['EXPOSURE'][pol]
      tau[pol] += self.header['EXPOSURE'][pol]
      ratio_spectrum[pol] /= n_recs
    return ratio_spectrum, n_recs, tau

if __name__ == "__main__":
  ########################
  #Ask user for file names:

  yearpath = select_files("/usr/local/RA_data/HDF5/dss43/*",
                               text="Select a year by index: ", single=True)
  year = basename(yearpath)

  daypath = select_files("/usr/local/RA_data/HDF5/dss43/"+year+"/*",
                               text="Select a day by index: ", single=True)
  day = basename(daypath)
 
  datafiles=select_files("*.pkl")
  #datafiles=select_files("/usr/local/RA_data/HDF5/dss43/"+year+"/"+day+"/*.h5")
  ########################
  for datafile in datafiles:

    print "Processing ",datafile
  
    filename=basename(datafile)
    #######################
    #Read data set

    #  hdf = SAOhdf5(datafile)
    #  data = hdf.to_dataset('/home/jpineda/PDR3/') # loads the data into the dataset
    #  data.save_pickle()

    ## process a pickled dataset
    fd = open(datafile, "rb")
    data = cPickle.load(fd)
    fd.close() 

    #################
    #update for X-band observations
    #  data.header['band']=8.1


    ####################################
    average_spec, tau = average_calibrated_spectrum(data)
    ####################################
    ############################################################
    index=0

    meantime = UnixTime_to_datetime((data.header['start']
                                    +data.header['end'])/2)
    frame = "RADI-LSR"
    #sidereal source
    spl = None
    #freq=8584.82 #H92alpha
    #freq=8309.38 #H91alpha
    freq=22364.17 #H67a 
    source='LMC12_SE'
    line='H67a'

    x, data.header['frame'], data.header['velocity'] = compute_X_axis(data,
                                                         frame, 43,
                                                         vspline=spl,
                                                         ref_freq=freq,
                                                         time=meantime)

    #  plot(x,average_spec)
    #  show()

    of = open(filename+line+"_pol1.asc","w")
  
    for i in range(0,len(average_spec[0])-1):
      print >> of, ("%10.3f  %8.4f " % (x[i],average_spec[0][i]))

    index=0

    of = open(filename+line+"_pol2.asc","w")
  
    for i in range(0,len(average_spec[1])-1):
      print >> of, ("%10.3f  %8.4f " % (x[i],average_spec[1][i]))

