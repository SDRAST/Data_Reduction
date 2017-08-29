"""
Get data from the .scr file

The relevant commands in the .scr file are::
  DDCLO       - selects the reatime signal processor (RSP) downconversion
                filter. Allowed values are AUTO, 80, 160, 240, 320, 400, 480,
                and 560. AUTO specifies that the predicts in current use are to
                be used to set the LO so that the channel filters are close to 
                centered in the 160 MHz pass band.
  FROV        - FRequency OVerride allows the user to override the sky
                frequency predicts with a fixed frequency. If the value is 0,
                then the predicts will be used instead.
  RF_TO_IF_LO - the RF to IF LO used to calculate baseband frequencies.
                If value is left at 0 the software sets it based on the 
                DOWNLINK_BAND from predicts and the following values::
                  L = 1200,
                  S = 2000,
                  X = 8100,
                  Ka = 31700 MHz 
  SFRO        - sets the specified WVSR channel frequency offset to the
                specified frequency value. The offset is relative to the
                predicts carrier plus any FRO frequency offset that has been
                specified, including any rate. The values  must be real with 
                optional K, M or G appended for kHz, MHz or GHz respectively.
"""
import logging

from glob import glob

logger = logging.getLogger(__name__)

def get_last_scr_file(obsdir, pol):
  """
  """
  logger.debug("get_last_scr_file: looking in %s", obsdir)
  files = glob(obsdir+"*"+pol+".scr")
  logger.debug("get_last_scr_file: found %s", files)
  if files == []:
    return None
  else:
    files.sort()
    return files[-1]  

def get_frequency(text):
  """
  """
  if text[-1].isdigit():
    freq = float(text[:-1])
  elif text[-1] == "M":
    freq = float(text[:-1])*1e6
  else:
    raise RuntimeError("Unknown FROV suffix %s", text[-1])
  return freq
  
def get_WVSR_parameters(scrdata):
  """
  """  
  bandwidth = {}
  bits_per_sample = {}
  sfro = {}
  for line in scrdata:
    parts = line.strip().split()
    if parts == []:
      continue
    if parts[0].upper() == "RF_TO_IF_LO":
      # the RF to IF LO used to calculate baseband frequencies.
      rf_to_if_lo = parts[-1]
      continue
    elif parts[0].upper() == "FROV":
      # sky frequency to be used instead of what is in the predicts
      freq_override = get_frequency(parts[-1])
      continue
    elif parts[0].upper() == "SFRO":
      # sets the specified WVSR channel frequency offset to the specified
      # frequency value.
      chan = int(parts[1])
      sfro[chan] = get_frequency(parts[-1])
      continue
    elif parts[0].upper() == "CHAN":
      # Allocates a channel to this client and configures its bandwidths and
      # bitrates 
      if parts[1].upper() == "ALL":
        continue
      chan = int(parts[1])
      bandwidth[chan] = get_frequency(parts[2])
      bits_per_sample[chan] = int(parts[3])
    elif parts[0].upper() == "DDCLO":
      chan = int(parts[1])
  ch_freq = {}
  if int(rf_to_if_lo)*1e6 != freq_override:
    logger.warning("get_WVSR_parameters: RF_TO_IF_LO does not match FROV")
  for chan in sfro.keys():
    ch_freq[chan] = freq_override + sfro[chan]
  return ch_freq, bandwidth, bits_per_sample

