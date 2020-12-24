import logging
logging.basicConfig(level=logging.DEBUG)
mylogger = logging.getLogger()
mylogger.setLevel(logging.DEBUG)
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)
plogger = logging.getLogger("parso")
plogger.setLevel(logging.WARNING)
alogger = logging.getLogger("Astronomy")
alogger.setLevel(logging.WARNING)

import warnings
warnings.filterwarnings("ignore")

from pylab import *

from Data_Reduction.GAVRT.plotter import SessionPlotter

print("\nGetting session plotter")
sp = SessionPlotter(None, 2020, 347)
mpl = sp.get_map_plotters()
sp.list_maps(save=True)
sp.make_bs_dir(save=True)
sp.summary()
print("\nMaking raster summary")
sp.raster_summary()

