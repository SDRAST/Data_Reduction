from matplotlib.pylab import date2num


def plot_tsys(fig, date_nums, tsys, label=None, picker=3):
  """
  Plot system temperatures

  @param fig : figure which contains current axes
  @type  fig : figure() instance

  @param date_nums : datetime for each point
  @type  date_nums : list of mpl date numbers

  @param tsys : system temperature counts
  @type  tsys : list of float
  """
  #plot_date(date_nums, tsys, '-')
  plot_date(date_nums, tsys, '-', picker=picker)
  grid()
  text(1.05,0.5,label,
       horizontalalignment='left',
       verticalalignment='center',
       transform = gca().transAxes)
  title = num2date(date_nums[0]).strftime("%Y %j")
  fig.autofmt_xdate()
  return title

