# -*- coding: utf-8 -*-
"""
hyperfine_fit -  Fit Gaussian lines to hyperfine components

Functions needed to fit data::
  ammonia_data()
  model()
  error_func()

Here is the basic procedure, assuming that the data are in arrays
x(in MHz) and y::
  plot(x,y)                          # to see what to expect
  set_model(line)                    # see ammonia_data
  parameter_guess = [.3,-0.6,width]  # Put in appropriate numbers
  result = fit_data(x,y,hpf_model,parameter_guess)
  best_fit = hpf_model(x,*result)    # model y-values
  plot(x,best_fit)                   # overlaid on the data

Problem
=======
Getting the ammonia hyperfine structure from the JPL catalog (v4)
and calculating A_ul gives the wrong relative intensities.
10**logint gives the right ratios.
"""
from Math import gaussian
from numpy import array, linspace, vectorize, zeros
Gaussian = vectorize(gaussian)

from pylab import draw, ion, title, plot, xlabel, ylabel
# This overrides the matplotlib.mlab function norm()
from scipy.stats.distributions import norm
from scipy.optimize import leastsq
#import time

import Physics.Radiation.Lines.Molec as M

def hpf_model(x, amp, x_0, w):
  """
  Builds multicomponent Gaussian profile

  @param x : numpy array of float::
    x-values

  @param amp : float::
    scale factor for line strenths

  @param x_0 : float::
    center position of the model

  @param w : float::
    width of each Gaussian component

  @return: numpy array of float::
    y-values for the model
  """
  global s, df
  result = zeros(len(x))
  for i in range(len(df)):
    result += amp*s[i]*Gaussian(x, x_0+df[i], w)
  return result

def error_func(guess, x, data, data_model):
  """
  @param guess : parameter list::
    an estimate of the model parameters of the final solution

  @param x : numpy array of float::
    x-values

  @param data : numpy array of y-values

  @param data_model : function to be fitted

  @return: numpy array of float::
    differences between the data and the model
  """
  return data - data_model(x, *guess)

def noise_data(SNR,width,span):
  """
  Create noisy data samples

  Noise is Gaussian (normal distribution) with a dispersion of 1.

  @param SNR : float::
    SNR of central component

  @param width : float::
    line width in MHz

  @param span : float::
    frquency range MHz for which data will be generated

  @return: numpy array of float::
    noisy y-values.
  """
  global nsamps
  x = linspace(-span/2,span/2,nsamps)
  noisy_data = hpf_model(x, SNR/s[0], 0, width) + norm.rvs(size=nsamps)
  return x, noisy_data

def ammonia_data_old(line):
  """
  hyperfine data for ammonia

  @type line : int
  @param line : line specifier
    1  - 23694.4955 MHz ammonia 1,1 line
    2  - 23722.6333 MHz ammonia 2,2
    -1 - 23694.4955 MHz ammonia 1,1 with very fine splitting

  @return: (name, freq_offset, rel_intensity)
  """
  if line == 1:
    # Parameters for 23694.4955 MHz ammonia 1,1 line
    name = r"$\mbox{NH}_3 ~~~1_1$"
    df = [-1.539, -0.599, 0.000,  0.599, 1.539] # frequency offsets
    s  = [ 0.110,  0.140, 0.500,  0.140, 0.110] # fraction of line strength
  elif line == -1:
    # Rydbeck et al, Ap.J. 215, L35, 1977
    name = r"$\mbox{NH}_3 ~~~1_1$"
    df = [-1.569, -1.525, -0.623, -0.587, -0.015, 0.010, 0.572, 0.615, 1.539]
    #ref [ 0.280,  0.480,  0.330,  0.560,  0.860, 1.000, 0.380, 0.520, 0.630]
    s  = [ 0.056,  0.095 , 0.065,  0.111,  0.171, 0.198, 0.075, 0.103, 0.125]
  elif line == 2:
    name = r"$\mbox{NH}_3 ~~~2_2$"
    # Parameters for 23722.6333 MHz ammonia 2,2
    df = [0.0, -1.315, 1.315, -2.046, 2.046] # frequency offsets
    s  = [0.796, 0.052, 0.052, 0.050, 0.050] # fraction of line strength
  return name,df,s

def ammonia_data(ju,ku,jl,kl):
  spec_lines = []
  nh3_file = open("/usr/local/line_cat/jpl/c017002.cat","r")
  data_lines = nh3_file.readlines()
  nh3_file.close()
  name, n_lines, part_fn = M.jpl.get_mol_metadata(17002)
  part_fn_300 = float(part_fn[300])
  first_line = True
  for line in data_lines[6408:6429]:
    #(ju,ku,V,
    line_data = M.jpl.parse_catalog_line(line)
    if first_line:
      print "Upper -> lower,",
      print M.jpl.quantum_label(17002,line_data['q lower'],
                                      line_data['qn format'],
                                      line_data['deg free'])[0]
      first_line = False
    freq   = float(line_data['freq'])
    intens = float(line_data['int'])
    g_up   = int(line_data['g upper'])
    E_lo   = float(line_data['E lower'])
    A_ul = M.jpl.einstein_a(intens, freq, part_fn_300, g_up, E_lo)
    print M.jpl.quantum_label(17002,line_data['q upper'],
                                    line_data['qn format'],
                                    line_data['deg free'])[1],
    print M.jpl.quantum_label(17002,line_data['q lower'],
                                    line_data['qn format'],
                                    line_data['deg free'])[1],
    print "%7.1f %5.2f" % ((freq-23694.4955)*1e3, (10.**intens)/6.4e-6)
    #print "%7.1f %5.2f" % ((freq-23694.4955)*1e3, A_ul/0.78e-7)
  return
  
def fit_data(x,data,model,parameter_guess):
  """
  fit noisy data to the model

  @type x : numpy array of float
  @param x : array of abcissa values

  @type data : numpy array of float
  @param data : array of ordinate values

  @type model : name of a function
  @param model : function to be fitted to data

  @type parameter_guess : list
  @param parameter_guess : [amplitude, center_x, width]

  @return: tuple
    output from scipy.optimize.leastsq
  """
  result, msg = leastsq(error_func,
                      x0 = parameter_guess,
                      args=(x,data,model))
  print msg
  return result

def set_model(line):
  """
  Make the model parameters global

  @type line : int
  @param line : line identifier
    see ammonia_data()

  @return: True if it worked
  """
  global name, df, s
  try:
    name, df, s = ammonia_data(line)
    return True
  except Exception,details:
    print Exception, details
    return False

def test():
  """
  Test with simulated and real data.
  """
  global nsamps
  width  = 0.1 # MHz
  ion()
  o = raw_input("Simulation (s) or real data (d)? ")
  if o.lower()[0] == 's':
    # make a set of noisy data samples
    SNR    = 6   # for central component
    nsamps = 200
    set_model(1)
    x, y = noise_data(SNR,width,10)
    parameter_guess = [1,0,width]
  else:
    x = array([-1.84357808, -1.83374336, -1.82390864, -1.81407392, -1.8042392 ,
       -1.79440448, -1.78456976, -1.77473504, -1.76490032, -1.7550656 ,
       -1.74523088, -1.73539616, -1.72556144, -1.71572672, -1.70589201,
       -1.69605729, -1.68622257, -1.67638785, -1.66655313, -1.65671841,
       -1.64688369, -1.63704897, -1.62721425, -1.61737953, -1.60754481,
       -1.59771009, -1.58787537, -1.57804065, -1.56820593, -1.55837121,
       -1.54853649, -1.53870177, -1.52886705, -1.51903233, -1.50919761,
       -1.49936289, -1.48952817, -1.47969345, -1.46985874, -1.46002402,
       -1.4501893 , -1.44035458, -1.43051986, -1.42068514, -1.41085042,
       -1.4010157 , -1.39118098, -1.38134626, -1.37151154, -1.36167682,
       -1.3518421 , -1.34200738, -1.33217266, -1.32233794, -1.31250322,
       -1.3026685 , -1.29283378, -1.28299906, -1.27316434, -1.26332962,
       -1.2534949 , -1.24366018, -1.23382546, -1.22399075, -1.21415603,
       -1.20432131, -1.19448659, -1.18465187, -1.17481715, -1.16498243,
       -1.15514771, -1.14531299, -1.13547827, -1.12564355, -1.11580883,
       -1.10597411, -1.09613939, -1.08630467, -1.07646995, -1.06663523,
       -1.05680051, -1.04696579, -1.03713107, -1.02729635, -1.01746163,
       -1.00762691, -0.99779219, -0.98795748, -0.97812276, -0.96828804,
       -0.95845332, -0.9486186 , -0.93878388, -0.92894916, -0.91911444,
       -0.90927972, -0.899445  , -0.88961028, -0.87977556, -0.86994084,
       -0.86010612, -0.8502714 , -0.84043668, -0.83060196, -0.82076724,
       -0.81093252, -0.8010978 , -0.79126308, -0.78142836, -0.77159364,
       -0.76175892, -0.75192421, -0.74208949, -0.73225477, -0.72242005,
       -0.71258533, -0.70275061, -0.69291589, -0.68308117, -0.67324645,
       -0.66341173, -0.65357701, -0.64374229, -0.63390757, -0.62407285,
       -0.61423813, -0.60440341, -0.59456869, -0.58473397, -0.57489925,
       -0.56506453, -0.55522981, -0.54539509, -0.53556037, -0.52572565,
       -0.51589093, -0.50605622, -0.4962215 , -0.48638678, -0.47655206,
       -0.46671734, -0.45688262, -0.4470479 , -0.43721318, -0.42737846,
       -0.41754374, -0.40770902, -0.3978743 , -0.38803958, -0.37820486,
       -0.36837014, -0.35853542, -0.3487007 , -0.33886598, -0.32903126,
       -0.31919654, -0.30936182, -0.2995271 , -0.28969238, -0.27985766,
       -0.27002295, -0.26018823, -0.25035351, -0.24051879, -0.23068407,
       -0.22084935, -0.21101463, -0.20117991, -0.19134519, -0.18151047,
       -0.17167575, -0.16184103, -0.15200631, -0.14217159, -0.13233687,
       -0.12250215, -0.11266743, -0.10283271, -0.09299799, -0.08316327,
       -0.07332855, -0.06349383, -0.05365911, -0.04382439, -0.03398968,
       -0.02415496, -0.01432024, -0.00448552,  0.0053492 ,  0.01518392,
        0.02501864,  0.03485336,  0.04468808,  0.0545228 ,  0.06435752,
        0.07419224,  0.08402696,  0.09386168,  0.1036964 ,  0.11353112,
        0.12336584,  0.13320056,  0.14303528,  0.15287   ,  0.16270472,
        0.17253944,  0.18237416,  0.19220888,  0.2020436 ,  0.21187831,
        0.22171303,  0.23154775,  0.24138247,  0.25121719,  0.26105191,
        0.27088663,  0.28072135,  0.29055607,  0.30039079,  0.31022551,
        0.32006023,  0.32989495,  0.33972967,  0.34956439,  0.35939911,
        0.36923383,  0.37906855,  0.38890327,  0.39873799,  0.40857271,
        0.41840743,  0.42824215,  0.43807687,  0.44791158,  0.4577463 ,
        0.46758102,  0.47741574,  0.48725046,  0.49708518,  0.5069199 ,
        0.51675462,  0.52658934,  0.53642406,  0.54625878,  0.5560935 ,
        0.56592822,  0.57576294,  0.58559766,  0.59543238,  0.6052671 ,
        0.61510182,  0.62493654,  0.63477126,  0.64460598,  0.6544407 ,
        0.66427542])
    y = array([1.55753875e-02,   9.20142978e-03,  -4.12695818e-02,
        -2.41837688e-02,  -5.67525066e-02,  -1.33656085e-01,
        -3.92932482e-02,  -9.93828475e-02,  -6.49584830e-02,
        -1.87308770e-02,  -8.49718973e-02,  -5.02231643e-02,
        -7.46907890e-02,  -1.06046379e-01,  -1.18749186e-01,
        -1.13565601e-01,  -4.63223010e-02,  -1.12377331e-01,
        -1.07838847e-01,  -7.19103441e-02,  -9.69987586e-02,
        -4.05682698e-02,  -5.52870296e-02,  -6.34647682e-02,
        -1.25110462e-01,  -1.33551285e-01,  -1.00283086e-01,
        -8.78261775e-02,  -1.22366831e-01,   3.73410434e-02,
        -2.16346290e-02,  -3.54011096e-02,  -5.69699146e-02,
        -1.50996730e-01,  -8.64437893e-02,  -1.06978618e-01,
        -7.91596845e-02,  -2.85348874e-02,  -4.85092625e-02,
        -7.40978718e-02,  -5.35184331e-03,  -1.41892433e-01,
        -1.09891705e-01,  -7.00225532e-02,   2.23670322e-02,
        -6.42345473e-02,  -1.14513963e-01,  -1.77867692e-02,
        -1.15476929e-01,  -7.50609040e-02,  -7.59665146e-02,
        -3.89640033e-02,  -6.71256706e-02,  -1.46708071e-01,
        -1.61539745e-02,  -5.59726842e-02,  -9.45299864e-02,
        -9.83446389e-02,  -9.92954075e-02,  -1.03172548e-01,
        -7.10593835e-02,  -5.33846729e-02,   3.38634439e-02,
         4.46895629e-01,   6.62468433e-01,   8.12499002e-02,
         3.69865060e-01,   9.81505096e-01,   5.99039435e-01,
         8.50967616e-02,  -2.67737526e-02,  -4.09759879e-02,
        -1.44490153e-01,  -2.30287295e-02,  -2.02018376e-02,
        -1.77040650e-03,  -8.88574719e-02,  -5.32815233e-02,
        -6.93681315e-02,  -7.17790946e-02,  -1.42520413e-01,
         1.50613356e-02,  -8.37516487e-02,  -9.93321687e-02,
        -2.07297225e-02,  -8.25655535e-02,   1.20367259e-02,
        -1.47362985e-02,  -4.14445363e-02,  -5.36099076e-02,
        -1.19483069e-01,  -8.75750855e-02,  -6.97698891e-02,
        -9.45113301e-02,  -5.86897917e-02,  -5.38971759e-02,
        -5.95922321e-02,   2.01958697e-02,  -5.67614287e-03,
        -4.83865663e-02,  -7.87640661e-02,  -1.30915985e-01,
        -1.45986080e-01,  -7.94370472e-02,  -5.61923422e-02,
        -9.50986519e-02,  -9.20939595e-02,  -4.11376543e-02,
        -1.29739568e-01,  -1.22105666e-01,  -8.70440751e-02,
        -8.68988112e-02,  -1.08260494e-02,  -7.90288299e-02,
        -5.83453253e-02,  -9.31360647e-02,  -6.06538579e-02,
        -3.26795094e-02,  -1.24720916e-01,   2.33035088e-02,
         4.42986703e-03,  -1.70680247e-02,  -1.65755842e-02,
         2.17673182e-01,   4.91183043e-01,   9.86441195e-01,
         1.60616803e+00,   9.89029050e-01,   1.42567885e+00,
         1.86783028e+00,   8.98995638e-01,   1.65190771e-01,
        -3.36427465e-02,  -9.43350866e-02,  -1.05553396e-01,
        -5.37899788e-03,   6.19346742e-03,  -7.22183101e-03,
        -6.04815148e-02,  -5.96636757e-02,  -6.51778141e-03,
        -8.12485069e-02,  -2.17945613e-02,  -6.93192706e-02,
        -1.69927523e-01,  -6.54176772e-02,  -6.80938214e-02,
        -1.08961679e-01,  -2.78380569e-02,  -6.92696646e-02,
        -7.72257894e-02,  -3.58553343e-02,  -8.55760425e-02,
        -5.15287071e-02,  -3.54854837e-02,  -1.05648793e-01,
        -1.01979360e-01,  -1.13662310e-01,  -5.91211058e-02,
        -4.10607755e-02,   3.95612381e-02,  -3.21216823e-04,
        -8.15489069e-02,  -7.26812184e-02,   3.38813802e-03,
         3.18101384e-02,   8.27607699e-03,  -7.05176294e-02,
        -1.20289661e-01,  -5.37291467e-02,  -4.78893109e-02,
        -8.00910443e-02,  -3.42484415e-02,  -9.23061371e-02,
        -6.11467026e-02,  -5.12490347e-02,   5.45026129e-03,
        -3.40601653e-02,  -5.34633473e-02,  -7.67978132e-02,
        -5.27321585e-02,  -7.51329362e-02,  -9.02341753e-02,
        -2.85653155e-02,  -2.81812195e-02,   1.88794062e-01,
         6.09963775e-01,   5.01888454e-01,   8.58971104e-02,
         9.57417712e-02,   5.89331269e-01,   9.81935441e-01,
         3.07364285e-01,  -3.65563519e-02,  -3.49376574e-02,
        -1.34642854e-01,  -1.64245758e-02,   6.07715966e-03,
        -3.28341946e-02,  -4.43529859e-02,  -5.27672656e-02,
         1.77064128e-02,   4.18064697e-03,   1.35755716e-02,
        -3.81845832e-02,  -4.23189811e-02,   2.52703223e-02,
        -7.12039247e-02,   5.16605303e-02,   7.01981178e-03,
        -6.71181753e-02,  -2.03371570e-02,  -1.20013859e-02,
        -1.14060365e-01,  -7.40282461e-02,  -3.78084294e-02,
        -1.20527424e-01,  -6.82442710e-02,  -1.02835357e-01,
         2.08887681e-02,  -1.96327586e-02,  -4.01197970e-02,
        -4.00166288e-02,   3.49126421e-02,  -3.74765843e-02,
        -3.12900096e-02,  -1.17622502e-02,   7.02238753e-02,
         1.81287788e-02,  -6.88833147e-02,  -8.13086852e-02,
        -8.02919865e-02,  -7.93092176e-02,  -8.38449318e-03,
        -1.22420341e-01,  -1.62812844e-02,  -9.11864787e-02,
        -1.47517873e-02,   2.48801224e-02,  -4.70457412e-02,
        -8.15037489e-02,  -6.75613731e-02,  -8.35428163e-02,
        -1.02822810e-01,  -4.38780636e-02,  -1.20214887e-01,
        -5.27682826e-02,  -1.31174894e-02,  -1.30414739e-01,
        -1.57103818e-02,  -4.95527051e-02,   2.20772102e-02,
        -1.40918205e-02,  -5.67496903e-02,   1.55445077e-02,
        -1.82226207e-02])
    parameter_guess = [.3,-0.6,width]
  plot(x,y)
  draw()  #  time.sleep(0.1)
  line = int(raw_input("Main lines (1) or all lines (-1)? "))
  set_model(line)

  result = fit_data(x,y,hpf_model,parameter_guess)
  best_fit = hpf_model(x,*result) # model y-values
  plot(x,best_fit)
  title(name)
  xlabel("Relative frequency (MHz)")
  ylabel("Amplitude (r.m.s = 1)")

  print "Amplitude: %6.3f" % (result[0]*s[0])
  print "Position:  %6.3f" % result[1]
  print "Width:     %6.3f" % result[2]
  
if __name__ == "__main__":
  test()
