import datetime
import logging
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import math

logger = logging.getLogger(__name__)

def plot(data):
    # Convert to time offset
    t0 = data['t0']

    # Get date string
    #dateStr = "%s - %s" % (data['times'][0].strftime('%Y-%m-%d-%H:%M:%S.%f')[:-4],
    #                       data['times'][-1].strftime('%Y-%m-%d-%H:%M:%S.%f')[:-4])
    dateStr = "Test"

    #tPlot = [(x - t0).microseconds/1000000.0 for x in data['times']]
    tPlot = []
    for x in data['times']:
        tPlot.append(x)
    duration = tPlot[-1] - tPlot[0]

    # Get x-limit
    xLimMag = [tPlot[0], tPlot[-1]]
    yLimMag = 2**(data['sampleSize']-1) # for +/- range
    if data['sampleSize'] > 4:
        yTicks = np.arange(-yLimMag, yLimMag+0.001, 2**(data['sampleSize']-2))
    else:
        yTicks = np.arange(-yLimMag, yLimMag+0.001, 1)

    # ----
    # Figure 1: Raw Data values
    # ----
    # Limit the number of points
    if len(data['values']) < 100000:
        plt.figure(figsize=(20,12))

        # Subplot 1: I-Signal
        ax1 = plt.subplot(211)
        plt.plot(tPlot, np.real(data['values']))
        ax1.grid(True)
        plt.ylabel('I Signal (counts)')
        plt.xlabel('Time (seconds)')
        plt.ylim([-yLimMag, yLimMag])
        plt.yticks(yTicks)
        plt.xlim(xLimMag)

        # Subplot 2: Q-Signal
        ax2 = plt.subplot(212)
        plt.plot(tPlot, np.imag(data['values']))
        ax2.grid(True)
        plt.ylabel('Q Signal (counts)')
        plt.xlabel('Time (seconds)')
        plt.ylim([-yLimMag, yLimMag])
        plt.yticks(yTicks)
        plt.xlim(xLimMag)

        # Add filename to figure?
        plt.suptitle('Raw Data\n%s\n%s' % (data['DATAFIL'], dateStr))

    # ----
    # Figure 2: Histogram
    # ----
    plt.figure(figsize=(20,12))
    if data['sampleSize'] > 8:
        nbins = 256;
        bins = np.linspace(-yLimMag-1, yLimMag+1, nbins)
    else:
        bins = np.arange(-yLimMag-2, yLimMag+3, 2) + 0.0001

    # Create Histogram 1 & normalize to percentage
    binCounts, binEdges = np.histogram(np.real(data['values']), bins)
    pctCounts = 100.0*binCounts/float(len(data['values']))
    widths = np.diff(binEdges)

    # Subplot 1: I-signal
    ax1 = plt.subplot(211)
    plt.bar(binEdges[:-1], pctCounts, widths, align='edge')
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.xlim([-yLimMag, yLimMag])
    plt.xlabel('I-Signal (counts)')
    plt.xticks(yTicks)
    plt.ylabel('Percentage Occurance')
    ax1.grid(True)

    # Create Histogram 1 & normalize to percentage
    binCounts, binEdges = np.histogram(np.imag(data['values']), bins)
    pctCounts = 100.0*binCounts/float(len(data['values']))
    widths = np.diff(binEdges)

    # Subplot 2: Q-signal
    ax2 = plt.subplot(212)
    plt.bar(binEdges[:-1], pctCounts, widths, align='edge')
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.xlim([-yLimMag, yLimMag])
    plt.xlabel('Q-Signal (counts)')
    plt.xticks(yTicks)
    plt.ylabel('Percentage Occurances')
    ax2.grid(True)

    # Get some stats
    Imean = np.mean(np.real(data['values']))
    Istd  = np.std(np.real(data['values']))
    Qmean = np.mean(np.imag(data['values']))
    Qstd  = np.std(np.imag(data['values']))


    # Add filename to figure?
    plt.suptitle('Histogram\n%s\n%s\nI-Signal: = %8.1f, sigma = %8.1f\nQ-Signal: = %8.1f, sigma = %8.1f' % (
        data['DATAFIL'], dateStr, Imean, Istd, Qmean, Qstd))

    # ----
    # Figure 3: FFT Spectrum
    # ----
    t1 = datetime.datetime.now()

    # Do FFT
    npad = int(2**math.ceil(math.log(len(data['values']), 2)) - len(data['values']))
    ytmp = np.pad(data['values'], (0, npad), mode='constant')

    t1 = datetime.datetime.now()
    fft = np.fft.fftshift(np.fft.fft(ytmp))/float(len(ytmp))

    fftPower = np.sqrt(fft*fft.conj()).real

    # Compute frequency range
    dt = duration/float(len(data['values']))
    freq = np.fft.fftshift(np.fft.fftfreq(len(fft), d=dt))
    freqLim = max(abs(freq))
    freqRes = np.mean(np.diff(freq))

    # Find peak frequency
    peakInd = np.argmax(fftPower)
    peakFreq = freq[np.argmax(fftPower)]
    xLim = [math.floor((peakFreq-100)/100)*100, math.ceil((peakFreq+100)/100)*100]


    # Reduce resolution for plotting
    Npts = 65536
    if len(freq)>(2*Npts):
        npts = np.ceil(len(freq)/float(Npts))
        thinFFT = []
        thinFreq = []
        for ind in range(0, Npts):
            ind0 = int(npts*(ind))
            ind1 = int(npts*(ind+1))

            if (ind0 >= len(freq)) or (ind0 >= len(fftPower)):
                break

            # Method 1: Max in bin
            thinFFT.append(np.max(fftPower[ind0:ind1]))
            indmax = np.argmax(fftPower[ind0:ind1])
            thinFreq.append(freq[ind0+indmax])

            # Metjhod 2: Sum in bin
            #thinFFT.append(np.sum(fftPower[ind0:ind1]))
            #thinFreq.append(np.mean(freq[ind0:ind1]))
    else:
        thinFFT = fftPower
        thinFreq = freq

    # Subplot 1: Full FFT
    plt.figure(figsize=(5,4))
    ax1 = plt.subplot(111)
    ax1.grid(True)
    plt.ylabel('Amplitude (dB)')
    plt.xlabel('Frequency (MHz)')
    plt.xlim([-freqLim/1000000.0, freqLim/1000000.0])
    if np.max(10*np.log10(fftPower)) < 25:
        plt.ylim([-40, 25])
    else:
        plt.ylim([-40, np.max(10*np.log10(fftPower))])
    plt.gca().ticklabel_format(useOffset=False, style='plain')
    #plt.plot(freq, 10*np.log10(fftPower), 'b-')
    plt.plot(np.array(thinFreq)/1000000.0, 10*np.log10(thinFFT), 'b-')

    # Subplot 2: Zoomed in FFT
    #ax2 = plt.subplot(212)
    #ax2.grid(True)
    #plt.ylabel('Amplitude (dB)')
    #plt.xlabel('Frequency (Hz)')
    #plt.xlim(xLim)
    #if np.max(10*np.log10(fftPower)) < 25:
    #    plt.ylim([-30, 25])
    #else:
    #    plt.ylim([-30, np.max(10*np.log10(fftPower))])
    #ax2.get_xaxis().get_major_formatter().set_useOffset(False)
    #ax2.get_yaxis().get_major_formatter().set_useOffset(False)
    #plt.plot(freq, 10*np.log10(fftPower), 'b-')

    # Add title
    #plt.suptitle('FFT\n%s\n%s\nPeak Freq = %1.2f Hz\nResolution = %1.2f Hz' % (
    #    data['DATAFIL'], dateStr, peakFreq, freqRes))

    # Show figure
    plt.show()

