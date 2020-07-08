"""
Module OLSR
===========

For open loop science recorder files.

Supports data files from the VSR, WVSR, PVSR, OLR, etc.

Recognizes RDEF, VDR, SFDU formats
"""
import numpy as np
import math
import sys
import datetime
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import FormatStrFormatter

# --------
# Lower Level Subfunctions
# --------
def bin2sint(bin):
    """
    binary to signed integer
    """
    if bin[0] == '1':
        # Negative
        integer = -(2**(len(bin)-1) - int(bin[1:], 2))
    else:
        integer = int(bin, 2)

    integer = int(2*integer + 1)
    return integer

def data2int(data, format):
    """
    data sample to integer
    """
    dec = 0
    if format == 'RDEF':
        for power in range(len(data)-1, -1, -1):
            dec = dec + data[power]*(2**(8*power))
    elif format == 'VDR' or format == 'SFDU':
        for power in range(0, len(data)):
            dec = dec + data[power]*2**(8*(len(data)-power-1.0))

    dec = int(dec)
    return dec

def data2float(data, format):
    """
    data sample to float
    """
    # Create binary
    binVec = '';
    if format == 'RDEF':
        for index in range(len(data)-1, -1, -1):
            binStr = "{0:b}".format(data[index])
            binStr = binStr.zfill(8)
            binVec = binVec + binStr
    elif format == 'VDR' or format == 'SFDU':
        for index in range(0, len(data), 1):
            #binStr = "{0:b}".format(data[index])
            binStr = np.binary_repr(data[index], width=8)
            binVec = binVec + binStr

    if len(data) == 8:
        # 64-bit Float
        sign = int(binVec[0], 2)
        expo = int(binVec[1:12], 2) - 2**10 + 1
        frac = 1 + int(binVec[12:], 2)/float(2**52)
        value = ((-1)**sign) * (2**expo) * frac
    elif len(data) == 4:
        # 32-bit Float
        sign = int(binVec[0], 2)
        expo = int(binVec[1:9], 2) - 2**7 + 1
        frac = 1 + int(binVec[9:], 2)/float(2**24)
        value = ((-1)**sign) * (2**expo) * frac
    else:
        raise NameError('Floating point number must contain 32 or 64 bits!')

    return value

def to_percent(y, position):
    """
    Convert a float to percent
    """
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] is True:
        return s + r'$\%$'
    else:
        return s + '%'

def readHeaders(f, options):
    """
    decode the file header
    """
    if options['format'] == 'RDEF':
        # Header Elements
        headers = ['RECORD_LABEL','RECORD_LENGTH', 'RECORD_VERSION_ID',
                'STATION_ID', 'SPACECRAFT_ID', 'SAMPLE_SIZE', 'SAMPLE_RATE',
                'VALIDITY_FLAG', 'AGENCY_FLAG', 'RF_TO_IF_DOWNCONV', 'IF_TO_CHANNEL_DOWNCONV',
                'TIME_TAG_YEAR', 'TIME_TAG_DOY', 'TIME_TAG_SECOND_OF_DAY',
                'TIMETAG_PICOSECONDS_OF_THE_SECOND', 'CHANNEL_ACCUMULATED_PHASE',
                'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT0', 'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT1',
                'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT2', 'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT3',
                'EMPTY_FIELDS1', 'PREDICT_PASS_NUMBER', 'UPLINK_BAND', 'DOWNLINK_BAND',
                'TRACK_MODE', 'UPLINK_DSS_ID', 'OLR_ID', 'OLR_SOFTWARE_VERSION', 'CHANNEL_POWER_CALIBRATION_FACTOR',
                'TOTAL_FREQUENCY_OFFSET', 'CHANNEL_NUMBER', 'EMPTY_FIELDS2', 'END_LABEL'];

        dict = {'RECORD_LABEL':{'Size':4, 'Type':'CHARACTER'},
                'RECORD_LENGTH':{'Size':4, 'Type':'UNSIGNED INTEGER'},
                'RECORD_VERSION_ID':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'STATION_ID':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'SPACECRAFT_ID':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'SAMPLE_SIZE':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'SAMPLE_RATE':{'Size':4, 'Type':'UNSIGNED INTEGER'},
                'VALIDITY_FLAG':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'AGENCY_FLAG':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'RF_TO_IF_DOWNCONV':{'Size':8, 'Type':'FLOATING POINT'},
                'IF_TO_CHANNEL_DOWNCONV':{'Size':8, 'Type':'FLOATING POINT'},
                'TIME_TAG_YEAR':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'TIME_TAG_DOY':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'TIME_TAG_SECOND_OF_DAY':{'Size':4, 'Type':'UNSIGNED INTEGER'},
                'TIMETAG_PICOSECONDS_OF_THE_SECOND':{'Size':8, 'Type':'FLOATING POINT'},
                'CHANNEL_ACCUMULATED_PHASE':{'Size':8, 'Type':'FLOATING POINT'},
                'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT0':{'Size':8, 'Type':'FLOATING POINT'},
                'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT1':{'Size':8, 'Type':'FLOATING POINT'},
                'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT2':{'Size':8, 'Type':'FLOATING POINT'},
                'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT3':{'Size':8, 'Type':'FLOATING POINT'},
                'EMPTY_FIELDS1':{'Size':36, 'Type':'EMPTY'},
                'PREDICT_PASS_NUMBER':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'UPLINK_BAND':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'DOWNLINK_BAND':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'TRACK_MODE':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'UPLINK_DSS_ID':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'OLR_ID':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'OLR_SOFTWARE_VERSION':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'CHANNEL_POWER_CALIBRATION_FACTOR':{'Size':4, 'Type':'FLOATING POINT'},
                'TOTAL_FREQUENCY_OFFSET':{'Size':8, 'Type':'FLOATING POINT'},
                'CHANNEL_NUMBER':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'EMPTY_FIELDS2':{'Size':19, 'Type':'EMPTY'},
                'END_LABEL':{'Size':4, 'Type':'UNSIGNED INTEGER'}}

        # Define number of header bytes
        headerSize = 176

    elif options['format'] == 'VDR':
        # Header Elements
        headers = ['RECORD_LABEL', 'RECORD_LENGTH', 'RECORD_VERSION',
                'WVSR_SOFTWARE_VERSION', 'SPC_ID', 'WVSR_ID', 'CHAN_ID',
                'SAMPLE_SIZE', 'SAMPLE_RATE', 'CHAN_DOPPLER_MODE',
                'PRDX_DSS_ID', 'PRDX_SC_ID', 'PRDX_PASS_NUMBER', 'PRDX_UL_BAND',
                'PRDX_DL_BAND', 'TRK_MODE', 'UL_DSS_ID', 'DDCLO_MHZ', 'RF_TO_IF_DOWNCONV_MHZ',
                'DATA_ERROR', 'TIME_TAG_YEAR', 'TIME_TAG_DOY', 'TIME_TAG_SECOND_OF_DAY',
                'DATA_TIME_OFFSET_NANOSECONDS', 'PREDICTS_FREQ_OVERRIDE',
                'PREDICTS_FREQ_OFFSET', 'PREDICTS_FREQ_RATE', 'CHANNEL_FREQ_OFFSET',
                'RF_FREQ_POINT', 'CHANNEL_ACCUMULATED_PHASE',
                'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT1',
                'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT2',
                'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT3',
                'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT4', 'CHAN_LABEL']

        dict = {'RECORD_LABEL':{'Size':4, 'Type':'CHARACTER'},
            'RECORD_LENGTH':{'Size':4, 'Type':'UNSIGNED INTEGER'},
            'RECORD_VERSION':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'WVSR_SOFTWARE_VERSION':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'SPC_ID':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'WVSR_ID':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'CHAN_ID':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'SAMPLE_SIZE':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'SAMPLE_RATE':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'CHAN_DOPPLER_MODE':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'PRDX_DSS_ID':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'PRDX_SC_ID':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'PRDX_PASS_NUMBER':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'PRDX_UL_BAND':{'Size':2, 'Type':'CHARACTER'},
            'PRDX_DL_BAND':{'Size':2, 'Type':'CHARACTER'},
            'TRK_MODE':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'UL_DSS_ID':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'DDCLO_MHZ':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'RF_TO_IF_DOWNCONV_MHZ':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'DATA_ERROR':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'TIME_TAG_YEAR':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'TIME_TAG_DOY':{'Size':2, 'Type':'UNSIGNED INTEGER'},
            'TIME_TAG_SECOND_OF_DAY':{'Size':4, 'Type':'UNSIGNED INTEGER'},
            'DATA_TIME_OFFSET_NANOSECONDS':{'Size':4, 'Type':'UNSIGNED INTEGER'},
            'PREDICTS_FREQ_OVERRIDE':{'Size':8, 'Type':'FLOATING POINT'},
            'PREDICTS_FREQ_OFFSET':{'Size':8, 'Type':'FLOATING POINT'},
            'PREDICTS_FREQ_RATE':{'Size':8, 'Type':'FLOATING POINT'},
            'CHANNEL_FREQ_OFFSET':{'Size':8, 'Type':'FLOATING POINT'},
            'RF_FREQ_POINT':{'Size':8, 'Type':'FLOATING POINT'},
            'CHANNEL_ACCUMULATED_PHASE':{'Size':8, 'Type':'FLOATING POINT'},
            'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT1':{'Size':8, 'Type':'FLOATING POINT'},
            'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT2':{'Size':8, 'Type':'FLOATING POINT'},
            'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT3':{'Size':8, 'Type':'FLOATING POINT'},
            'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT4':{'Size':8, 'Type':'FLOATING POINT'},
            'CHAN_LABEL':{'Size':16, 'Type':'CHARACTER'}}

        # Define number of header bytes
        headerSize = 152

    elif options['format'] == 'SFDU':
        # Header Elements
        headers = ['CONTROL_AUTHORITY_ID', 'VERSION_ID', 'CLASS_ID', 'RESERVED1',
                'DATA_DESCRIPTION_ID', 'LENGTH', 'HACHDO_TYPE', 'HACHDO_LENGTH',
                'PHCHDO_TYPE', 'PHCHDO_LENGTH', 'DATA_CLASS', 'DATA_SOURCE',
                'MISSION_ID', 'FORMAT_CODE', 'SHCHDO_TYPE', 'SHCHDO_LENGTH',
                'ORIGINATOR_ID', 'LAST_MODIFIED_ID', 'RSR_SOFTWARE_ID',
                'RECORD_SEQUENCE_NUMBER', 'SPC_ID', 'DSS_ID', 'RSR_ID', 'SCCHAN_ID',
                'RESERVED2', 'SPACECRAFT_ID', 'PRDX_PASS_NUMBER', 'UL_BAND', 'DL_BAND',
                'TRK_MODE', 'UL_DSS_ID', 'FGAIN_PX_NO', 'FGAIN_IF_BANDWIDTH', 'FROV_FLAG',
                'ATTENUATION', 'ADC_RMS', 'ADC_PEAK', 'ADC_TIME_TAG_YEAR', 'ADC_TIME_TAG_DOY',
                'ADC_TIME_TAG_SECOND_OF_DAY', 'SAMPLE_SIZE', 'DATA_ERROR', 'SAMPLE_RATE',
                'DDC_LO', 'RF_IF_LO', 'SFDU_TIME_TAG_YEAR', 'SFDU_TIME_TAG_DOY',
                'SFDU_TIME_TAG_SECOND_OF_DAY', 'PREDICTS_TIME_SHIFT', 'PREDICTS_FREQ_OVERRIDE',
                'PREDICTS_FREQ_RATE', 'PREDICTS_FREQ_OFFSET', 'SUBCHANNEL_FREQ_OFFSET',
                'RF_FREQ_POINT_1', 'RF_FREQ_POINT_2', 'RF_FREQ_POINT_3', 'SCHAN_FREQ_POINT_1',
                'SCHAN_FREQ_POINT_2', 'SCHAN_FREQ_POINT_3', 'SCHAN_FREQ_POLY_COEF_1',
                'SCHAN_FREQ_POLY_COEF_2', 'SCHAN_FREQ_POLY_COEF_3', 'SCHAN_ACCUM_PHASE',
                'SCHAN_PHASE_POLY_COEF_1', 'SCHAN_PHASE_POLY_COEF_2', 'SCHAN_PHASE_POLY_COEF_3',
                'SCHAN_PHASE_POLY_COEF_4', 'SCHAN_MULT', 'RESERVED3']

        dict = {'CONTROL_AUTHORITY_ID':{'Size':4, 'Type':'CHARACTER'},
                'VERSION_ID':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'CLASS_ID':{'Size':1, 'Type':'CHARACTER'},
                'RESERVED1':{'Size':2, 'Type':'EMPTY'},
                'DATA_DESCRIPTION_ID':{'Size':4, 'Type':'CHARACTER'},
                'LENGTH':{'Size':8, 'Type':'UNSIGNED INTEGER'},
                'HACHDO_TYPE':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'HACHDO_LENGTH':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'PHCHDO_TYPE':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'PHCHDO_LENGTH':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'DATA_CLASS':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'DATA_SOURCE':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'MISSION_ID':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'FORMAT_CODE':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'SHCHDO_TYPE':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'SHCHDO_LENGTH':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'ORIGINATOR_ID':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'LAST_MODIFIED_ID':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'RSR_SOFTWARE_ID':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'RECORD_SEQUENCE_NUMBER':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'SPC_ID':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'DSS_ID':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'RSR_ID':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'SCCHAN_ID':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'RESERVED2':{'Size':1, 'Type':'EMPTY'},
                'SPACECRAFT_ID':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'PRDX_PASS_NUMBER':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'UL_BAND':{'Size':1, 'Type':'CHARACTER'},
                'DL_BAND':{'Size':1, 'Type':'CHARACTER'},
                'TRK_MODE':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'UL_DSS_ID':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'FGAIN_PX_NO':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'FGAIN_IF_BANDWIDTH':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'FROV_FLAG':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'ATTENUATION':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'ADC_RMS':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'ADC_PEAK':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'ADC_TIME_TAG_YEAR':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'ADC_TIME_TAG_DOY':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'ADC_TIME_TAG_SECOND_OF_DAY':{'Size':4, 'Type':'UNSIGNED INTEGER'},
                'SAMPLE_SIZE':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'DATA_ERROR':{'Size':1, 'Type':'UNSIGNED INTEGER'},
                'SAMPLE_RATE':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'DDC_LO':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'RF_IF_LO':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'SFDU_TIME_TAG_YEAR':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'SFDU_TIME_TAG_DOY':{'Size':2, 'Type':'UNSIGNED INTEGER'},
                'SFDU_TIME_TAG_SECOND_OF_DAY':{'Size':8, 'Type':'FLOATING POINT'},
                'PREDICTS_TIME_SHIFT':{'Size':8, 'Type':'FLOATING POINT'},
                'PREDICTS_FREQ_OVERRIDE':{'Size':8, 'Type':'FLOATING POINT'},
                'PREDICTS_FREQ_RATE':{'Size':8, 'Type':'FLOATING POINT'},
                'PREDICTS_FREQ_OFFSET':{'Size':8, 'Type':'FLOATING POINT'},
                'SUBCHANNEL_FREQ_OFFSET':{'Size':8, 'Type':'FLOATING POINT'},
                'RF_FREQ_POINT_1':{'Size':8, 'Type':'FLOATING POINT'},
                'RF_FREQ_POINT_2':{'Size':8, 'Type':'FLOATING POINT'},
                'RF_FREQ_POINT_3':{'Size':8, 'Type':'FLOATING POINT'},
                'SCHAN_FREQ_POINT_1':{'Size':8, 'Type':'FLOATING POINT'},
                'SCHAN_FREQ_POINT_2':{'Size':8, 'Type':'FLOATING POINT'},
                'SCHAN_FREQ_POINT_3':{'Size':8, 'Type':'FLOATING POINT'},
                'SCHAN_FREQ_POLY_COEF_1':{'Size':8, 'Type':'FLOATING POINT'},
                'SCHAN_FREQ_POLY_COEF_2':{'Size':8, 'Type':'FLOATING POINT'},
                'SCHAN_FREQ_POLY_COEF_3':{'Size':8, 'Type':'FLOATING POINT'},
                'SCHAN_ACCUM_PHASE':{'Size':8, 'Type':'FLOATING POINT'},
                'SCHAN_PHASE_POLY_COEF_1':{'Size':8, 'Type':'FLOATING POINT'},
                'SCHAN_PHASE_POLY_COEF_2':{'Size':8, 'Type':'FLOATING POINT'},
                'SCHAN_PHASE_POLY_COEF_3':{'Size':8, 'Type':'FLOATING POINT'},
                'SCHAN_PHASE_POLY_COEF_4':{'Size':8, 'Type':'FLOATING POINT'},
                'SCHAN_MULT':{'Size':4, 'Type':'FLOATING POINT'},
                'RESERVED3':{'Size':12, 'Type':'EMPTY'}
                };

        headerSize = 256

    # Read header
    headerBytes = np.fromfile(f, dtype=np.uint8, count=headerSize)

    # Initialize output
    output = {'bytes':headerBytes}

    # Check for errors
    if len(headerBytes) != headerSize:
        output = {'EOF':True}
        return output
    
    startByte = 0
    for header in headers:
        # Extract data bytes
        endByte = startByte + dict[header]['Size']
        dataBytes = headerBytes[startByte:endByte]

        # Convert
        if dict[header]['Type'] == 'CHARACTER':
            data = "".join(map(chr, dataBytes))
        elif dict[header]['Type'] == 'UNSIGNED INTEGER':
            data = data2int(dataBytes, options['format'])
        elif dict[header]['Type'] == 'FLOATING POINT':
            data = data2float(dataBytes, options['format'])

        # Save data
        if dict[header]['Type'] != 'EMPTY':
            output[header] = data

        # Increment start byte
        startByte = endByte

    # Print Headers
    if options['toPrintHeaders']:
        print('\n')
        for header in headers:
            if header in output:
                print(header.ljust(40) + ': ' + str(output[header]))
        print('\n')

    return output

def checkFormat(DATAFIL):
    """
    check the file format
    """
    # Open file, read first 4 bytes as string
    f = open(DATAFIL, 'rb')
    dataBytes = np.fromfile(f, dtype=np.uint8, count=4)
    data = "".join(map(chr, dataBytes))
    f.close()

    # If string matches a particular format, return format
    if data == 'RDEF':
        output = 'RDEF'
    elif data == 'VSRD':
        output = 'VDR'
    elif data == 'NJPL':
        output = 'SFDU'
    else:
        raise NameError('Data file is not a recognized format')
        return

    return output

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
    yLimMag = 2**(data['sampleSize'])
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

# --------
# Main File I/O Subfunctions
# --------

# ----
# RDEF
#   This script reads in RDEF data formatted files
# ----
def RDEF(DATAFIL, options):
    # Default Inputs
    # DATAFIL = '/Users/cvolk/Documents/MATLAB/RSR_Test/olr4.gdscc_stab_test_c1.18_088_035414'
    # DATAFIL = ''
    if 'startTime' not in options:
        timeRange = [0, float('inf')]
    else:
        timeRange = [options['startTime'], float('inf')]

    if 'startDate' not in options:
        options['startDate'] = [datetime.datetime.strptime('1970-01-01', '%Y-%m-%d')]

    if 'file' not in options:
        # Open file
        f = open(DATAFIL, 'rb')
        curRecord = 1
    else:
        f = options['file']['f']
        curRecord = options['file']['curRecord']

    if 'leaveFileOpen' not in options:
        leaveFileOpen = False
    else:
        leaveFileOpen = True

    # Read headers
    tmpOptions = {'toPrintHeaders':False, 'format':options['format']}
    headers = readHeaders(f, tmpOptions)
    if 'EOF' in headers:
        raise NameError('Data file does not contain valid headers!')
        return

    # Alter time range if looking for date
    if options['lookForDate']:
        tStr = '%04d-%03d' % (headers['TIME_TAG_YEAR'], headers['TIME_TAG_DOY'])
        t0 = datetime.datetime.strptime(tStr, '%Y-%j')
        t0 = (t0 + datetime.timedelta(0,headers['TIME_TAG_SECOND_OF_DAY'])
                + datetime.timedelta(0,headers['TIMETAG_PICOSECONDS_OF_THE_SECOND'])/1000000000000)

        td = (options['startDate'] - t0)
        deltaT = ((td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) // 10.0**6)
        if deltaT > 0:
            timeRange[0] = deltaT
    else:
        tStr = '%04d-%03d' % (headers['TIME_TAG_YEAR'], headers['TIME_TAG_DOY'])
        t0 = (datetime.datetime.strptime(tStr, '%Y-%j') + datetime.timedelta(0,headers['TIME_TAG_SECOND_OF_DAY'])
                                   + datetime.timedelta(0,headers['TIMETAG_PICOSECONDS_OF_THE_SECOND'])/1000000000000)
        options['startDate'] = (t0 + datetime.timedelta(0,timeRange[0]))

    # Apply duration (if applicable)
    if 'duration' in options:
        timeRange[1] = timeRange[0] + options['duration']
    options['endDate'] = options['startDate'] + datetime.timedelta(0, timeRange[1]-timeRange[0])

    # Calculate data rate, time step
    bytesPerSecond = int(2*headers['SAMPLE_SIZE']*headers['SAMPLE_RATE']/8.0)
    timeStep = 1/float(headers['SAMPLE_RATE'])

    # Rewind file
    f.seek(-176, 1)

    # Determine start and stop records
    startRecord = int(math.floor(timeRange[0]) + 1)
    stopRecord = int(math.floor(timeRange[1]) + 1)

    # Determine start and stop bytes (to nearest word)
    startByte = int(math.floor(bytesPerSecond*(timeRange[0]%1)//4.0)*4.0)
    stopByte  = int(math.ceil(bytesPerSecond*(timeRange[1]%1)//4.0)*4.0)

    if stopByte == 0:
        stopRecord = stopRecord - 1
        stopByte = int(bytesPerSecond)

    # Start reading data
    output = {'values':[], 'times':[], 'phase':[], 'phCoeff':[],
                  'LO': headers['RF_TO_IF_DOWNCONV'],
                  'DDCLO':headers['IF_TO_CHANNEL_DOWNCONV'],
                  'sampleSize':headers['SAMPLE_SIZE'],
                  'DATAFIL':DATAFIL,
                  't0':t0,
                  'file': {'f':[], 'curRecord':[]}}

    if headers['SAMPLE_SIZE'] == 16:
        dt = np.dtype('uint8')
    else:
        dt = np.dtype('uint8')

    while curRecord <= stopRecord:
        # Clear memory
        dataBytes = []
        if curRecord < startRecord:
            # Skip this section
            f.seek(176, 1)
            f.seek(bytesPerSecond, 1)

        elif (curRecord == startRecord) and (curRecord == stopRecord):
            # This is the only section we need to read
            headers = readHeaders(f, options)
            f.seek(startByte, 1)
            dataBytes = np.fromfile(f, dtype=dt, count=int(stopByte-startByte))

        elif curRecord == startRecord:
            # This is the first section to read
            headers = readHeaders(f, options)
            f.seek(startByte, 1)
            dataBytes = np.fromfile(f, dtype=dt, count=int(bytesPerSecond-startByte))

        elif (curRecord == stopRecord) and (stopByte > 0):
            # This is the last section we need to read
            headers = readHeaders(f, options)
            dataBytes = np.fromfile(f, dtype=dt, count=int(stopByte))

        elif (curRecord > startRecord) and (curRecord < stopRecord):
            # This is one of the sections we need to read
            headers = readHeaders(f, options)
            dataBytes = np.fromfile(f, dtype=dt, count=int(bytesPerSecond))

        # Check for EOF
        if 'EOF' in headers:
            break

        if len(dataBytes) > 0:
            # Get samples per word, etc
            nSamples = 16/headers['SAMPLE_SIZE']


            if headers['SAMPLE_SIZE'] == 2:
                # Split into 2-bit samples
                tmp = np.empty( (len(dataBytes),4) )
                for i in range(0,4):
                    tmp[:,i] = dataBytes/(2**(2*i)) % 2**2
                    dataBytes = dataBytes - tmp[:,i]*(2**(2*i))

                # Change datatype before 2's complement
                data = np.array(tmp).astype(np.int8)

            elif headers['SAMPLE_SIZE'] == 4:
                # Split into 4-bit samples
                tmp = np.empty( (len(dataBytes),2) )
                for i in range(0,2):
                    tmp[:,i] = dataBytes/(4**(2*i)) % 2**4
                    dataBytes = dataBytes - tmp[:,i]*(4**(2*i))

                # Change datatype before 2's complement
                data = np.array(tmp).astype(np.int8)

            elif headers['SAMPLE_SIZE'] == 8:
                # Change datatype before 2's complement
                data = np.array(dataBytes).astype(np.int16)

            elif headers['SAMPLE_SIZE'] == 16:
                # Reshape into 16-bit rows
                dataBytes = np.reshape(dataBytes, (-1, 2))

                # Combine row into one 16-bit unsigned integer
                #   Convert to 32-bit signed integer for conversion to engineering values
                data = np.array(dataBytes[:,0] + (2**8)*dataBytes[:,1]).astype(np.int32)

            else:
                print('')
                print('ERROR: SAMPLE SIZE %d is not currently supported!!' % headers['SAMPLE_SIZE'])
                print('')
                stopHere

            # Reshape into samples (real, imag)
            data = np.reshape(data, (-1, 2))

            # Perform 2's complement
            b = (data >= 2**(headers['SAMPLE_SIZE']-1))
            data[b] = data[b] - 2**headers['SAMPLE_SIZE']

            # Convert to complex numbers
            output['values'] = np.append(output['values'],
                                         (2*data[:,0]+1) + 1j*(2*data[:,1]+1))
            N = len(data[:,0])

            # Create time vector
            if curRecord == startRecord:
                tVec = startByte/float(bytesPerSecond) + timeStep*np.arange(0, N)
            else:
                tVec = timeStep*np.arange(0, N)
            output['times'] = np.append(output['times'], (curRecord-1) + tVec)

            # Get phase coefficients
            c0 = (headers['CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT0']
                + headers['CHANNEL_ACCUMULATED_PHASE'])
            c1 = headers['CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT1']
            c2 = headers['CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT2']
            c3 = headers['CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT3']

            output['phCoeff'].append([c0, c1, c2, c3])
            output['phase'] = np.append(output['phase'], c0 + c1*tVec + c2*(tVec**2) + c3*(tVec**3))

            if options['toPrintData']:
                printStr = ''
                for val in output['values']:
                    if headers['SAMPLE_SIZE'] <= 2:
                        printStr = printStr + ( '%3s%2d %2d' % ('', val.real, val.imag) )
                    elif headers['SAMPLE_SIZE'] == 4:
                        printStr = printStr + ( '%4s%3d %3d' % ('', val.real, val.imag) )
                    elif headers['SAMPLE_SIZE'] == 8:
                        printStr = printStr + ( '%5s%4d %4d' % ('', val.real, val.imag) )
                    elif headers['SAMPLE_SIZE'] == 16:
                        printStr = printStr + ( '%6s%6d %6d' % ('', val.real, val.imag) )

                    if len(printStr) > 108:
                        print(printStr)
                        printStr = ''

        # Increment record
        if (curRecord == stopRecord) and leaveFileOpen:
            # Rewind to start of record
            f.seek(-stopByte, 1)
            f.seek(-176, 1)
            break

        curRecord = curRecord + 1

    # Close file
    if not leaveFileOpen:
        f.close()
    else:
        output['file']['curRecord'] = curRecord
        output['file']['f'] = f


    # Return outputs
    return output

# ----
#     VDR
# ----
def VDR(DATAFIL, options):
    # Default Inputs
    # DATAFIL = '/Users/cvolk/Documents/MATLAB/RSR_Test/olr4.gdscc_stab_test_c1.18_088_035414'
    # DATAFIL = ''
    if 'startTime' not in options:
        timeRange = [0, float('inf')]
    else:
        timeRange = [options['startTime'], float('inf')]

    if 'startDate' not in options:
        options['startDate'] = [datetime.datetime.strptime('1970-01-01', '%Y-%m-%d')]

    if 'file' not in options:
        # Open file
        f = open(DATAFIL, 'rb')
        curRecord = 1
    else:
        f = options['file']['f']
        curRecord = options['file']['curRecord']

    if 'leaveFileOpen' not in options:
        leaveFileOpen = False
    else:
        leaveFileOpen = True

    # Read headers
    tmpOptions = {'toPrintHeaders':False, 'format':options['format']}
    headers = readHeaders(f, tmpOptions)
    if 'EOF' in headers:
        raise NameError('Data file does not contain valid headers!')
        return

    # Determine y-limit
    yLimMag = 2**headers['SAMPLE_SIZE']

    # Alter time range if looking for date
    if options['lookForDate']:
        tStr = '%04d-%03d' % (headers['TIME_TAG_YEAR'], headers['TIME_TAG_DOY'])
        t0 = datetime.datetime.strptime(tStr, '%Y-%j')
        t0 = (t0 + datetime.timedelta(0,headers['TIME_TAG_SECOND_OF_DAY'])
                + datetime.timedelta(0,headers['DATA_TIME_OFFSET_NANOSECONDS'])/1000000000)

        td     = (options['startDate'] - t0)
        deltaT = ((td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) // 10.0**6)
        if deltaT > 0:
            timeRange[0] = deltaT
    else:
        tStr = '%04d-%03d' % (headers['TIME_TAG_YEAR'], headers['TIME_TAG_DOY'])
        t0 = (datetime.datetime.strptime(tStr, '%Y-%j') + datetime.timedelta(0,headers['TIME_TAG_SECOND_OF_DAY'])
                                   + datetime.timedelta(0,headers['DATA_TIME_OFFSET_NANOSECONDS'])/1000000000)
        options['startDate'] = (t0 + datetime.timedelta(0,timeRange[0]))

    # Apply duration (if applicable)
    if 'duration' in options:
        timeRange[1] = timeRange[0] + options['duration']

    options['endDate'] = options['startDate'] + datetime.timedelta(0, timeRange[1]-timeRange[0])

    # Calculate data rate, time step
    bytesPerSecond = 2000*headers['SAMPLE_SIZE']*headers['SAMPLE_RATE']/8.0
    timeStep = 1/(1000.0*float(headers['SAMPLE_RATE']))

    # Rewind file
    f.seek(0, 0)

    # Determine start and stop records
    startRecord = math.floor(timeRange[0]) + 1
    stopRecord = math.floor(timeRange[1]) + 1
    curRecord = 1

    # Determine start and stop bytes (to nearest word)
    startByte = math.floor(bytesPerSecond*(timeRange[0]%1)//4.0)*4.0
    stopByte  = math.ceil(bytesPerSecond*(timeRange[1]%1)//4.0)*4.0

    # Start reading data
    output = {'values':[], 'times':[], 'phase':[], 'phCoeff':[],
                  'LO':1000000*headers['RF_TO_IF_DOWNCONV_MHZ'],
                  'DDCLO':1000000*headers['DDCLO_MHZ'],
                  'sampleSize':headers['SAMPLE_SIZE'],
                  'DATAFIL':DATAFIL,
                  't0':t0,
                  'file': {'f':[], 'curRecord':[]},
                  'headers':[]}

    while curRecord <= stopRecord:
        # Clear memory
        dataBytes = []

        if curRecord < startRecord:
            # Skip this section
            f.seek(152, 1)
            f.seek(bytesPerSecond, 1)

        elif (curRecord == startRecord) and (curRecord == stopRecord):
            # This is the only section we need to read
            headers = readHeaders(f, options)
            f.seek(startByte, 1)
            dataBytes = np.fromfile(f, dtype=np.uint8, count=int(stopByte-startByte))

        elif curRecord == startRecord:
            # This is the first section to read
            headers = readHeaders(f, options)
            f.seek(startByte, 1)
            dataBytes = np.fromfile(f, dtype=np.uint8, count=int(bytesPerSecond-startByte))

        elif (curRecord == stopRecord) and (stopByte > 0):
            # This is the last section we need to read
            headers = readHeaders(f, options)
            dataBytes = np.fromfile(f, dtype=np.uint8, count=int(stopByte))

        elif (curRecord > startRecord) and (curRecord < stopRecord):
            # This is one of the sections we need to read
            headers = readHeaders(f, options)
            dataBytes = np.fromfile(f, dtype=np.uint8, count=int(bytesPerSecond))

        # Check for EOF
        if 'EOF' in headers:
            break

        if len(dataBytes) > 0:
            # Convert 8-bit to sample size... worried about memory size here...
            if headers['SAMPLE_SIZE'] == 2:
                # Convert to binary
                data = np.unpackbits(dataBytes)

                # Convert binary to 2-bit
                data = 2*data[::2] + data[1::2]

                # Get ready to 2's complement and 2*n + 1
                data = np.array(data).astype(np.int8)

            elif headers['SAMPLE_SIZE'] == 8:
                # Get ready to 2's complement and 2*n + 1
                data = np.array(dataBytes).astype(np.int16)

            elif headers['SAMPLE_SIZE'] == 16:
                # Reshape into 16-bit rows
                dataBytes = np.reshape(dataBytes, (-1, 2))

                # Combine row into one 16-bit unsigned integer
                #   Convert to 32-bit signed integer for conversion to engineering values
                data = np.array(dataBytes[:,0] + (2**8)*dataBytes[:,1]).astype(np.int32)

            else:
                print('')
                print('ERROR: SAMPLE SIZE %d is not currently supported!!' % headers['SAMPLE_SIZE'])
                print('')
                stopHere

            # Add headers to output
            output['headers'].append(headers)

            # Perform 2's complement
            b = (data >= 2**(headers['SAMPLE_SIZE']-1))
            data[b] = data[b] - 2**headers['SAMPLE_SIZE']

            # Reshape into words, and fliplr (I thought?)
            Nsamps = 32 / headers['SAMPLE_SIZE']
            data = np.reshape(data, (-1, Nsamps))
            data = np.fliplr(data)

            # Reshape into complex numbers
            data = np.reshape((2*data[:, 0:(Nsamps/2)]+1) + 1j*(2*data[:, (Nsamps/2):]+1), (-1, 1))

            # Convert to complex numbers
            output['values'] = np.append( output['values'], data )
            N = len(data[:,0])

            # Create time vector
            if curRecord == startRecord:
                tVec = startByte/float(bytesPerSecond) + timeStep*np.arange(0, N)
            else:
                tVec = timeStep*np.arange(0, N)
            output['times'] = np.append(output['times'], (curRecord-1) + tVec)

            # Get phase coefficients
            c0 = (headers['CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT1']
                + headers['CHANNEL_ACCUMULATED_PHASE'])
            c1 = headers['CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT2']
            c2 = headers['CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT3']
            c3 = headers['CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT4']

            output['phCoeff'].append([c0, c1, c2, c3])
            output['phase'] = np.append(output['phase'], c0 + c1*tVec + c2*(tVec**2) + c3*(tVec**3))

            if options['toPrintData']:
                printStr = ''
                for val in data:
                    if headers['SAMPLE_SIZE'] <= 2:
                        printStr = printStr + ( '%3s%2d %2d' % ('', val.real, val.imag) )
                    elif headers['SAMPLE_SIZE'] == 4:
                        printStr = printStr + ( '%4s%3d %3d' % ('', val.real, val.imag) )
                    elif headers['SAMPLE_SIZE'] == 8:
                        printStr = printStr + ( '%5s%4d %4d' % ('', val.real, val.imag) )
                    elif headers['SAMPLE_SIZE'] == 16:
                        printStr = printStr + ( '%6s%6d %6d' % ('', val.real, val.imag) )

                    if len(printStr) > 108:
                        print(printStr)
                        printStr = ''

        # Increment record
        if (curRecord == stopRecord) and (stopByte > 0) and leaveFileOpen:
            # Rewind to start of record
            f.seek(-stopByte, 1)
            f.seek(-152, 1)
            break

        curRecord = curRecord + 1

    # Close file
    if not leaveFileOpen:
        f.close()
    else:
        output['file']['curRecord'] = curRecord
        output['file']['f'] = f

    # Return outputs
    return output

# ----
# RSR Format
# ----

def SFDU(DATAFIL, options):
    # Open file
    f = open(DATAFIL, 'rb')

    # Look for start date
    if 'startTime' not in options:
        timeRange = [0, float('inf')]
    else:
        timeRange = [options['startTime'], float('inf')]

    if 'startDate' not in options:
        options['startDate'] = [datetime.datetime.strptime('1970-01-01', '%Y-%m-%d')]

    # Open file
    f = open(DATAFIL, 'rb')

    # Read headers
    tmpOptions = {'toPrintHeaders':False, 'format':options['format']}
    headers = readHeaders(f, tmpOptions)
    if 'EOF' in headers:
        raise NameError('Data file does not contain valid headers!')
        return

    # Determine y-limit
    yLimMag = 2**headers['SAMPLE_SIZE']

    # Alter time range if looking for date
    if options['lookForDate']:
        tStr = '%04d-%03d' % (headers['SFDU_TIME_TAG_YEAR'], headers['SFDU_TIME_TAG_DOY'])
        t0 = datetime.datetime.strptime(tStr, '%Y-%j')
        t0 = (t0 + datetime.timedelta(0,headers['SFDU_TIME_TAG_SECOND_OF_DAY']))

        td     = (options['startDate'] - t0)
        deltaT = ((td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) // 10.0**6)
        if deltaT > 0:
            timeRange[0] = deltaT
    else:
        tStr = '%04d-%03d' % (headers['SFDU_TIME_TAG_YEAR'], headers['SFDU_TIME_TAG_DOY'])
        t0 = datetime.datetime.strptime(tStr, '%Y-%j') + datetime.timedelta(0,headers['SFDU_TIME_TAG_SECOND_OF_DAY'])
        options['startDate'] = (t0 + datetime.timedelta(0,timeRange[0]))

    # Apply duration (if applicable)
    if 'duration' in options:
        timeRange[1] = timeRange[0] + options['duration']
    options['endDate'] = options['startDate'] + datetime.timedelta(0, timeRange[1]-timeRange[0])


    # Calculate data rate, time step
    bytesPerSecond = 2000*headers['SAMPLE_SIZE']*headers['SAMPLE_RATE']/8.0
    timeStep = 1/float(headers['SAMPLE_RATE']*1000)

    # Rewind file
    f.seek(0, 0)

    # Determine start and stop records
    startRecord = math.floor(timeRange[0]) + 1
    stopRecord = math.floor(timeRange[1]) + 1
    curRecord = 1

    # Determine start and stop bytes (to nearest word)
    #startByte = math.floor(bytesPerSecond*(timeRange[0]%1)/2.0)*2.0
    #stopByte  = math.ceil(bytesPerSecond*(timeRange[1]%1)/2.0)*2.0
    # Some floating point issues cause this to not work always... see floor(8.2-0.2)=7
    startByte = 0
    stopByte = bytesPerSecond

    # Initialize output dictionary
    output = {'values':[], 'times':[], 'phase':[],
              'LO':1000000*headers['RF_IF_LO'],
              'DDCLO':1000000*headers['DDC_LO'],
              'sampleSize':headers['SAMPLE_SIZE'],
              'DATAFIL':DATAFIL,
              't0':t0}

    # Start reading data
    while curRecord <= stopRecord:
        # Clear memory
        dataBytes = []

        if curRecord < startRecord:
            # Skip this section
            f.seek(256, 1)
            f.seek(4, 1)
            f.seek(bytesPerSecond, 1)

        elif (curRecord == startRecord) and (curRecord == stopRecord):
            # This is the only section we need to read
            headers = readHeaders(f, options)
            f.seek(4, 1)
            f.seek(startByte, 1)
            dataBytes = np.fromfile(f, dtype=np.uint8, count=int(stopByte-startByte))

        elif curRecord == startRecord:
            # This is the first section to read
            headers = readHeaders(f, options)
            f.seek(4, 1)
            f.seek(startByte, 1)
            dataBytes = np.fromfile(f, dtype=np.uint8, count=int(bytesPerSecond-startByte))

        elif (curRecord == stopRecord) and (stopByte > 0):
            # This is the last section we need to read
            headers = readHeaders(f, options)
            f.seek(4, 1)
            dataBytes = np.fromfile(f, dtype=np.uint8, count=int(stopByte))

        elif (curRecord > startRecord) and (curRecord < stopRecord):
            # This is one of the sections we need to read
            headers = readHeaders(f, options)
            f.seek(4, 1)
            dataBytes = np.fromfile(f, dtype=np.uint8, count=int(bytesPerSecond))

        # Check for EOF
        if 'EOF' in headers:
            break

        if len(dataBytes) > 0:
            # Get time
            tStr = '%04d-%03d' % (headers['SFDU_TIME_TAG_YEAR'], headers['SFDU_TIME_TAG_DOY'])
            t0 = datetime.datetime.strptime(tStr, '%Y-%j')
            t0 = t0 + datetime.timedelta(0,headers['SFDU_TIME_TAG_SECOND_OF_DAY'])

            # Get phase coefficients
            c0 = -(headers['SCHAN_PHASE_POLY_COEF_1']
                + headers['SCHAN_ACCUM_PHASE'])
            c1 = -headers['SCHAN_PHASE_POLY_COEF_2']
            c2 = -headers['SCHAN_PHASE_POLY_COEF_3']

            # Add bias if not starting at start of a second
            if curRecord == startRecord:
                dt0 = startByte/float(bytesPerSecond)
                t0  = t0 + datetime.timedelta(0,dt0)
            else:
                dt0 = 0

            # Determine number of words & samples
            nWords = len(dataBytes)/2.0
            nSamples = 16/(headers['SAMPLE_SIZE'])

            # Initialize data
            startByteData = 0
            colIndex = 1
            rowIndex = 1
            printStr = ''
            for ii in range(0, int(nWords)):
                endByteData = startByteData + 2 - 1
                word = dataBytes[startByteData:endByteData+1]

                # Convert to binary
                binVec = ''
                for index in range(0, len(word), 1):
                    binStr = "{0:b}".format(word[index])
                    binStr = binStr.zfill(8)
                    binVec = binVec + binStr

                for index in range(int(nSamples), 0, -1):
                    startBit = (index-1)*headers['SAMPLE_SIZE']
                    stopBit  = startBit + headers['SAMPLE_SIZE']

                    dataVal = bin2sint(binVec[startBit:stopBit]);
                    if colIndex == 1:
                        imagVal = dataVal
                        tVal    = t0 + datetime.timedelta(0,timeStep*(rowIndex-1))
                        td      = (tVal - t0)
                        dt      = ((td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10.0**6) + dt0
                        ph      = c0 + c1*dt + c2*(dt**2)
                        colIndex = 2

                    elif colIndex == 2:
                        realVal = dataVal
                        colIndex = 1
                        rowIndex = rowIndex + 1;

                        if (tVal >= options['startDate']) and (tVal < options['endDate']):
                            printStr = printStr + str(int(realVal)).rjust(8) + str(int(imagVal)).rjust(8)
                            output['times'].append(tVal)
                            output['values'].append(complex(realVal, imagVal))
                            output['phase'].append(ph)

                    if (len(printStr)>=96) and options['toPrintData']:
                        print(printStr)
                        printStr = ''

                # Increment start byte
                startByteData = endByteData + 1

            if (len(printStr)>0) and options['toPrintData']:
                print(printStr)

        # Increment record
        curRecord = curRecord + 1

    # Close file
    f.close()

    # Return outputs
    return output

# --------------- functions to run as stand-alone program ----------------------

def getopts(argv):
    fileName = ''
    opts = {}  # Empty dictionary to store key-value pairs.
    index = 0
    while index < len(argv):  # While there are arguments left to parse...
        if argv[index] == '--help':
            opts['--help'] = ''
            index = index + 1
            return fileName, opts
        elif argv[index] == '-p':
            opts['-p'] = ''
            index = index + 1
        elif argv[index] == '-v':
            opts['-v'] = ''
            index = index + 1
        elif argv[index] == '-q':
            opts['-q'] = ''
            index = index + 1
        elif argv[index] == '-d':
            opts['-d'] = argv[index+1]
            index = index + 2
        elif argv[index] == '-s':
            opts['-s'] = argv[index+1]
            index = index + 2
        else:
            fileName = argv[index]
            index = index + 1
    return fileName, opts


def main(DATAFIL, options):
    # Check for file type
    options['format'] = checkFormat(DATAFIL)

    #Run as main program
    if options['format'] == 'RDEF':
        data = RDEF(DATAFIL, options)
    elif options['format'] == 'VDR':
        data = VDR(DATAFIL, options)
    elif options['format'] == 'SFDU':
        data = SFDU(DATAFIL, options)

    return data

if __name__ == "__main__":
    # Default inputs
    options = {'toPlot':         False, 
               'toPrintData':    True, 
               'toPrintHeaders': False, 
               'lookForDate':    False}

    # Parse inputs
    [DATAFIL, myargs] = getopts(sys.argv)

    # If help mode prompted
    if ('--help' in myargs) or ('-h' in myargs) or (DATAFIL == '--help') or (DATAFIL == '-h'):
        print('')
        print('readFile.py:')
        print('\tThis function reads in a RDEF/VDR/SFDU data file.')
        print('')
        print('Syntax:')
        print('\tpython readFile.py fileName [ --help ] [ -v ] [ -p ] [ -q ] [ -d duration ]')
        print('\t                            [ -s "YYYY-MM-DD-HH:MM:SS.fff" ] [ -h ]')
        print('')
        print('Command Line Options:')
        print('\t--help : Displays help dialog')
        print('\t-v     : Includes header output at start of each record')
        print('\t-p     : Plots raw data values, histogram, and FFT spectrum')
        print('\t-q     : Suppresses all text output')
        print('\t-d     : Allows user to enter duration to parse (in seconds)')
        print('\t-s     : Allows user to enter start time (format: SS.fff or YYYY-MM-DD-HH:MM:SS.fff)')
        print('')

    else:
        if '-p' in myargs:
            options['toPlot'] = True

        if '-v' in myargs:
            options['toPrintHeaders'] = True

        if '-q' in myargs:
            options['toPrintHeaders'] = False
            options['toPrintData'] = False

        if '-s' in myargs:
            if myargs['-s'].find('-') != -1:
                options['startDate'] = datetime.datetime.strptime(myargs['-s'], '%Y-%m-%d-%H:%M:%S.%f')
                options['lookForDate'] = True
            else:
                try:
                    options['startTime'] = float(myargs['-s'])
                except:
                    raise ValueError('Start time input not recognized, use %Y-%m-%d-%H:%M:%S.%f or %S')

        if '-d' in myargs:
            options['duration'] = float(myargs['-d'])
        else:
            options['duration'] = float(99999)

        # Check for file type
        options['format'] = checkFormat(DATAFIL)

        #Run as main program
        if options['format'] == 'RDEF':
            data = RDEF(DATAFIL, options)
        elif options['format'] == 'VDR':
            data = VDR(DATAFIL, options)
        elif options['format'] == 'SFDU':
            data = SFDU(DATAFIL, options)

        # Plot if desired
        if options['toPlot']:
            plot(data)


