"""FITS spectral line analysis
A directory window provides a list of the FITS files in the currently
selected directory.  A subset can be selected for plotting.  The x-axis
units can be selected by the user.
The functions defined here are:
plotFile        - plot the selected scans
saveFigure      - create a toplevel to save a figure to a file
terminate       - to exit this program
save_plot       - save the plot to a file
cancel_save     - to cancel saving a plot and remove the plot save toplevel
show_par        - print a parameter being changed
update_data_dir - change the data directory"""

fits_gui_debug = 1

from Tkinter import *

import pyfits
import os
import os.path as P

# The next two are part of this package
import Data_Reduction.FITS.fits_support as FI
import file_selection as FS
from Astronomy         import *
from datetime          import datetime
from time              import strptime
from matplotlib        import rc
from matplotlib.ticker import FormatStrFormatter
from SimpleDialog      import SimpleDialog

import pylab as M # must be after importing matplotlib.rc
import webbrowser
  
x_default = "RF"
y_offset_default = 0.0

helpfile = "/usr/local/lib/python2.4/site-packages/Data_Reduction/index.html"
aboutfile = "/usr/local/lib/python2.4/site-packages/Data_Reduction/about.html"

def check_browser():
  """This checks to see which browser is available on this host."""
  webbrowser._tryorder = \
    ['iceweasel %s &','firefox %s &','mozilla-firefox %s &','mozilla %s &']
  browseropen = False

def select_files():
    """This select scans from current data directory and adds them to the
    selected scan list."""
    #----------------
    global list_of_scans
    global scans
    global selected_listbox

    scans = []         # indices of selected scans in upper window
    if fits_gui_debug:
        print "\nfits_gui.select_files entered with list_of_scans =\n",\
              list_of_scans
    # 'scans' is an integer list of the indices of the selected scans to be
    # plotted
    for s in FS.whichSelected():
        scans.append(int(s))
    if fits_gui_debug:
        print "\nselect_files; scans selected: ",scans
    # Trap the case of no-scans-selected
    try:
        # get header data from the first selected scan
        scan,date,time,source,restfreq,vlsr,integ,instrument,f = \
                                                        list_of_scans[scans[0]]
        fitsfile = f
        hdulist = pyfits.open(f)
        h = hdulist[0].header
        try:
            source = h['object']
        except KeyError:
            source = "Unknown source"
            print "plotFile, no source; Header:\n",h
        hdulist.close()
    except IndexError:
        print "plotFile: Select at least one scan"
        return
    except KeyError:
        print "plotfile: Header:\n",h
        hdulist.close()
        raise KeyError, 'Keyword %s not found.' % `key`
    for s in scans:
        scan,date,time,source,restfreq,vlsr,integ,instrument,f \
            = list_of_scans[s]
        scan_files.append(f)
        selected_listbox.insert (END, \
                    str(date)+' '+str(time)+' '+ \
                    str(scan)+source+' '+str(restfreq)+' '+ \
                    str(vlsr)+' '+str(integ)+' '+instrument)
    win.update_idletasks()

def clear_selected():
    """Clears the scan list""" 
    global scans
    scans_files = []
    selected_listbox.delete(0,END)
    
def plotFile():
    """This plots the selected scans with the specified x-axis units and an
    optional offset between scans in a window. The window has built-in
    methods for panning, zooming, recalling previous plots, etc."""
    #----------------
    global list_of_scans
    global scans
    
    if fits_gui_debug:
        print "\nfits_gui.plotFile entered with scans =\n",scans
    x_par = x_parameter.get()
    offset = int(float(vertOffset.get()))
    ax_bottom = M.subplot(111)
    if fits_gui_debug:
        print "plotFile: selected scans:",scans
    scan_labels = []
    s = -1 # plotted scans labelled 0, 1, ...
    for fits_file in scan_files:
        s += 1
        if fits_gui_debug:
            print "plotFile processing scan",s
        hdulist = pyfits.open(fits_file)
        # Get header data
        h = hdulist[0].header
        if s == 0:
            source, restfreq, vlsr, scan, integ, date, time, instrument \
                = FS.get_FITS_header(fits_file)
        # See if BANDWID exists
        try:
            bandwidth = h["bandwid"]
        except KeyError:
            # At CSO, and maybe others, the BANDWID keyword is not used
            if h['origin'] == 'CLASS-Grenoble':
                n_chan = int(h['naxis1'])
                deltax = float(h['cdelt1'])
                x0 = float(h['crval1'])
                pix0 = float(h['crpix1'])
                if h['telescop'] == 'CSO 4GHZ IF3':
                    bandwidth = 1000
                elif h['telescop'] == 'CSO 50MHZ':
                    bandwidth = 50
                else:
                    print "plotFile: Don't know about",h['telescop']
            else:
                print "plotFile: Don't know how to handle files from",h['origin']
                sys.exit()
        try:
            scan_labels.append(str(h['scan-num']))
        except KeyError:
            scan_labels.append('')
                               
        # Get spectral data
        scidata = hdulist[0].data
        spectrum = scidata[0][0][0]
        spec = map((lambda x:x+s*offset), spectrum)
        x = []
        if x_par == 'RF':
            # h['velo-lsr'] is in m/s
            dopfac=1.0-h['velo-lsr']/3.e8
            #print "plotFile: Vlsr =%6.1f, dopfac=%f" % (h['velo-lsr'],dopfac)
        for i in range(len(spectrum)):
            if x_par == 'BF':
                # This is the baseband frequency in MHz, at least at CSO
                value = (x0+(i-pix0)*deltax)/1e6
            elif x_par == 'RF':
                # This is the frequency in the LSR in GHz
                # for the upper sideband.  For the lower sideband
                value = (h['restfreq']+x0+(i-pix0)*deltax)*dopfac/1e9
            elif x_par == 'VL':
                # This is the velocity in the LSR in km/s
                value = (h['velo-lsr']+(i-pix0)*h['deltav'])/1e3
            x.append(value)
        hdulist.close()
        M.plot(x,spec)
        locs,labels = M.xticks()
    # end for
    if len(scan_files) == 1:
        head = source+" with "+h['telescop']+" at "+h['date-obs']
    else:
        head = source+" with "+h['telescop']
    # M.setp(labels, "rotation", "vertical")
    if x_par == 'BF':
        M.xlabel('Frequency (MHz)')
        majorFormatter = FormatStrFormatter('%d')
    elif x_par == 'RF':
        M.xlabel('Frequency (GHz)')
        majorFormatter = FormatStrFormatter('%6.3f')            
    elif x_par == 'VL':
        M.xlabel('LSR Velocity (km/s)')
        majorFormatter = FormatStrFormatter('%5.1f')
    ax_bottom.xaxis.set_major_formatter(majorFormatter)
    labels = []
    for s in scans:
        labels.append(str(s))
    M.legend(scan_labels)
    M.ylabel('Antenna Temperature (K)')
    M.figtext(0.15,0.85,head,horizontalalignment='left') # custom title location
    M.grid(True)
    # M.connect('button_press_event',click)
    if x_par == 'RF':
        rest_freq=h['RESTFREQ']/1e9
        try:
            imag_freq=h['IMAGFREQ']/1e9
            band_sep = rest_freq-imag_freq
            if fits_gui_debug:
                print rest_freq, imag_freq, band_sep
            ax_bottom_limits = M.xlim()
            if fits_gui_debug:
                print "Axis limits are",ax_bottom_limits
            left_rest_freq = ax_bottom_limits[0]
            right_rest_freq = ax_bottom_limits[1]
            plot_width = right_rest_freq-left_rest_freq
            if fits_gui_debug:
                print plot_width
            left_imag_freq = left_rest_freq-band_sep
            right_imag_freq = right_rest_freq - band_sep - 2.*plot_width
            if fits_gui_debug:
                print left_imag_freq, right_imag_freq
            ax_top = twiny()
            ax_top_limits = M.xlim([left_imag_freq, right_imag_freq])
            M.xlabel('Image Frequency (GHz)')
            majorFormatter = FormatStrFormatter('%5.2f')            
            ax_top.xaxis.set_major_formatter(majorFormatter)
            ax_top.xaxis.tick_top()
        except:
            pass
    M.show()

    
def print_header():
    """Prints header data summary of all the current scans in a pretty
    format."""
    scans = whichSelected()
    for s in scans:
        scan,date,time,source,restfreq,vlsr,integ,instrument,f = \
                                                    list_of_scans[int(s)]
        fitsfile = f
        hdulist = pyfits.open(f)
        h = hdulist[0].header
        print h['OBJECT'],\
              decimal_hrs_to_hms(h['CRVAL2']/15,1),\
              decimal_deg_to_dms(h['CRVAL3'],1),h['EQUINOX']
        print "%6.1f s integration at %s" % (h['OBSTIME'],h['DATE-OBS'])
        print h['LINE'],"freq=",h['RESTFREQ']/1e9,"  imag=",h['IMAGFREQ']/1e9
        print "Tsys=%4f, image gain=%4.2f" % (h['TSYS'],h['GAINIMAG'])
        print "az=%5.1f el=%5.1f tau=%5.3f" % \
              (h['AZIMUTH'],h['ELEVATIO'],h['TAU-ATM'])
    
def closePlot():
    """Closes the plot window."""
    M.close()

def saveFigure():
    """This creates a toplevel to specify plot output options and
    file name and location."""
    global save
    global plotfilePath
    global plotfile
    global orientation
    global figFormat
    save = Toplevel()
    save.wm_title("Save Plot")

    p = Frame(save)
    p.pack(side=TOP)
    plbl = Label(p, text="Path:")
    plbl.pack(side=LEFT)
    plotfilePath = StringVar()
    plotfilePath.set("/tmp/")
    pval = Entry(p,textvariable=plotfilePath, width=50)
    pval.pack(side=LEFT)

    f = Frame(save)
    f.pack(side=TOP,anchor=W)
    flbl = Label(f, text="File:")
    flbl.pack(side=LEFT)
    plotfile = StringVar()
    defaultName = datetime.strftime(datetime.now(),"%Y%m%d-%H%M%S")
    plotfile.set(defaultName)
    fval = Entry(f, textvariable=plotfile, width=20)
    fval.pack(side=LEFT)

    o = Frame(save)
    o.pack(side=TOP,anchor=W)
    orientation = StringVar()
    orientation.set("landscape")
    r1 = Radiobutton(o, text="Portrait", variable=orientation, \
                     value="portrait")
    r1.pack(side=LEFT)
    r2 = Radiobutton(o, text="Landscape", variable=orientation, \
                     value="landscape")
    r2.pack(side=LEFT)

    # Agg can handle svg, ps, pdf
    figFormat = StringVar()
    figFormat.set("eps")
    r3 = Radiobutton(o, text="EPS", variable=figFormat, value="eps")
    r3.pack(side=LEFT)
    r4 = Radiobutton(o,text="PS",variable=figFormat,value="ps")
    r4.pack(side=LEFT)
    r5 = Radiobutton(o,text="SVG",variable=figFormat,value="svg")
    r5.pack(side=LEFT)
    r6 = Radiobutton(o,text="PDF",variable=figFormat,value="pdf")
    r6.pack(side=LEFT)

    c = Frame(save)
    c.pack()
    bs = Button(c,text="Save", command=save_plot)
    bs.pack(side=LEFT)
    bc = Button(c,text="Cancel", command=cancel_save)
    bc.pack(side=LEFT)
    win.update_idletasks()
    # end saveFigure

def terminate():
    """Ends the program."""
    global win
    M.close()
    try:
        save.destroy()
    except:
        pass
    win.destroy()
    win.quit()

def save_plot():
    """Saves the current plot to a file."""
    plotname=plotfilePath.get()+plotfile.get()+'.'+figFormat.get()
    orient = orientation.get()
    paper = 'letter'
    print "Saving in",plotname,"as",orient,"on",paper
    M.savefig(plotname,orientation=orient,papertype=paper)
    save.destroy()

def cancel_save():
    """Cancel saving a plot and remove the plot save window"""
    if fits_gui_debug:
        print "Cancelling Save Plot"
    save.destroy()
    
def show_par(*args):
    """print a parameter whose value is being changed"""
    print "new parameter value =",x_parameter.get()

def update_data_dir(*args):
    """Invoked by change of data_directory"""
    global list_of_scans
    if fits_gui_debug:
        print "\nfits_gui.update_data_dir: Old list of scans =\n",\
          list_of_scans[0],'...'
    FS.NEW_DATA_DIR()
    list_of_scans = FS.report_scanlist()
    if fits_gui_debug:
        print "\nfits_gui.update_data_dir: New list of scans =\n",\
          list_of_scans[0],'...'
    return list_of_scans
    
def makeWindow (DataDir) :
    """This makes the main user interface"""
    global x_parameter
    global vertOffset
    global data_directory
    global x_default
    global y_offset_default
    global list_of_scans
    global helpfile
    global wbfd
    global selected_listbox
    
    # win = Tk()
    win.wm_title("Main User Interface")

    menubar = Menu(win)
    
    filemenu = Menu(menubar, tearoff=0)
    filemenu.add_command(label="Exit", command=terminate)
    menubar.add_cascade(label="File", menu=filemenu)
    
    helpmenu = Menu(menubar, tearoff=0)
    action1 = lambda url=helpfile: get_help(url)
    helpmenu.add_command(label="Overview", command=action1)
    action2 = lambda url=aboutfile: get_help(url)
    helpmenu.add_command(label="About", command=action2)
    menubar.add_cascade(label="Help", menu=helpmenu)

    win.config(menu=menubar)

    # Frame to show/select data directory, from module file_selection
    data_directory = FS.MAKE_DIRECTORY_FRAME(win, DataDir)
    # update list_of_scans when data_directory changes
    data_directory.trace("w",update_data_dir)

    frame2 = Frame(win)       # Row of buttons
    frame2.pack()
    b0 = Button(frame2,text="Header",command=print_header)
    b1 = Button(frame2,text="Select",command=select_files)
    b3 = Button(frame2,text="Clear Selection",command=FS.clearSelection)
    b0.pack(side=LEFT)
    b1.pack(side=LEFT)
    b3.pack(side=LEFT)

    selectboxhead = Frame(win)
    selectboxhead.pack()
    select_head = Label(selectboxhead, 
          text="   Date     Time   Scan      Source   Freq   Vlsr    sec  Instrument  ",
          font="Courier", anchor=W)
    select_head.pack()
    select_frame = FS.MAKE_FILE_FRAME(win,"selection")
    
    # Widget to select X coordinate type
    frameXtype = Frame(win)
    frameXtype.pack(anchor=W)

    xpar = Label(frameXtype,text="X axis:")
    xpar.grid(row=0,column=1)
    
    x_parameter = StringVar()
    x_parameter.set(x_default)
    
    xrb1 = Radiobutton(frameXtype, text="Baseband", variable=x_parameter, 
                       value="BF", command=show_par)
    xrb2 = Radiobutton(frameXtype, text="Rest Frequency", variable=x_parameter, 
                       value="RF", command=show_par)        
    xrb3 = Radiobutton(frameXtype, text="LSR Velocity", variable=x_parameter, 
                       value="VL", command=show_par)           
    xrb1.grid(row=0,column=2, sticky=W)
    xrb2.grid(row=0,column=3, sticky=W)
    xrb3.grid(row=0,column=4, sticky=W)

    frameY= Frame(win)
    frameY.pack(side=TOP,anchor=W)
    ylbl = Label(frameY,text="Y axis:")
    ylbl.pack(side=LEFT,anchor=W)
    ypar = Label(frameY, text="Vertical Offset")
    ypar.pack(side=LEFT,anchor=W)
    vertOffset= StringVar()
    offset= Entry(frameY, textvariable=vertOffset,width=5)
    offset.pack(side=LEFT,anchor=W)
    vertOffset.set(y_offset_default)

    frame3 = Frame(win)       # Another row of buttons
    frame3.pack()
    b1 = Button(frame3,text=" Plot  ",command=plotFile)
    b4 = Button(frame3,text="Save Plot",command=saveFigure)
    b0 = Button(frame3,text="Clear Selected",command=clear_selected)
    b1.pack(side=LEFT)
    b4.pack(side=LEFT)
    b0.pack(side=LEFT)

    selectedboxhead = Frame(win)
    selectedboxhead.pack()
    selected_head = Label(selectedboxhead, 
          text="                               Selected Files                         ",
          font="Courier", anchor=W)
    selected_head.pack()
    selected_listbox = FS.MAKE_FILE_FRAME(win,"selected")
    return win 

def click(event):
    """Acts on a mouse click"""
    print event.name,"in",event.canvas
    print "button=",event.button,"x =",event.x,"y =",event.y
    print "Data values: x =",event.xdata,"y =",event.ydata
    print "key:",event.key
    [xmin,xmax,ymin,ymax] = M.axis()
    if event.button == 1:
        xmin=event.xdata
        ymin=event.ydata
        if fits_gui_debug:
            print "Left lower corner is now %f,%f" % (xmin,ymin)
    if event.button == 2:
        if fits_gui_debug:
            print "The middle button was clicked"
    if event.button == 3:
        xmax=event.xdata
        ymax=event.ydata
        if fits_gui_debug:
            print "Upper right corner is now %f,%f" % (xmax,ymax)
    # M.axis([xmin,xmax,ymin,ymax])

def get_help(url):
    """This opens a web browser at the specified URL or else changes the
    URL of an already open browser."""
    SimpleDialog(win,
             text="If there is a browser open software,\n\
        the help text will be there.",
             buttons=["OK"],
             default=0,
             title="Help Notice").go()
    webbrowser.open(url)
    
        
def FITS_GUI(defaultDataDir):
    """This creates the user interface in the main window"""
    global list_of_scans
    global win
    win = makeWindow(defaultDataDir)
    list_of_scans = FI.list_FITS_files(defaultDataDir) # from module fits_support
    FS.FILL_FILE_LISTBOX(list_of_scans) # from module file_selection
    win.mainloop()

def start_GUI():
  """This initializes various lists, creates a Tk window, and makes a
  FITS GUI in it.  For the initial data directory it uses
  /home/kuiper/Projects/Methanol/Observations/CSO/FITS/"""
  global list_of_scans
  global scan_files
  list_of_scans = [] # scans in scan selection (upper) window
  scan_files = []    # List filenames of selected scans (lower window)

  win = Tk()

  FITS_GUI('/home/kuiper/Projects/Methanol/Observations/CSO/FITS/')
