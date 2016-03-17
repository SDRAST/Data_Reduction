#     FITS spectral line analysis
#
from Tkinter import *

import pyfits
import glob
import os
import os.path as P
from datetime import datetime
from time import strptime
import pylab as M
from matplotlib.ticker import FormatStrFormatter

x_default = "BF"
y_offset_default = 0.0

def isotime_to_object(iso):
    return datetime(*strptime(iso, "%Y-%m-%dT%H:%M:%S")[0:6])

def make_scan_list():
  global scanlist
  filenames = \
        glob.glob(os.getcwd()+'/*.fits')
  path,dummy = P.split(filenames[0])

  scanlist = []
  for f in filenames:
    hdulist = pyfits.open(f)
    h = hdulist[0].header
    source = "%12s" % (h['object'].strip())
    restfreq = "%7.3f" % (float(h['restfreq'])/1e9)
    vlsr = "%5.1f" % (float(h['velo-lsr'])/1e3)
    scan = int(h['scan-num'])
    integ = "%7.1f" % float(h['obstime'])
    [date_obs, fraction] = h['date-obs'].split('.')
    T = isotime_to_object(date_obs)
    date = T.date()
    time = T.time()
    if h['origin'] == 'CLASS-Grenoble':
        instrument = h['telescop'][4:]
    else:
        instrument = h['instrume']
    scanlist.append([scan,date,time,source,restfreq,vlsr,integ,instrument,f])
    
def whichSelected () :
    return select.curselection()

def plotFile() :
    # Get header info from first file
    x_par = x_parameter.get()
    offset = int(float(vertOffset.get()))
    scans = []
    for s in whichSelected():
        scans.append(int(s))
    try:
        scan,date,time,source,restfreq,vlsr,integ,instrument,f = \
                                                        scanlist[scans[0]]
        fitsfile = f
        hdulist = pyfits.open(f)
        h = hdulist[0].header
        head = h['object']+" with "+h['telescop']+" at "+h['date-obs']
        hdulist.close()
        ax = M.subplot(111)
        scan_list = []
        for s in scans:
            scan,date,time,source,restfreq,vlsr,integ,instrument,f = scanlist[s]
            scan_list.append(str(scan))
            fitsfile = f
            hdulist = pyfits.open(f)
            h = hdulist[0].header
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
                    print "Don't know about",h['telescop']
            else:
                print "Don't know how to handle files from",h['origin']
                sys.exit()
            scidata = hdulist[0].data
            spectrum = scidata[0][0][0]
            spec = map((lambda x:x+s*offset), spectrum)
            x = []
            for i in range(len(spectrum)):
                if x_par == 'BF':
                    value = (x0+(i-pix0)*deltax)/1e6
                elif x_par == 'RF':
                    value = (h['restfreq']+x0+(i-pix0)*deltax)/1e9
                elif x_par == 'VL':
                    value = (h['velo-lsr']+(i-pix0)*h['deltav'])/1e3
                x.append(value)
            hdulist.close()
            M.plot(x,spec)
            locs,labels = M.xticks()
            # M.setp(labels, "rotation", "vertical")
        if x_par == 'BF':
            M.xlabel('Frequency (MHz)')
            majorFormatter = FormatStrFormatter('%d')
        elif x_par == 'RF':
            M.xlabel('Frequency (GHz)')
            majorFormatter = FormatStrFormatter('%5.2f')            
        elif x_par == 'VL':
            M.xlabel('LSR Velocity (km/s)')
            majorFormatter = FormatStrFormatter('%5.1f')
        ax.xaxis.set_major_formatter(majorFormatter)
        M.legend(scan_list)
        M.ylabel('Antenna Temperature (K)')
        M.title(head)
        M.grid(True)
        M.connect('button_press_event',click)
        M.show()
    except IndexError:
        print "Select at least one scan"

def closePlot():
    M.close()

def clearSelection() :
    select.selection_clear(0,END)

def saveFigure():
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
    r1 = Radiobutton(o, text="Portrait",      variable=orientation,  value="portrait")
    r1.pack(side=LEFT)
    r2 = Radiobutton(o, text="Landscape", variable=orientation,  value="landscape")
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
    #r2.select()

    c = Frame(save)
    c.pack()
    bs = Button(c,text="Save", command=save_plot)
    bs.pack(side=LEFT)
    bc = Button(c,text="Cancel", command=cancel_save)
    bc.pack(side=LEFT)
    win.update_idletasks()

def terminate():
    M.close()
    try:
        save.destroy()
    except:
        pass
    win.destroy()
    win.quit()

def save_plot():
    # M.savefig(plotfile.get(),orientation=orientation,papertype='letter',format=figFormat)
    plotname=plotfilePath.get()+plotfile.get()+'.'+figFormat.get()
    orient = orientation.get()
    paper = 'letter'
    print "Saving in",plotname,"as",orient,"on",paper
    M.savefig(plotname,orientation=orient,papertype=paper)
    save.destroy()

def cancel_save():
    print "Cancelling Save Plot"
    save.destroy()
    
def show_par(*args):
    print "X parameter =",x_parameter.get()
    
def makeWindow () :
    global select
    global x_parameter
    global vertOffset
    
    win = Tk()
    win.wm_title("Main User Interface")

    frame1 = Frame(win)
    frame1.pack(anchor=W)

    xpar = Label(frame1,text="X axis:")
    xpar.grid(row=0,column=1)
    
    x_parameter = StringVar()
    x_parameter.set(x_default)
    
    xrb1 = Radiobutton(frame1, text="Baseband", variable=x_parameter, 
                       value="BF", command=show_par)
            
    xrb2 = Radiobutton(frame1, text="Rest Frequency", variable=x_parameter, 
                       value="RF", command=show_par)
            
    xrb3 = Radiobutton(frame1, text="LSR Velocity", variable=x_parameter, 
                       value="VL", command=show_par)
            
    xrb1.grid(row=0,column=2, sticky=W)
    xrb2.grid(row=0,column=3, sticky=W)
    xrb3.grid(row=0,column=4, sticky=W)

    Label(frame1, text="Vertical Offset").grid(row=1, column=0, sticky=W)
    vertOffset= StringVar()
    offset= Entry(frame1, textvariable=vertOffset)
    offset.grid(row=1, column=1, sticky=W)
    vertOffset.set(y_offset_default)

    frame2 = Frame(win)       # Row of buttons
    frame2.pack()
    b1 = Button(frame2,text=" Plot  ",command=plotFile)
    b2 = Button(frame2,text="Close Plot",command=closePlot)
    b3 = Button(frame2,text="Clear Selection",command=clearSelection)
    b4 = Button(frame2,text="Save Plot",command=saveFigure)
    b5 = Button(frame2,text=" Exit ",command=terminate)
    b1.pack(side=LEFT); b2.pack(side=LEFT)
    b3.pack(side=LEFT); b4.pack(side=LEFT)
    b5.pack(side=RIGHT)

    listboxhead = Frame(win)
    listboxhead.pack()
    head = Label(listboxhead, 
          text="   Date     Time   Scan      Source   Freq   Vlsr    sec  Instrument  ",
          font="Courier", anchor=W)
    head.pack()
    
    frame3 = Frame(win)       # select of names
    frame3.pack()
    scroll = Scrollbar(frame3, orient=VERTICAL)
    select = Listbox(frame3, yscrollcommand=scroll.set, 
                     width=70,height=12,font="Courier", 
                     selectmode=MULTIPLE)
    scroll.config (command=select.yview)
    scroll.pack(side=RIGHT, fill=Y)
    select.pack(side=LEFT,  fill=BOTH, expand=1)
    return(win)

def setSelect () :
    """Fill up the listbox"""
    scanlist.sort()
    select.delete(0,END)
    for scan,date,time,source,restfreq,vlsr,integ,instrument,f in scanlist :
        select.insert (END, \
                       str(date)+' '+str(time)+' '+ \
                       str(scan)+source+' '+str(restfreq)+' '+ \
                       str(vlsr)+' '+str(integ)+' '+instrument)

def click(event):
    """Act on a mouse click"""
    [xmin,xmax,ymin,ymax] = M.axis()
    if event.button == 1:
        xmin=event.xdata
        ymin=event.ydata
    if event.button == 2:
        print "The middle button was clicked"
    if event.button == 3:
        xmax=event.xdata
        ymax=event.ydata
    M.axis([xmin,xmax,ymin,ymax])
    
def make_plot_GUI():
  win = makeWindow()
  setSelect ()
  win.mainloop()

