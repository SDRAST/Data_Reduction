"""
Tk GUIs for file selection

GUIs to select directories and files for data reduction.
This is an obsolete module since we've switched to Qt.
"""

import glob
import os
from Tkinter      import *
from tkMessageBox import *
from Data_Reduction.FITS.fits_support import *
from tkFileDialog import askdirectory

file_selection_debug = 0


def whichSelected () :
    """
    Returns the cursor selection(s)
    """
    # file_select is the name of the listbox with the files from which
    # selections are made
    global file_select
    return file_select.curselection()

def clearSelection():
    """
    Clears the file selections
    """
    # file_select is the name of the listbox with the files from which
    # selections are made
    global file_select
    file_select.selection_clear(0,END)

def NEW_DATA_DIR():
    """
    Open a new data directory
    
    Invoked by the [New] button on the main GUI.  Also by update_data_dir()
    in module fits_gui which is a trace function for variable
    'data_directory'.
    'data_directory' is returned from MAKE_DIRECTORY_FRAME.
    Opens a directory selection window; returns as globals (accessible only
    to this module) the new directory -dataDir- and the list of FITS files
    -scanlist- in the directory
    """
    global dataDir
    global scanlist

    newDir=dataDir.get() # initialize newDir
    if file_selection_debug:
        print "\nNEW_DATA_DIR: Starting from",newDir
    n_files = 0
    # will loop until num. FITS files in newDir is > 0.
    while n_files < 1:
        newDir = askdirectory(initialdir=newDir)
        if file_selection_debug:
            print "\nfile_selection.NEW_DATA_DIR: Checking\n",newDir
        if newDir != ():
            filenames = glob.glob(newDir+'/*.fits')
            if file_selection_debug:
                print "\nfile_selection.NEW_DATA_DIR: found files:\n",\
                      filenames[0],'...'
            n_files = len(filenames)
        else:
            n_files = 0
        if n_files < 1:
            showinfo(title="No FITS",message="No FITS files available")
            os.chdir(newDir+"/..")
    if file_selection_debug:
        print "\nfile_selection.NEW_DATA_DIR: directory:\n",newDir
        print "\nfile_selection.NEW_DATA_DIR: files:\n",filenames[0],'...'
    showinfo(title="FITS found",message=str(n_files)+" FITS files")
    os.chdir(newDir)
    dataDir.set(newDir)
    scanlist = list_FITS_files(newDir)
    if file_selection_debug:
        print "\nfile_selection.NEW_DATA_DIR: scanlist=\n",scanlist[0],'...'
    FILL_FILE_LISTBOX(scanlist)

def report_scanlist():
    """
    Get the scan list

    This simply returns the current scan list.  Reference to
    'file_selection.scanlist' would work just as well.
    """
    global scanlist
    if file_selection_debug:
        print "\nfile_selection.report_scanlist reports scanlist:\n",scanlist[0],'...'
    return scanlist

def MAKE_DIRECTORY_FRAME(win, defaultDataDir):
    """
    Make a current directory frame

    Creates a frame showing the current data directory and a button
    for selecting a new data directory
    """
    global dataDir # for local use in this module
    
    frameDir = Frame(win)
    frameDir.pack(anchor=W)

    dataDirList = Button(frameDir,text="List FITS dirs",command=SELECT_FITS_DIR)
    dataDirList.pack(side=LEFT,anchor=W)
    
    dataDirLbl = Label(frameDir,text="Data directory:")
    dataDirLbl.pack(side=LEFT,anchor=W)
    
    dataDirNew = Button(frameDir,text="New",command=NEW_DATA_DIR)
    dataDirNew.pack(side=LEFT,anchor=W)
    
    dataDir = StringVar(win)
    dataDir.set(defaultDataDir)
    dataDirName = Entry(frameDir,width=50,textvariable=dataDir)
    dataDirName.pack(side=LEFT,anchor=W)
    return dataDir
    
def MAKE_FILE_FRAME(win,contents):
    """
    Creates a frame to hold the names of the FITS files.
     
    This hold names of FITS files in the current data
    directory. FILL_FILE_LISTBOX is used to populate the listbox
    with the data
    in global variable 'scanlist'
    """
    global file_select
    global main_win
    main_win = win
    file_frame = Frame(win)       # select of names
    file_frame.pack()
    scroll = Scrollbar(file_frame, orient=VERTICAL)
    if contents == "selection":
        height=12
    elif contents == "selected":
        height=6
    else:
        print "MAKE_FILE_FRAME called for unknown contents:",contents
        return
    box = Listbox(file_frame, yscrollcommand=scroll.set, 
                     width=70,height=height,font="Courier",bg="white",
                     selectmode=MULTIPLE)
    scroll.config (command=box.yview)
    scroll.pack(side=RIGHT, fill=Y)
    box.pack(side=LEFT,  fill=BOTH, expand=1)
    if contents == "selection":
        file_select = box
    else:
        selected_file_box = box
    return box
    
def MAKE_FITS_SELECT_DIR():
    """
    Create a toplevel Listbox for selecting a data directory
    """
    global dirwin
    global dir_select
    global directories
    dirwin = Toplevel()
    dirwin.wm_title("Select FITS directory")
    dir_frame = Frame(dirwin)       # select of names
    dir_frame.pack()
    dir_scroll = Scrollbar(dir_frame, orient=VERTICAL)
    dir_select = Listbox(dir_frame, yscrollcommand=dir_scroll.set, 
                     width=70,height=12,font="Courier")
    dir_scroll.config (command=dir_select.yview)
    dir_scroll.pack_configure(side=RIGHT, fill=Y)
    dir_select.pack_configure(side=LEFT,  fill=BOTH, expand=1)
    done = Button(dirwin,text="Select",command=SELECT_NEW_FITS_DIR)
    done.pack()
    
    # Find all directories with FITS files
    directories = find_FITS_dirs()
    # Fill up the listbox
    dir_select.delete(0,END)
    for directory in directories:
        # print "Inserting",directory
        dir_select.insert (END,directory)
    if file_selection_debug:
        print "MAKE_FITS_DIR_SELECT: dir_select loaded"
    # end MAKE_FITS_SELECT_DIR
    
def SELECT_FITS_DIR():
    """
    This just invokes MAKE_FITS_SELECT_DIR().

    It's not clear what
    the point of this is.
    """
    global dirwin
    global dir_select
    global main_win
    MAKE_FITS_SELECT_DIR()
    if file_selection_debug:
        print "SELECT_FITS_DIR: mainloop started"
    main_win.update_idletasks()
    if file_selection_debug:
        print "MAKE_FITS_DIR_SELECT: main_win updated"

def SELECT_NEW_FITS_DIR():
    """Another mysterious function.

    It's used to select a directory,
    not a file."""
    global dir_select
    global dirwin
    global dataDir
    global directories
    selected = dir_select.curselection()
    if file_selection_debug:
        print "SELECT_NEW_FITS_DIR: selected",selected
    # selected is a sequence
    index = int(selected[0])
    new_dir = directories[index]
    print "SELECT_NEW_FITS_DIR: name",new_dir
    dataDir.set(new_dir)
    dirwin.destroy()
    return
    # end SELECT_NEW_FITS_DIR
    
def FILL_FILE_LISTBOX(scans):
    """
    Fill up the file selection listbox
    """
    global file_select
    global scanlist
    scanlist = scans # makes argument from other module global here
    if file_selection_debug:
            print "\nfile_selection.FILL_FILE_LISTBOX entered with scanlist:\n",\
                  scanlist[0],'...'
    if len(scanlist) > 0:
        if file_selection_debug:
            print "\nfile_selection.FILL_FILE_LISTBOX processing scanlist:\n",\
                  scanlist[0],'...'
        file_select.delete(0,END)
        for scan,date,time,source,restfreq,vlsr,integ,instrument,f in scanlist :
            file_select.insert (END, \
                       str(date)+' '+str(time)+' '+ \
                       str(scan)+source+' '+str(restfreq)+' '+ \
                       str(vlsr)+' '+str(integ)+' '+instrument)
    else:
        pass
