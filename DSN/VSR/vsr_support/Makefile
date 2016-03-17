SRCDIR = .

FFTWDIR = /usr/include
FFTWLIB = /usr/lib/libfftw.a

SERVER = /usr/local/bin

INCS = -I$(SRCDIR) -I/usr/include/python2.6

CC = gcc -Wall
CC_SWITCHES =   -g -c -fPIC
SHLIB_LD =      /usr/bin/ld -G -z text
SHLIB_SUFFIX =  .so

CFLAGS = -ggdb -DUNIX -DCDATE="\"`date`\"" $(INCS) 
LD = gcc

PYINC = /usr/include/python2.6
PYLIB = /usr/bin


vdr_support.so: vdr_support.c
	gcc vdr_support.c  ${CC_SWITCHES} \
	-I$(PYINC) -L$(PYLIB)
	${SHLIB_LD} vdr_support.o \
	-lpython2.6 -o vdr_support.so

clean: 
	rm -f *.o;
	rm -f *.a;
	rm -f *.exe;
	rm -f *.FITS;
	rm -f $(EXECS);

$(UTILLIB): FORCE
	cd $(UTILDIR); $(MAKE) all

FORCE:

# use gcc -MM to generate non-system include file dependencies
# use sed to add ".depend" as a target which depends on ALL included files
.depend: $(SRCS) Makefile
	$(CC) -MM $(CFLAGS) $(SRCS) | sed 's/^\([^# ]\)/.depend \1/' >$@

include .depend
