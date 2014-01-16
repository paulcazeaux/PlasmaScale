##
##	ES1 Makefile 
##
FILE_EXT = 
##
EXEC= plasmascale
##
CC = g++
##
##      C compiler.  Normally cc; is scc on Unicos.
#
#  This program needs libtk.a, libtcl.a, libXpm.a, and libXG20.a.
LIBDIRS = -L/opt/X11/lib -L/usr/local/opt/tcl-tk/lib
##
CFLAGS= -O3 -I/opt/X11/include -std=c++11 -lfftw
##
##      Flags used for compiling.  -O for optimization, which is optional.
##	If the X11 include files are not located in /usr/include/X11, you
##	must specify the directory as a flag prefixed by -I.  For instance, 
##	if the include files for X11 are located in /usr/local/include, the 
##	flag -I/usr/local/include must also be passed.
##
LIBS  = $(LIBDIRS) -lXGC250 -ltk8.6 -ltcl8.6 -lXpm -lX11 -lm -ldl
##
##      Libraries and their directories used in loading.  Normal is -lm and -lX11
##	for the Math and X11 libraries, which is NOT optional.  On Unicos, -lnet
##	is also used.  If the X11 libraries are not located in /usr/lib/X11, you
##	must also specify the directory as a flag prefixed by -L.  For instance,
##	if the include files for X11 are located in /usr/local/lib, the flag
##	-L/usr/local/lib must also be passed.  On Sun, use -Bstatic to avoid
##	problems with differing versions of OS.  Man cc for more information on flags.
##
##
PS1OBJ=	tools/Tools.o tools/RandomTools.o tools/Diagnostic.o tools/CurveDiagnostic.o tools/ScatterDiagnostic.o fields/Field.o fields/PlasmaFields.o particles/PopulationOfPopulations.o particles/CollectionOfParticles.o plasma/Plasma.o plasma/State.o plasma/History.o plasma/Simulation.o tools/XGrafixInterface.o PlasmaScale.o

all: $(EXEC)

.cpp.o:	PlasmaScale.h
	$(CC) -c $(CFLAGS) fields/*.cpp particles/*.cpp plasma/*.cpp tools/*.cpp $*.cpp 

$(EXEC):	$(PS1OBJ)
		$(CC) $(CFLAGS) -o $(EXEC) $(ES1OBJ) $(LIBS)

clean:
	@rm *.o
