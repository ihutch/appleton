
FORTRAN=gfortran
XLIBS = -L. -laccisX -lXt -lX11
MINGW=i686-w64-mingw32-gfortran
MINGW-SWITCHES= -H -mwindows -mconsole --static
MINGWLIBS=-L. -laccisWin
# -02 provoked aggressive optimization warnings and undefined behavior.
COMPILE-SWITCHES = -Wall -O1

##########################################################################
% : %.f90 makefile;
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f90  $(XLIBS)

%.exe : %.f90 makefile;
	$(MINGW) -o $*.exe $(MINGW-SWITCHES) $*.f90 $(MINGWLIBS)
##########################################################################

default : coldplas.exe libaccisWin.a

LIBS: libaccisX.a libaccisWin.a

libaccisX.a :
	cp /home/hutch/accis/libaccisX.a .

libaccisWin.a :
	cp /home/hutch/accis/drivers/win64/libaccisWin.a .
