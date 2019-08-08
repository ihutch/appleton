
FORTRAN=gfortran
XLIBS = -L. -laccisX -lXt -lX11
MINGW=i686-w64-mingw32-gfortran
MINGW-SWITCHES= -H -mwindows -mconsole --static
MINGWLIBS=-L. -laccisWin
COMPILE-SWITCHES = -Wall -O2

##########################################################################
% : %.f90 makefile;
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f90  $(XLIBS)

%.exe : %.f90 makefile;
	$(MINGW) -o $*.exe $(MINGW-SWITCHES) $*.f90 $(MINGWLIBS)
##########################################################################

default : coldplaslog.exe libaccisWin.a

LIBS: libaccisX.a libaccisWin.a

libaccisX.a :
	cp /home/hutch/accis/libaccisX.a .

libaccisWin.a :
	cp /home/hutch/accis/drivers/win64/libaccisWin.a .
