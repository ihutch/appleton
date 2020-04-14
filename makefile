
FORTRAN=gfortran
XLIBS = -L. -laccisX -lX11
MINGW=i686-w64-mingw32-gfortran
MINGW-SWITCHES= -H -mwindows -mconsole --static
MINGWLIBS=-L. -laccisWin
# -02 provoked aggressive optimization warnings and undefined behavior.
COMPILE-SWITCHES = -Wall -O1 

##########################################################################
% : %.f90 makefile libaccisX.a;
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f90  $(XLIBS) $(LIBS)

%.exe : %.f90 makefile;
	$(MINGW) -o $*.exe $(MINGW-SWITCHES) $*.f90 $(MINGWLIBS) $(LIBS)
##########################################################################

default : coldplas libaccisX.a

windows : coldplas.exe libaccisWin.a

renewlibs :
	cp /home/hutch/accis/libaccisX.a .
	cp /home/hutch/accis/drivers/win64/libaccisWin.a .

