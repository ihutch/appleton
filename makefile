# Put the required graphics library libaccisX.a in directory ./accis/ 
ACCISPARENT=./
export ACCISPARENT
# Then this will get and compile, libaccisX.a in ./accis/
include ACCIS.mk
# $(LIBPATH) $(LIBRARIES) and $(LIBDEPS) are inherited for below.

XLIBS = $(LIBPATH)
LIBS = $(LIBRARIES)
MINGW=i686-w64-mingw32-gfortran
MINGW-SWITCHES= -H -mwindows -mconsole --static
MINGWLIBS=-L. -laccisWin
# -02 provoked aggressive optimization warnings and undefined behavior.
COMPILE-SWITCHES = -Wall -O1 

##########################################################################
% : %.f90 makefile;
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f90  $(XLIBS) $(LIBS)

# The windows version depends on having libaccisWin.a
# It does not normally need to be remade because the executable just works. 
%.exe : %.f90 makefile;
	$(MINGW) -o $*.exe $(MINGW-SWITCHES) $*.f90 $(MINGWLIBS)
##########################################################################

default : coldplas $(LIBDEPS)

windows : coldplas.exe libaccisWin.a

renewlibs :
#	cp /home/hutch/accis/libaccisX.a .
	cp /home/hutch/accis/drivers/win64/libaccisWin.a .

clean :
	rm -f coldplas
	rm -fr accis
