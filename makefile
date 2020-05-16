# Put the required graphics library libaccisX.a in directory ./accis/ 
ACCISPARENT=./
export ACCISPARENT
# Then this will get and compile, libaccisX.a in ./accis/
include ACCIS.mk
# $(LIBPATH) $(LIBRARIES) and $(LIBDEPS) are inherited for below.

XLIBS = $(LIBPATH)
LIBS = $(LIBRARIES)
# -02 provoked aggressive optimization warnings and undefined behavior.
COMPILE-SWITCHES = -Wall -Wno-integer-division -O1

# Cross-compiler only. Needs libaccisWin.a in this directory.
MINGW=i686-w64-mingw32-gfortran
MINGW-SWITCHES= -H -mwindows -mconsole --static
MINGWLIBS=-L. -laccisWin
##########################################################################
# Make syntax $(if ..., ..., ...) used to opt out if no compiler.
% : %.f90 makefile ACCIS.mk;
	$(if $(FORTRAN),\
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f90  $(XLIBS) $(LIBS),\
	@echo "No FORTRAN compiler found"; exit 1;)
# Compiling the windows version depends on having libaccisWin.a
# It does not normally need to be remade because the executable just works. 
%.exe : %.f90 makefile;
	$(MINGW) -o $*.exe $(MINGW-SWITCHES) $*.f90 $(MINGWLIBS)
##########################################################################

default : coldplas $(LIBDEPS)

windows : coldplas.exe libaccisWin.a

renewlibs :
	cp /home/hutch/accis/drivers/win64/libaccisWin.a .

clean :
	rm -f coldplas
	rm -fr accis
