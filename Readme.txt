This code interactively plots the Appleton-Hartree dispersion relation
governing waves in a cold plasma. Its purpose is primarily
pedagogical.  It is controlled from the keyboard, for which help can
be obtained when it is running by hitting the key 'h'.

To create an executable program, a fortran compiler must be available.
The makefile is set up for the mingw cross-compiler to produce a
Windows executable coldplaslog.exe, or for gfortran to produce a linux
executable coldplaslog. Edit the makefile to replace MINGW or FORTRAN
with references to your compiler. But it must be compatible with the
gfortran/gcc compilers that produced the required and included
graphics libaries libaccisX.a or libaccisWin.a.

Alternatively just download the Windows executable: coldplaslog.exe
and run it.

No guarantee of fitness for any purpose whatever is given, and anyone
using it does so at their own risk.
Copyright (C) Ian H Hutchinson 2019.
