This code interactively plots the Appleton-Hartree dispersion relation
governing waves in a cold plasma. Its purpose is primarily
pedagogical.  It is controlled from the keyboard, for which help can
be obtained when it is running by hitting the key 'h'.

Easiest is just to download and run the Windows complete executable:
coldplas.exe. It runs perfectly in native Windows and under the 'Wine'
system.

On a Mac, to run coldplas.exe you can either install Wine free (see e.g
https://www.davidbaumgold.com/tutorials/wine-mac/) or WineBottler
(http://winebottler.kronenberg.org/) or use a packaged app for running
Windows programs such as WinOnX
(https://apps.apple.com/us/app/winonx-64/) or Crossover
(https://www.codeweavers.com/products).

To create an executable program from source on linux type make, but
you need gfortran and the X11 libraries accessible.  The makefile with
'make windows' is set up to use the mingw-32 cross-compiler to produce
a Windows executable coldplas.exe.  Edit the makefile to replace MINGW
or FORTRAN with references to your compiler. But it must be compatible
with the gfortran/gcc compilers that produced the required and
included graphics libaries libaccisX.a or libaccisWin.a.



No guarantee of fitness for any purpose whatever is given, and anyone
using the program does so at their own risk.
Copyright (C) Ian H Hutchinson 2019.
