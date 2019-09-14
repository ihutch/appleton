This code interactively plots the Appleton-Hartree dispersion relation
governing waves in a cold plasma. Its purpose is primarily
pedagogical.  It is controlled from the keyboard, for which help can
be obtained when it is running by hitting the key 'h'.

To create an executable program on linux or Windows, a fortran
compiler must be available.  The makefile is set up for the mingw
cross-compiler to produce a Windows executable coldplas.exe, or for
gfortran to produce a linux executable coldplaslog. Edit the makefile
to replace MINGW or FORTRAN with references to your compiler. But it
must be compatible with the gfortran/gcc compilers that produced the
required and included graphics libaries libaccisX.a or libaccisWin.a.

Alternatively just download and run the Windows executable:
coldplas.exe. It runs perfectly in native Windows and the 'wine'
system.

On a Mac, to run coldplas.exe you can either install Wine free (se e.g
https://www.davidbaumgold.com/tutorials/wine-mac/) or WineBottler
(http://winebottler.kronenberg.org/) or use a packaged app for running
Windows programs such as WinOnX
(https://apps.apple.com/us/app/winonx-64/) or Crossover
(https://www.codeweavers.com/products).



No guarantee of fitness for any purpose whatever is given, and anyone
using the program does so at their own risk.
Copyright (C) Ian H Hutchinson 2019.
