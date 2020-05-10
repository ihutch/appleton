This code interactively plots the Appleton-Hartree dispersion relation
governing waves in a cold plasma. Its purpose is primarily
pedagogical.  It is controlled from the keyboard, for which help can
be obtained when it is running by hitting the key 'h'.

On MSWindows, download just the complete executable: coldplas.exe. It
runs perfectly in native Windows and under the 'wine' system.

On MacOS, you can run coldplas.exe if you install wine by
$ brew cask install wine-stable
(Or see e.g https://www.davidbaumgold.com/tutorials/wine-mac/)

To create an executable program from source on MacOS or linux type
$ make
it will git clone and compile the required graphics library.
But you need gfortran and the X11 libraries accessible.
[On MacOS $ brew install gcc; brew cask install xquartz ]
It is usually best run from the command line by: $ ./coldplas
The native compiled version is considerably more responsive than
running under wine. 

For the adventurous only: the makefile with 'make windows' is set up
to use the mingw-32 cross-compiler to produce a Windows executable
coldplas.exe.  Edit the makefile to replace MINGW or FORTRAN with
references to your compiler. But it must be compatible with the
gfortran/gcc compilers that produced the required and included
graphics libaries libaccisX.a or libaccisWin.a.

No guarantee of fitness for any purpose whatever is given, and anyone
using the program does so at their own risk.
Copyright (C) Ian H Hutchinson 2017-2020.
