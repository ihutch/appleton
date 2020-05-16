This code interactively plots the Appleton-Hartree dispersion relation
governing waves in a cold plasma. 

Its purpose is primarily pedagogical, providing rapid evolution of the
display, which animates the variation of parameters. Playing with it helps
to develop an intuitive understanding of the complex physics of plasma
waves. It is controlled from the keyboard. Help can be obtained when
it is running by hitting the key 'h'.
_______________________________________________________________________

On MSWindows, download and run just the complete executable:
coldplas.exe.  That executable also runs perfectly under the 'wine'
system. On MacOS, if you have wine installed then 'wine coldplas.exe'
will work fine.

To create an executable program from source, probably the best option
on MacOS and certainly on linux, at the commandline type 'make' and it
will git clone and compile the required graphics library and compile
coldplas.  The command './coldplas' will run the program.

You need git, gfortran and the X11 libraries accessible for compilation.
On MacOS, fulfilling those needs involves something like either:

/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
In Finder double click on Homebrew to install, then in Terminal
brew install gcc
brew cask install xquartz
make

or

From https://www.macports.org/install.php
download the install file of MacPorts appropriate for your macOS version.
In Finder double click on it to install it, then in Terminal
sudo port install xorg-libX11
sudo port install gcc5
sudo port select --set gcc mp-gcc5

Don't try to 'make coldplas.exe' using a cross-compiler unless you
have one, really know what you are doing, and can hack the makefile.

No guarantee of fitness for any purpose whatever is given, and anyone
using the program does so at their own risk.
Copyright (C) Ian H Hutchinson 2017-2020.
