# This is a sample automatic makefile to compile ROOT code with a typical C++
# compiler. It requires that $ROOTSYS/bin be added to $PATH and that 
# $ROOTSYS/lib be added to $LD_LIBRARY_PATH.
#
# Author : Pier-Olivier DeViveiros
# May 22nd 2007
# Thanks to Francois-Xavier Gentit for his help in creating this.

# -----------------------------------------------------------------

# Define the name of the project
# Make sure that this is the same name as the source file

NAME        = lsRoot
SUTILS        = stringUtils
LSUTILS        = lsUtils
LSCONFIG	= lsConfig
LSPRINT		= lsPrint
LSDRAW		= lsDraw

# -----------------------------------------------------------------

# Define Suffixes for diferent filetypes
ObjSuf        = o
SrcSuf        = cpp
ExeSuf        = 
OutPutOpt     = -o # keep whitespace after "-o"

# -----------------------------------------------------------------

# Grab the flags required by ROOT when compiling using the
# root-config shell commands
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
# If root-config --cflags contains -m32 then add it to the --libs flags!
ROOTLIBS32	=  $(shell root-config --cflags | grep -i '\(m32\|m64\)') 
ROOTLIBS	+= $(ROOTLIBS32)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTINCS     = $(shell root-config --incdir)

# -----------------------------------------------------------------

# Compiler Information
CXX           = g++
CXXFLAGS      = -g -Wall -fPIC   # Debug
#CXXFLAGS     = -O -Wall -fPIC   # Optimized

# Linker Information
LD            = g++
LDFLAGS       = -g      # Debug
#LDFLAGS       = -O      # Optimized

SOFLAGS       = -c

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#------------------------------------------------------------------------------

SRCS          = $(NAME).$(SrcSuf) $(SUTILS).$(SrcSuf) $(LSUTILS).$(SrcSuf) $(LSCONFIG).$(SrcSuf) $(LSPRINT).$(SrcSuf) $(LSDRAW).$(SrcSuf)
OBJS          = $(NAME).$(ObjSuf) $(SUTILS).$(ObjSuf) $(LSUTILS).$(ObjSuf) $(LSCONFIG).$(ObjSuf) $(LSPRINT).$(ObjSuf) $(LSDRAW).$(ObjSuf)

PROGRAM        = $(NAME)$(ExeSuf)

# -----------------------------------------------------------------

all: $(PROGRAM)

.$(SrcSuf).$(ObjSuf):

#	The real Compile Command
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -I$(ROOTINCS) $^ $(OutPutOpt)$@

$(PROGRAM): $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$(PROGRAM)
	@echo "$(PROGRAM) done!"
	@echo "cp .lsRoot_colors ~/ for color formatting in $(PROGRAM)"

clean:   
	rm *.$(ObjSuf)
#	rm $(NAME)$(ExeSuf)
# 	rm $(NAME).$(ObjSuf)
# 	rm $(SUTILS).$(ObjSuf)
# 	rm $(LSUTILS).$(ObjSuf)
# 	rm $(LSCONFIG).$(ObjSuf)
# 	rm $(LSPRINT).$(ObjSuf)
# 	rm $(NAME)$(ExeSuf)

cleanall:

	rm *.$(ObjSuf)
	rm $(NAME)$(ExeSuf)

