ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
EVENTLIB      = $(EVENTO)
OutPutOpt     = -o
 
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
 
# Linux with egcs
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

# Solaris
#CXX           = CC
#CXXFLAGS      = -O -KPIC
#LD            = CC
#LDFLAGS       = -g
#SOFLAGS       = -G
 
 
CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS)
 
#------------------------------------------------------------------------------
JETFINDO     = JetFinder.$(ObjSuf) JetFindDict.$(ObjSuf)
JETFINDS     = JetFinder.$(SrcSuf)
JETFINDSO    = libJetFind.$(DllSuf) 

EVENTSHAPEO  = EventShape.$(ObjSuf) EventShapeDict.$(ObjSuf)
EVENTSHAPES  = EventShape.$(SrcSuf)
EVENTSHAPESO = libEventShape.$(DllSuf) 


OBJS          = $(JETFINDO)  $(EVENTSHAPEO)

PROGRAMS      = $(JETFINDSO) $(EVENTSHAPESO)

all: $(PROGRAMS)


$(JETFINDSO): $(JETFINDO) 
	      $(LD) $(SOFLAGS) $(LDFLAGS) $(JETFINDO) $(LIBS) \
              $(OutPutOpt) $(JETFINDSO)

$(EVENTSHAPESO): $(EVENTSHAPEO) 
	         $(LD) $(SOFLAGS) $(LDFLAGS) $(EVENTSHAPEO) $(LIBS) \
                 $(OutPutOpt) $(EVENTSHAPESO)

clean:
		@rm -f $(OBJS) *Dict.* core

.SUFFIXES: .$(SrcSuf)

###

JetFindDict.$(SrcSuf):JetFinder.h
	@echo "Generating dictionary ..."
	@$(ROOTSYS)/bin/rootcint -f $@ -c $+

EventShapeDict.$(SrcSuf):EventShape.h
	@echo "Generating dictionary ..."
	@$(ROOTSYS)/bin/rootcint -f $@ -c $+

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
