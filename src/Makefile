# CXX = g++
# CXXFLAGS = -std=c++14

OBJDIR = ../obj/
MODULES := FragFunct COL photon_flux
SOURCES := $(MODULES).cpp
OBJS := $(OBJDIR)/$(MODULES).o


ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)

# CXXFLAGS += $(ROOTCFLAGS)

# ROOTFLAGS = -m64 -I${HOME}/research/root-6.18.04/obj/include
# ROOTLIBS = -L/home/charlie/research/root-6.18.04/obj/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic

$(OBJDIR)FCN.o: FCN.h FCN.cpp
	$(CXX) $(CXXFLAGS) FCN.cpp -c
	@mv FCN.o $(OBJDIR)
	@echo "................FCN done.\n"
	
	
	
$(OBJDIR)FragFunct.o: FragFunct.h FragFunct.cpp
	$(CXX) $(CXXFLAGS) FragFunct.cpp -c -L/st100-gr4/flore/libs/ -lff -lgfortran
	@mv FragFunct.o $(OBJDIR)
	@echo "................FragFunct done.\n"
	
$(OBJDIR)photon_flux.o: photon_flux.h photon_flux.cpp
	$(CXX) $(CXXFLAGS) photon_flux.cpp -c
	@mv photon_flux.o $(OBJDIR)
	@echo "................photon_flux done.\n"

$(OBJDIR)COL.o: COL.h COL.cpp
	$(CXX) $(CXXFLAGS) COL.cpp -c
	@mv COL.o $(OBJDIR)
	@echo "................COL done.\n"
		
	

sources:
	$(MAKE) $(OBJDIR)FragFunct.o
	$(MAKE) $(OBJDIR)COL.o
	$(MAKE) $(OBJDIR)photon_flux.o
