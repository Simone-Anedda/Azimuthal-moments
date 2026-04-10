# CXX = g++
# CXXFLAGS = -std=c++14

.DEFAULT_GOAL := sources

OBJDIR = ../obj/
MODULES := FragFunct COL photon_flux PhysicsCalculator
SOURCES := $(MODULES).cpp
OBJS := $(OBJDIR)/$(MODULES).o


ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)

.PHONY: sources sources_base sources_full

$(OBJDIR):
	@mkdir -p $(OBJDIR)

# CXXFLAGS += $(ROOTCFLAGS)

# ROOTFLAGS = -m64 -I${HOME}/research/root-6.18.04/obj/include
# ROOTLIBS = -L/home/charlie/research/root-6.18.04/obj/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic

$(OBJDIR)FCN.o: FCN.h FCN.cpp | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) FCN.cpp -c
	@mv FCN.o $(OBJDIR)
	@echo "................FCN done.\n"
	
$(OBJDIR)FragFunct.o: FragFunct.h FragFunct.cpp | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) FragFunct.cpp -c -L/st100-gr4/flore/libs/ -lff -lgfortran
	@mv FragFunct.o $(OBJDIR)
	@echo "................FragFunct done.\n"
	
$(OBJDIR)photon_flux.o: photon_flux.h photon_flux.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) photon_flux.cpp -c
	@mv photon_flux.o $(OBJDIR)
	@echo "................photon_flux done.\n"

$(OBJDIR)PhysicsCalculator.o: PhysicsCalculator.h PhysicsCalculator.cpp photon_flux.h | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) PhysicsCalculator.cpp -c
	@mv PhysicsCalculator.o $(OBJDIR)
	@echo "................PhysicsCalculator done.\n"

$(OBJDIR)COL.o: COL.h COL.cpp | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) COL.cpp -c
	@mv COL.o $(OBJDIR)
	@echo "................COL done.\n"
		
	

sources:
	$(MAKE) $(OBJDIR)FragFunct.o
	$(MAKE) $(OBJDIR)COL.o
	$(MAKE) $(OBJDIR)photon_flux.o
	$(MAKE) $(OBJDIR)PhysicsCalculator.o
