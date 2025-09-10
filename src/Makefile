# CXX = g++
# CXXFLAGS = -std=c++14

OBJDIR = ../obj/
MODULES := FragFunct COL 
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
	$(CXX) $(CXXFLAGS) FragFunct.cpp -c -lff -lgfortran
	@mv FragFunct.o $(OBJDIR)
	@echo "................FragFunct done.\n"
	
# $(OBJDIR)COL_z.o: COL_z.h COL_z.cpp
# 	$(CXX) $(CXXFLAGS) COL_z.cpp -c
# 	@mv COL_z.o $(OBJDIR)
# 	@echo "................COL_z done.\n"
# 	
# $(OBJDIR)TRANSV_x.o: TRANSV_x.h TRANSV_x.cpp
# 	$(CXX) $(CXXFLAGS) TRANSV_x.cpp -c
# 	@mv TRANSV_x.o $(OBJDIR)
# 	@echo "................TRANSV_x done.\n"

$(OBJDIR)COL.o: COL.h COL.cpp
	$(CXX) $(CXXFLAGS) COL.cpp -c
	@mv COL.o $(OBJDIR)
	@echo "................COL done.\n"
		
	

sources:
# 	$(MAKE) $(OBJDIR)Data.o
# 	$(MAKE) $(OBJDIR)FCN.o
# 	$(MAKE) $(OBJDIR)Hessian_engine.o
# 	$(MAKE) $(OBJDIR)MC_engine.o
# 	$(MAKE) $(OBJDIR)ParVec.o
# 	$(MAKE) $(OBJDIR)CovMat.o
# 	$(MAKE) $(OBJDIR)CollPDF.o
	$(MAKE) $(OBJDIR)FragFunct.o
	$(MAKE) $(OBJDIR)COL.o
# 	$(MAKE) $(OBJDIR)TRANSV.o
# 	$(MAKE) $(OBJDIR)auxiliary_funcs.o
