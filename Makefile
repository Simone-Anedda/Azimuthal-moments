REQUESTED_CXX := $(strip $(CXX))
REQUESTED_FF := $(strip $(FF))

CXX := $(or $(strip $(shell [ -n "$(REQUESTED_CXX)" ] && command -v "$(REQUESTED_CXX)" 2>/dev/null)), \
	$(strip $(shell command -v g++-14 2>/dev/null)), \
	$(strip $(shell command -v g++ 2>/dev/null)), \
	$(strip $(shell command -v clang++ 2>/dev/null)), \
	g++)
FF := $(or $(strip $(shell [ -n "$(REQUESTED_FF)" ] && command -v "$(REQUESTED_FF)" 2>/dev/null)), \
	$(strip $(shell command -v gfortran 2>/dev/null)), \
	gfortran)

.DEFAULT_GOAL := predictions_all
 
#CXXFLAGS = -w -std=gnu++0x -O2 -g
CXXFLAGS = -w -O2 -std=c++20 -g 
CPPFLAGS =
LDFLAGS =
FFOPT = -w -ffixed-line-length-0 -std=legacy

HOMEDIR = ~/research/Azimuthal_moments
SRCDIR = src
OBJDIR = obj
#OBJDIR = src
DSSDIR = /hpe-gr4/flore/libs
# TRANSVDIR = transversity
EXEDIR = exe
INPUTDIR = input

MAINEXE = Azimuthal_moments
PREDEXE = Azimuthal_moments_predictions
PREDEXE2 = Azimuthal_moments_predictions_2
PREDEXE2_OLD = Azimuthal_moments_predictions_2_old

MAIN = $(SRCDIR)/Azimuthal_moments.cpp
PRED_MAIN = $(SRCDIR)/Azimuthal_moments_predictions.cpp
PRED_MAIN2 = $(SRCDIR)/Azimuthal_moments_predictions_2.cpp
PRED_MAIN2_OLD = $(SRCDIR)/Azimuthal_moments_predictions_2_old.cpp

INPUT = Azimuthal_moments.input
PREDINPUT = Azimuthal_moments_predictions.input

# Root is the same for both systems...
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
LHAPDF_CPPFLAGS = $(shell lhapdf-config --cppflags 2>/dev/null || pkg-config --cflags lhapdf 2>/dev/null)
LHAPDF_LDFLAGS = $(shell lhapdf-config --ldflags 2>/dev/null || pkg-config --libs lhapdf 2>/dev/null)

CPPFLAGS += $(ROOTCFLAGS) $(LHAPDF_CPPFLAGS)

# HOPPETDIR1 = /st100-gr4/codes/hoppet-1.1.5-modified_cf-build/lib
# HOPPETDIR2 = /st100-gr4/codes/hoppet-1.1.5-modified_cff-build/lib
# LHAPDFDIR = /u/lib
# CUBADIR = /st100-gr4/codes/Cuba-4.2.2-build/lib

# CXXFLAGS += $(ROOTCFLAGS)


# DSSOBJ = $(DSSDIR)/dlib.o $(DSSDIR)/fDSS17LO_separated.o $(DSSDIR)/fDSS17NLO_separated.o \
	 $(DSSDIR)/fDSS_HESSIAN_modified.o $(DSSDIR)/fDSS_modified_pk_separated.o \
	 $(DSSDIR)/grille_had_charged.o $(DSSDIR)/kkp.o $(DSSDIR)/locate.o $(DSSDIR)/pkhff.o \
	 $(DSSDIR)/polin2.o $(DSSDIR)/polint.o
DSSOBJ = /usr/local/lib/libff.a
	 

# GRVOBJ = $(GRVDIR)/grv98.o
# MINUIT2OBJ = $(OBJDIR)/FCN.o
# HESSEOBJ = $(OBJDIR)/Hessian_engine.o
# MINUIT2OBJ = $(OBJDIR)/FCN_dev.o
# DATAOBJ = $(OBJDIR)/getData.o $(OBJDIR)/CovMat.o $(OBJDIR)/ParVec.o
# DATAOBJ = $(OBJDIR)/Data.o $(OBJDIR)/CovMat.o $(OBJDIR)/ParVec.o
# DISTROBJ = $(OBJDIR)/CollPDF.o 
DISTROBJ = $(OBJDIR)/FragFunct.o $(OBJDIR)/COL.o $(OBJDIR)/photon_flux.o $(OBJDIR)/PhysicsCalculator.o
#TRANSVOBJ="transversity/transv_pdf.o transversity/DSSV.o transversity/grv98.o $obj/TRANSV_x.o"
# TRANSVOBJ = $(OBJDIR)/TRANSV.o
HOPPETOBJ="/usr/local/lib/libhoppet_v1.a /usr/local/lib/libhoppet_v1_collins.a /usr/local/lib/libhoppet_v1_collins2.a"
# GSLOBJ = -I${HOME}/gsl-2.6-build/include ${HOME}/gsl-2.6-build/lib/libgsl.so ${HOME}/gsl-2.6-build/lib/libgslcblas.so
GSLOBJ = -I /usr/local/include /usr/local/lib/libgsl.so /usr/local/lib/libgslcblas.so
CUBAOBJ = /usr/local/lib/libcuba.a
# FBTOBJ = $(OBJDIR)/FBT_gsl.o
# FBTOBJ =  $(OBJDIR)/FBT_boost.o

# OBJS   = $(MINUIT2OBJ) $(DATAOBJ)
OBJS = $(DISTROBJ)
OLD_OBJS = $(BASEOBJ)
#$(FBTOBJ)
          
# MCOBJS    = $(FITOBJS)  $(OBJDIR)/MC_engine.o

# BANDSOBJS = $(MINUIT2OBJ) $(OBJDIR)/Data.o $(DISTROBJ) $(TRANSVOBJ) $(COLOBJ) 
#$(AUXFUNCSOBJ)

# BANDSHESSEOBJS = $(FITOBJS)

# TENSORCHARGEOBJS = $(OBJDIR)/CollPDF.o $(OBJDIR)/TRANSV.o 

# LIBS = $(ROOTLIBS) -ltbb -lMinuit2 -lgomp -lff -lgrv -lhoppet_v1 -lgfortran -lLHAPDF -lm   
# LIBS = -lMinuit2 -lgomp -lff -lgrv -ltransv -lcuba -lhoppet_v1 -lhoppet_v1_collins -lhoppet_v1_collins2 -lgfortran -lLHAPDF -lm # -lceres   
# LIBS = -lMinuit2 -fopenmp -lstdc++fs -lboost_iostreams -lboost_system -lboost_filesystem -lff -lgrv -ltransv -lcuba -lhoppet_v1 -lhoppet_v1_collins -lhoppet_v1_collins2 -lgfortran -lLHAPDF -lm 
# LIBS = -lLHAPDF -L/st100-gr4/flore/libs -lff -lcuba -L$(HOPPETDIR1) -lhoppet_v1_collins -L$(HOPPETDIR2) -lhoppet_v1_collins2 -lgfortran  -lm
FF_LIBDIR ?= /st100-gr4/flore/libs
FF_LIBS ?= -lff
CUBA_LIBS ?= -lcuba
HOPPET_LIBS ?= -lhoppet_v1_collins -lhoppet_v1_collins2
FORTRAN_LIBS ?= -lgfortran
MATH_LIBS ?= -lm

LIBS = $(LHAPDF_LDFLAGS) -L$(FF_LIBDIR) $(FF_LIBS) $(CUBA_LIBS) $(HOPPET_LIBS) $(FORTRAN_LIBS) $(MATH_LIBS)
# BANDSLIBS = -lMinuit2 -fopenmp -lstdc++fs -lboost_iostreams -lboost_system -lboost_filesystem -lff -lgrv -ltransv -lcuba -lhoppet_v1 -lhoppet_v1_collins -lhoppet_v1_collins2 -lgfortran -lLHAPDF -lgsl -lgslcblas -lm 
# TENSORCHARGELIBS = -lstdc++fs -lboost_iostreams -lboost_system -lboost_filesystem -lgrv -ltransv -lcuba -lhoppet_v1 -lgfortran -lLHAPDF -lm
# -lceres   

.PHONY: main predictions predictions2 predictions2_old predictions_all \
	sources sources_old all clean clean_all

predictions_all: predictions predictions2

main: sources | $(EXEDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(MAIN) \
	$(OBJS)  $(LIBS) -lc -o $(MAINEXE)
	@if [ -f $(MAINEXE) ]; then \
		mv $(MAINEXE) $(EXEDIR); echo "...................Azimuthal_moments done.\n"; \
	else \
		echo "could not create executable for '$(MAINEXE)'";\
	fi
	

predictions: sources | $(EXEDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(PRED_MAIN) \
	$(OBJS)  $(LIBS) -lc -o $(PREDEXE)
	@if [ -f $(PREDEXE) ]; then \
		mv $(PREDEXE) $(EXEDIR); echo "...................Azimuthal_moments_predictions done.\n"; \
	else \
		echo "could not create executable for '$(PREDEXE)'";\
	fi

predictions2: sources | $(EXEDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(PRED_MAIN2) \
	$(OBJS)  $(LIBS) -lc -o $(PREDEXE2)
	@if [ -f $(PREDEXE2) ]; then \
		mv $(PREDEXE2) $(EXEDIR); echo "...................Azimuthal_moments_predictions_2 done.\n"; \
	else \
		echo "could not create executable for '$(PREDEXE2)'";\
	fi

predictions2_old: sources_old | $(EXEDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(PRED_MAIN2_OLD) \
	$(OLD_OBJS) $(LIBS) -lc -o $(PREDEXE2_OLD)
	@if [ -f $(PREDEXE2_OLD) ]; then \
		mv $(PREDEXE2_OLD) $(EXEDIR); echo "...................Azimuthal_moments_predictions_2_old done.\n"; \
	else \
		echo "could not create executable for '$(PREDEXE2_OLD)'";\
	fi

$(EXEDIR):
	@mkdir -p $(EXEDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)


sources: | $(OBJDIR)
	$(MAKE) --directory=src sources CXX="$(CXX)" CPPFLAGS="$(CPPFLAGS)" CXXFLAGS="$(CXXFLAGS)"

sources_old: | $(OBJDIR)
	$(MAKE) --directory=src sources_base CXX="$(CXX)" CPPFLAGS="$(CPPFLAGS)" CXXFLAGS="$(CXXFLAGS)"
	
all: predictions_all

# run_fit:
# 	$(MAKE) fit_main
# 	@time ./exe/fit_main $(cat $(INPUTDIR)/$(FITINPUT))
	
	
clean:
	@rm -f ./obj/*.o
	
clean_all:
	@rm -f ./obj/*.o
	@rm -f ./exe/*
# test.o:	test.cpp $(OBJECTS)
# 	$(CXX) $(OPT) $(CFLAGS) -o $@ test.cpp
# 	@echo "..................done Test."
# 
# 
# clean:
# 	rm *.o
# 	rm *.exe
# 	$(MAKE) --directory=Particles clean
# 	$(MAKE) --directory=Processes clean
# 
# all:
# 	$(MAKE) --directory=Particles
# 	$(MAKE) --directory=Processes
# 	make test.exe
