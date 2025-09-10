CXX = g++
FF = gfortran
 
#CXXFLAGS = -w -std=gnu++0x -O2 -g
CXXFLAGS = -w -O2 -std=c++14 -g 
FFOPT = -w -ffixed-line-length-0 -std=legacy

HOMEDIR = ~/research/Azimuthal_moments
SRCDIR = src
OBJDIR = obj
#OBJDIR = src
DSSDIR = dss
# TRANSVDIR = transversity
EXEDIR = exe
INPUTDIR = input

MAINEXE = Azimuthal_moments

MAIN = $(SRCDIR)/Azimuthal_moments.cpp

INPUT = Azimuthal_moments.input

# Root is the same for both systems...
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)

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
DISTROBJ = $(OBJDIR)/FragFunct.o $(OBJDIR)/COL.o
#TRANSVOBJ="transversity/transv_pdf.o transversity/DSSV.o transversity/grv98.o $obj/TRANSV_x.o"
# TRANSVOBJ = $(OBJDIR)/TRANSV.o
HOPPETOBJ="/usr/local/lib/libhoppet_v1.a /usr/local/lib/libhoppet_v1_collins.a /usr/local/lib/libhoppet_v1_collins2.a"
# GSLOBJ = -I${HOME}/gsl-2.6-build/include ${HOME}/gsl-2.6-build/lib/libgsl.so ${HOME}/gsl-2.6-build/lib/libgslcblas.so
GSLOBJ = -I /usr/local/include /usr/local/lib/libgsl.so /usr/local/lib/libgslcblas.so
CUBAOBJ = /usr/local/lib/libcuba.a
# FBTOBJ = $(OBJDIR)/FBT_gsl.o
# FBTOBJ =  $(OBJDIR)/FBT_boost.o

# OBJS   = $(MINUIT2OBJ) $(DATAOBJ) 
OBJS = $(DISTROBJ) $(COLOBJ) $(TRANSVOBJ)
#$(FBTOBJ)
          
# MCOBJS    = $(FITOBJS)  $(OBJDIR)/MC_engine.o

# BANDSOBJS = $(MINUIT2OBJ) $(OBJDIR)/Data.o $(DISTROBJ) $(TRANSVOBJ) $(COLOBJ) 
#$(AUXFUNCSOBJ)

# BANDSHESSEOBJS = $(FITOBJS)

# TENSORCHARGEOBJS = $(OBJDIR)/CollPDF.o $(OBJDIR)/TRANSV.o 

# LIBS = $(ROOTLIBS) -ltbb -lMinuit2 -lgomp -lff -lgrv -lhoppet_v1 -lgfortran -lLHAPDF -lm   
# LIBS = -lMinuit2 -lgomp -lff -lgrv -ltransv -lcuba -lhoppet_v1 -lhoppet_v1_collins -lhoppet_v1_collins2 -lgfortran -lLHAPDF -lm # -lceres   
# LIBS = -lMinuit2 -fopenmp -lstdc++fs -lboost_iostreams -lboost_system -lboost_filesystem -lff -lgrv -ltransv -lcuba -lhoppet_v1 -lhoppet_v1_collins -lhoppet_v1_collins2 -lgfortran -lLHAPDF -lm 
LIBS = -lff -lcuba -lhoppet_v1 -lhoppet_v1_collins -lhoppet_v1_collins2 -lgfortran -lLHAPDF -lm 
# BANDSLIBS = -lMinuit2 -fopenmp -lstdc++fs -lboost_iostreams -lboost_system -lboost_filesystem -lff -lgrv -ltransv -lcuba -lhoppet_v1 -lhoppet_v1_collins -lhoppet_v1_collins2 -lgfortran -lLHAPDF -lgsl -lgslcblas -lm 
# TENSORCHARGELIBS = -lstdc++fs -lboost_iostreams -lboost_system -lboost_filesystem -lgrv -ltransv -lcuba -lhoppet_v1 -lgfortran -lLHAPDF -lm
# -lceres   

main: $(OBJS) $(EXEDIR)/$(MAINEXE)
	@make sources
	$(CXX) $(CXXFLAGS) $(MAIN) \
	$(OBJS)  $(LIBS) -lc -o $(MAINEXE)
	@if [ -f $(MAINEXE) ]; then \
		mv $(MAINEXE) $(EXEDIR); echo "...................Azimuthal_moments done.\n"; \
	else \
		echo "could not create executable for '$(MAINEXE)'";\
	fi
	

sources:
	$(MAKE) --directory=src sources
	
all:
	$(MAKE) --directory=src sources

# run_fit:
# 	$(MAKE) fit_main
# 	@time ./exe/fit_main $(cat $(INPUTDIR)/$(FITINPUT))
	
	
clean:
	@rm ./obj/*.o
	
clean_all:
	@rm ./obj/*.o
	@rm ./exe/*
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
