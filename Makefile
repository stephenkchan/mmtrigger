Target  = MMtriggerSim
Objects = newvecdic.o MMT_struct.o MMT_Loader.o MMT_Finder.o MMT_Fitter.o MM_Trigger.o #AutoDict_vector_vector_float___.o AutoDict_vector_vector_int___.o

ROOTFLAGS = $(shell root-config --cflags) -fPIC
ROOTLIBS = $(shell root-config --libs) -lTMVA -lMinuit 
#lxplus pattern:
#INCBOOST = -I/cvmfs/sft.cern.ch/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc48-opt/include/boost-1_53/boost 
#LIBBOOST = -L/cvmfs/sft.cern.ch/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc48-opt/lib
#herophysics pattern:
BOOST_STEM = /cvmfs/atlas.cern.ch/repo/sw/software/x86_64-slc6-gcc48-opt/19.1.1/sw/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc48-opt/
#/cvmfs/atlas.cern.ch/repo/sw/software/17.7.90/sw/lcg/external/Boost/1.48.0_python2.6/i686-slc5-gcc43-opt/
INCBOOST = -I$(BOOST_STEM)include/boost-1_53/
LIBBOOST = -L$(BOOST_STEM)lib/ 

INCHERE = -I./
LIBHERE = -I./ -L./ 

all:$(Target)

MMtriggerSim: MMtriggerSim.cxx $(Objects)
	g++ $(INCHERE) $(INCBOOST) -o $@ $^ $(ROOTFLAGS) $(ROOTLIBS) $(LIBBOOST) $(LIBHERE) 

#vectorDict.o: vectorDict.cxx
#	g++ -c vectorDict.cxx $(ROOTFLAGS) $(ROOTLIBS)

newvecdic.o: newvecdic.cxx
	g++ -c newvecdic.cxx $(ROOTFLAGS) $(ROOTLIBS)

MMT_struct.o: MMT_struct.cxx
	g++ -c MMT_struct.cxx $(INCHERE) $(INCBOOST) $(ROOTFLAGS) $(ROOTLIBS) $(LIBBOOST) $(LIBHERE) 

MMT_Loader.o: MMT_Loader.cxx 
	g++ $(INCHERE) $(INCBOOST) -c MMT_Loader.cxx $(ROOTFLAGS) $(ROOTLIBS) $(LIBBOOST) $(LIBHERE) 

MMT_Fitter.o: MMT_Fitter.cxx MMT_struct.o
	g++ $(INCHERE) $(INCBOOST) -c MMT_Fitter.cxx $(ROOTFLAGS) $(ROOTLIBS) $(LIBBOOST) $(LIBHERE) 

MMT_Finder.o: MMT_Finder.cxx MMT_struct.o MMT_Fitter.o
	g++ $(INCHERE) $(INCBOOST) -c MMT_Finder.cxx $(ROOTFLAGS) $(ROOTLIBS) $(LIBBOOST) $(LIBHERE) 

MM_Trigger.o: MM_Trigger.cxx MMT_struct.o MMT_Loader.o MMT_Fitter.o MMT_Finder.o
	g++ $(INCHERE) $(INCBOOST) -c MM_Trigger.cxx $(ROOTFLAGS) $(ROOTLIBS) $(LIBBOOST) $(LIBHERE) 

clean: 
	rm *.o *~

