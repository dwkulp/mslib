############################################################################################
#
#  Users can set some enviromental variables to inform MSL about the presence of external 
#  libraries and to use debug level output
#
#   GSL (Gnu Scientific Library)
#     The GSL library is highly reccommended. It is used for alignment and minimization.
#     Required library files are libgsl.a libgslcblas.a (package libgsl0-dev in Ubuntu)
#      
#     If installed, set the environmental variable $MSL_GSL to "T" (default), else to "F"
#     How to set $MSL_GSL in your .cshrc: add
#          setenv MSL_GSL T
#     How to set $MSL_GSL in your .bash: add
#          export MSL_GSL=T
#
#   BOOST C++ libraries
#     The BOOST library is optional.  
#     Required libraries: libboost_serialization.a libboost_regex.a
#     If installed, set the environmental variable $MSL_BOOST to "T", else to "F" (default)
#
#   GLPK (Gnu Linear Programming Kit)
#     The GLPK library is normally optional and ***currently disabled***.  
#     This option is useful only if one is using the LinearProgrammingOptimization object
#     and the programs that support it
#     Required libraries: libglpk.a
#     If installed, set the environmental variable $MSL_GLPK to "T", else to "F" (default)
#
#   R libraries
#     The R library is optional, allowing interfacing with the statitical package R  
#     Required libraries: ?
#     If installed, set the environmental variable $MSL_R to "T", else to "F" (default)
#
#   Debug mode
#     Set $MSL_DEBUG to "T" to compile in "debug" mode, or else to "F" (default)
# 
############################################################################################

# Define compiler command
#CC  = g++ 
CCOPTIM = g++ -Wall -Wno-parentheses -Wno-sign-compare -O3 -msse3 -mfpmath=sse -funroll-loops 
CCDEBUG = g++ -Wall -Wno-parentheses -Wno-sign-compare -msse3 -mfpmath=sse -funroll-loops -g


GSLDEFAULT = F
GSLOLDDEFAULT = F
GLPKDEFAULT = F
BOOSTDEFAULT = F
#OPENMPDEFAULT = T
ARCH32BITDEFAULT = F
FFTWDEFAULT = F
RDEFAULT = F
MACOSDEFAULT = F
DEBUGDEFAULT = F
TESTINGDEFAULT = F
MSLOUT_DEBUG_OFFDEFAULT = F
STATICDEFAULT = T
EXTERNAL_LIB_DIR_DEFAULT=/usr/lib
EXTERNAL_INCLUDE_DIR_DEFAULT=/usr/include

VPATH = src


SOURCE  = ALNReader Atom Atom3DGrid AtomAngleRelationship AtomContainer AtomDihedralRelationship AtomDistanceRelationship \
          AtomGeometricRelationship AtomGroup AtomSelection AtomPointerVector CartesianGeometry \
          BaselineEnergyBuilder BaselineInteraction BBQTable BBQTableReader BBQTableWriter CartesianPoint\
          Chain CharmmAngleInteraction CharmmBondInteraction CharmmDihedralInteraction \
          CharmmElectrostaticInteraction CharmmEnergy CharmmIMM1Interaction CharmmIMM1RefInteraction CharmmImproperInteraction CharmmParameterReader CharmmEEF1ParameterReader \
          CharmmSystemBuilder CharmmTopologyReader CharmmTopologyResidue CharmmUreyBradleyInteraction \
          CharmmVdwInteraction CharmmEEF1Interaction CharmmEEF1RefInteraction ChiStatistics CoiledCoils CrystalLattice DeadEndElimination EnergySet EnergeticAnalysis Enumerator EnvironmentDatabase \
          EnvironmentDescriptor File FormatConverter FourBodyInteraction Frame FuseChains Helanal HydrogenBondBuilder IcEntry IcTable Interaction \
          InterfaceResidueDescriptor Line LogicalParser MIDReader Matrix Minimizer MoleculeInterfaceDatabase \
          MslOut MslTools OptionParser CRDFormat PDBFormat PDBReader PDBWriter PDBTopology CRDReader CRDWriter PolymerSequence PSFReader \
          Position PotentialTable Predicate PrincipleComponentAnalysis PyMolVisualization Quaternion Reader Residue ResiduePairTable \
          ResiduePairTableReader ResidueSelection ResidueSubstitutionTable ResidueSubstitutionTableReader RotamerLibrary \
          RotamerLibraryReader SidechainOptimizationManager SelfPairManager SasaAtom SasaCalculator Scwrl4HBondInteraction SphericalPoint SurfaceSphere Symmetry System SystemRotamerLoader TBDReader \
          ThreeBodyInteraction Timer Transforms Tree TwoBodyDistanceDependentPotentialTable OneBodyInteraction TwoBodyInteraction Writer UserDefinedInteraction  UserDefinedEnergy \
          UserDefinedEnergySetBuilder HelixGenerator RotamerLibraryBuilder RotamerLibraryWriter AtomBondBuilder LogicalCondition MonteCarloManager \
	  SelfConsistentMeanField PhiPsiReader PhiPsiStatistics RandomNumberGenerator \
	  BackRub CCD MonteCarloOptimization Quench SpringConstraintInteraction SurfaceAreaAndVolume VectorPair VectorHashing PDBTopologyBuilder SysEnv \
	  FastaReader PSSMCreator PrositeReader PhiPsiWriter ConformationEditor DegreeOfFreedomReader OnTheFlyManager CharmmEnergyCalculator EZpotentialInteraction EZpotentialBuilder \
	  OptimalRMSDCalculator DSSPReader StrideReader CIFWriter CIFReader CIFFormat



HEADER = Hash.h MslExceptions.h Real.h Selectable.h Tree.h release.h 

# Quick test that might or might not work for you
SANDBOX = testAtomGroup testAtomSelection testAtomPointerVector testBBQ testBBQ2 \
          testCharmmTopologyReader testCoiledCoils testEnergySet testEnergeticAnalysis testEnvironmentDatabase \
          testEnvironmentDescriptor testFrame testFormatConverter testGenerateCrystalLattice testIcBuilding testLoopOverResidues \
          testMolecularInterfaceDatabase testMslToolsFunctions testCRDIO testPDBIO testPhiPsi testPolymerSequence testPSFReader \
          testResiduePairTable testResidueSubstitutionTable testSasaCalculator testSymmetry testSystemCopy \
          testSystemIcBuilding testTransforms testHelixGenerator testRotamerLibraryWriter testALNReader \
	  testAtomAndResidueId testAtomBondBuilder testTransformBondAngleDiheEdits testAtomContainer testCharmmEEF1ParameterReader \
	  testResidueSelection testMslOut testMslOut2 testRandomNumberGenerator \
	  testPDBTopology testVectorPair testSharedPointers2 testTokenize testSaveAtomAltCoor testPDBTopologyBuild testSysEnv \
	  testConformationEditor testDeleteBondedAtom testOptimalRMSDCalculator testRosettaScoredPDBReader testClustering testBebl

# These tests need to be compile before a commit can be contributed to the repository
LEAD =    

# These tests need to be passed before a commit can be contributed to the repository
GOLD =    testEZpotential testCharmmEnergies testCharmmBuild

PROGRAMS = getSphericalCoordinates fillInSideChains generateCrystalLattice getDihedrals energyTable analEnergy \
	   getSelection calculateSasa printSequence \
           insertLoopIntoTemplate setConformation coiledCoilBuilder findClashes mutate calculateDistanceOrAngle \
	   repackSideChains backrubPdb renumberResidues getChiRecovery createEnergyTable createEBL \
	   designSideChains generateCoiledCoils trimConformerLibrary pdb2crd glycanTransfer

# PROGRAMS/SANDBOX_THAT_DO_NOT_COMPLILE = testBoost testRInterface  testLinkedPositions testEEF1 testEEF1_2  testAddCharmmIdentity testNonBondedCutoff 
# PROGRAMS/SANDBOX_THAT_COMPILE_BUT_SEGFAULT =  testTree
# PROGRAMS/SANDBOX_WITHOUT_SOURCE_FILES =  testBoostSpriit testBoostSpirit2 testLogicalCondition testDistanceHashing 


# To ever-ride the defaults, set the $MSL_GSL $MSL_GLPK $MSL_BOOST $MSL_DEBUG and $MSL_TESTING environmental variables
ifndef MSL_DEBUG
   MSL_DEBUG = DEBUGDEFAULT
endif
ifndef MSL_TESTING
   MSL_TESTING = TESTINGDEFAULT
endif
ifndef MSL_GSL
   MSL_GSL=${GSLDEFAULT}
endif
ifndef MSL_GSL_OLD
   MSL_GSL_OLD=${GSLOLDDEFAULT}
endif
ifndef MSL_GLPK
   MSL_GLPK=${GLPKDEFAULT}
endif
ifndef MSL_BOOST
   MSL_BOOST=${BOOSTDEFAULT}
endif
#ifndef MSL_OPENMP
#   MSL_OPENMP=${OPENMPDEFAULT}
#endif
ifndef MSL_STATIC
   MSL_STATIC=${STATICDEFAULT}
endif

ifndef ARCH32BIT
   ARCH32BIT=${ARCH32BITDEFAULT}
endif
ifndef FFTW
   FFTW=${FFTWDEFAULT}
endif

ifndef MSL_R
    MSL_R=${RDEFAULT}
endif

ifndef MSL_EXTERNAL_LIB_DIR
   MSL_EXTERNAL_LIB_DIR=${EXTERNAL_LIB_DIR_DEFAULT}
endif

ifndef MSL_EXTERNAL_INCLUDE_DIR
   MSL_EXTERNAL_INCLUDE_DIR=${EXTERNAL_INCLUDE_DIR_DEFAULT}
endif

ifndef MSL_MACOS
    MSL_MACOS=${MACOSDEFAULT}
endif

ifndef MSL_MSLOUT_DEBUG_OFF
   MSL_MSLOUT_DEBUG_OFF=${MSLOUT_DEBUG_OFFDEFAULT}
endif

ifeq ($(MSL_DEBUG),T)
    CC = ${CCDEBUG}
#   FLAGS =   -Wall -Wno-sign-compare -msse3 -mfpmath=sse -funroll-loops -g 
#   LINKFLAGS =
else
    CC= ${CCOPTIM}
#   FLAGS =  -Wall -Wno-sign-compare -O3 -msse3 -mfpmath=sse -funroll-loops
#   LINKFLAGS =
endif


# compilation with alternative testing code
ifeq ($(MSL_TESTING),T)
    FLAGS          += -D__TESTING__
    SOURCE         += 
    SANDBOX        += 
    GOLD           += 
    LEAD           +=
    PROGRAMS       += 
    STATIC_LIBS    += 
endif

# GLPK Libraries
ifeq ($(MSL_GLPK),T)
    FLAGS          += -D__GLPK__
    SOURCE         += LinearProgrammingOptimization
    SANDBOX        += 
    GOLD           += 
    LEAD           += 
    PROGRAMS       += optimizeLP
    STATIC_LIBS    += ${MSL_EXTERNAL_LIB_DIR}/libglpk.a

# DO_NOT_COMPILE = testRotamerOptimization
endif


# GSL Libraries 
ifeq ($(MSL_GSL),T)
    FLAGS          += -D__GSL__
    SOURCE         += GSLMinimizer HelixFusion CoiledCoilFitter Clustering
    SANDBOX        += testDerivatives testCCD testBackRub testSurfaceAreaAndVolume testHelixFusion testMinimization
    GOLD           += testRMSDalignment
    LEAD           +=
    PROGRAMS       += tableEnergies runQuench runKBQuench optimizeMC alignMolecules searchFragmentDatabase getSurroundingResidues minimize 
    STATIC_LIBS    += ${MSL_EXTERNAL_LIB_DIR}/libgsl.a ${MSL_EXTERNAL_LIB_DIR}/libgslcblas.a

#  DO_NOT_COMPILE = testQuench 
endif
# GSL OLD remove some functions in the Minimizer when the GSL version available is pre-v1.14
ifeq ($(MSL_GSL_OLD),T)
    FLAGS          += -D__GSL_OLD__
endif

# BOOST Libraries
ifeq ($(MSL_BOOST),T)
    FLAGS          += -D__BOOST__ -DBOOST_DISABLE_THREADS

    SOURCE         +=  RegEx RandomSeqGenerator RosettaScoredPDBReader 
    ifeq ($(MSL_GSL),T)
        SOURCE         +=  PDBFragments
    endif
    SANDBOX        += testRegEx testRandomSeqGenerator testBoost
    GOLD           +=
    LEAD           +=
    PROGRAMS       +=  createFragmentDatabase 
    ifeq ($(MSL_GSL),T)
	    PROGRAMS       +=  grepSequence
    endif
    STATIC_LIBS    += ${MSL_EXTERNAL_LIB_DIR}/libboost_serialization.a   
    # For MAC I only compiled the non-multithreaded library, sometime I'll figure it out, but we do not use multi-threading so for now this is ok.
    ifeq ($(MSL_MACOS),T)
        STATIC_LIBS    += ${MSL_EXTERNAL_LIB_DIR}/libboost_regex.a
    else
        STATIC_LIBS    += ${MSL_EXTERNAL_LIB_DIR}/libboost_regex.a
    endif
endif

#ifeq ($(MSL_OPENMP),T)
#    FLAGS          += -fopenmp -D__OPENMP__
#endif

ifeq ($(FFTW),T)
    STATIC_LIBS    += ${MSL_EXTERNAL_LIB_DIR}/libfftw3.a
endif

# R Libraries
ifeq ($(MSL_R),T)

#     R_HOME         = /opt/applications/R/2.12.1/gnu/
     FLAGS         += -D__R__ 

ifeq ($(MSL_STATIC),T)     
     STATIC_LIBS   += ${MSL_EXTERNAL_LIB_DIR}/libRcpp.a ${MSL_EXTERNAL_LIB_DIR}/libRInside.a
else

ifeq ($(MSL_MACOS),T)
     LINKFLAGS     += -framework R
else

     R_HOME := 		$(shell R RHOME | grep -v WARNING)
     RCPPFLAGS := 	$(shell $(R_HOME)/bin/R CMD config --cppflags)
     RLDFLAGS := 	$(shell $(R_HOME)/bin/R CMD config --ldflags)
     RBLAS := 		$(shell $(R_HOME)/bin/R CMD config BLAS_LIBS)
     RLAPACK := 	$(shell $(R_HOME)/bin/R CMD config LAPACK_LIBS)

     ## if you need to set an rpath to R itself, also uncomment
     #RRPATH :=		-Wl,-rpath,$(R_HOME)/lib

     ## include headers and libraries for Rcpp interface classes
     RCPPINCL := 		$(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
     RCPPLIBS := 		$(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)


     ## include headers and libraries for RInside embedding classes
     RINSIDEINCL := 		$(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
     RINSIDELIBS := 		$(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

     FLAGS         += $(RCPPFLAGS) $(RCPPINCL) $(RINSIDEINCL) $(shell $(R_HOME)/bin/R CMD config CXXFLAGS | grep -v WARN)
     LINKFLAGS     += $(RLDFLAGS) $(RRPATH) $(RBLAS) $(RLAPACK) $(RCPPLIBS) $(RINSIDELIBS)

endif  # IF MACOS


endif # IF MSL_STATIC

endif # IF MSL_R


ifeq ($(MSL_MSLOUT_DEBUG_OFF),T)
    FLAGS += -D__MSL_MSLOUT_DEBUG_OFF__
endif
# Generic Includes,Flags.  Static compile.  
# NOTE IS THE FOLLOWING STILL NECESSARY?
INCLUDE  = src -Iprograms -I${MSL_EXTERNAL_INCLUDE_DIR} 



# Add a MACOS flag for certain code breaks (see bottom of Tree.h, Selectable.h ... templated classes don't need pre-instantiations?)
ifeq ($(MSL_MACOS),T)
    FLAGS += -D__MACOS__ 
    MSL_STATIC=F
endif

ifeq ($(MSL_STATIC),T)
    FLAGS   += -static
endif

FLAGS   +=  -DUSE_REAL_EQ_DOUBLE 

# Include local Makefile
MYDIR=myProgs
-include myProgs/myProgs.mk

# -include Makefile.local

# Include local Makefile
-include examples/examples.mk
# -include Makefile.local

# Add poroper suffix
OBJECTS       = $(patsubst %,objs/%.o, $(SOURCE)) 
BINARIES      = $(patsubst %,bin/%, $(PROGRAMS)) 
EXAMPLEBINS   = $(patsubst %,bin/%, $(EXAMPLES)) 
SANDBOXBIN    = $(patsubst %,bin/%, $(SANDBOX)) 
GOLDBIN       = $(patsubst %,bin/%, $(GOLD)) 
LEADBIN       = $(patsubst %,bin/%, $(LEAD)) 
PHEADERS      = $(patsubst %,programs/%.h, $(PROGRAMS_HEADERS))

# Include myProg subdirectories
MYBINS        = $(patsubst %,bin/%, $(MYPROGS)) 
MYOBJS        = $(patsubst %,objs/%.o, $(MYSOURCE)) 
MYHEADERFILES = $(patsubst, %,myProgs/%, $(MYHEADERS))


# Compile/Link commands
# all: ${BINARIES} ${MYBINS} ${SANDBOXBIN} ${GOLDBIN} ${EXAMPLEBINS} ${LEADBIN}
all: ${BINARIES} ${MYBINS} ${GOLDBIN} ${EXAMPLEBINS}
programs: ${BINARIES}
objects: ${OBJECTS}
sandbox: ${SANDBOXBIN}
gold: ${GOLDBIN}
lead: ${LEADBIN}
examples: ${EXAMPLEBINS}
mybins: ${MYBINS}
test: ${GOLDBIN}
	bin/testEZpotential 
#	bin/testCharmmEnergies 
	bin/testCharmmBuild

${OBJECTS}: objs/%.o : src/%.cpp src/%.h 
	${CC} ${FLAGS} -I${INCLUDE} ${SYMBOLS} -c $< -o $@ 

${SANDBOXBIN}: bin/% : tests/sandbox/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS}
	${CC} ${FLAGS} ${LINKFLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread

${GOLDBIN}: bin/% : tests/gold/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS}
	${CC} ${FLAGS} ${LINKFLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread

${LEADBIN}: bin/% : tests/lead/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS}
	${CC} ${FLAGS} ${LINKFLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread

${BINARIES}: bin/% : programs/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS} ${PHEADERS}
	${CC} ${FLAGS} ${LINKFLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread

${EXAMPLEBINS}: bin/% : examples/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS} ${PHEADERS}
	${CC} ${FLAGS} ${LINKFLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread

${MYOBJS}: objs/%.o : ${MYDIR}/%.cpp ${MYDIR}/%.h 
	${CC} ${FLAGS} -I${INCLUDE} ${SYMBOLS} -c $< -o $@

${MYBINS}: bin/% : ${MYDIR}/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS} ${MYHEADERFILES}
	${CC} ${FLAGS} ${LINKFLAGS} -I${INCLUDE} -o $@  ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS}  -lpthread

.PHONY : clean
clean :
	-rm -f ${OBJECTS} ${BINARIES} ${EXAMPLEBINS} ${SANDBOXBIN} ${GOLDBIN} ${LEADBIN} ${MYOBJS} ${MYBINS}


pythonLin:
	gcc -fpic -c src/PythonMSL.cpp -o objs/PythonMSL.o -Wall -I${INCLUDE} -I/usr/include/python2.6 -I/usr/include  -I/Library/Frameworks/Python.framework/Versions/2.6/Headers/
	g++ -lm -shared objs/PythonMSL.o ${OBJECTS} ${STATIC_LIBS} -o PythonMSL.so 
	cp PythonMSL.so /usr/share/python-support/pymol/pymol/


pythonMac:
	gcc -fPIC -c src/PythonMSL.cpp -o objs/PythonMSL.o -Wall -I${INCLUDE} -I/usr/include/python2.6 -I/usr/include  -I/System/Library/Frameworks/Python.framework/Versions/2.5/Headers/
	g++ -bundle -undefined dynamic_lookup objs/PythonMSL.o ${OBJECTS} ${STATIC_LIBS} -o PythonMSL.so 
