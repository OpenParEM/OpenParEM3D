DEBUG ?= 1

ifeq ($(DEBUG), 1)
   CFLAGS=-Wall -g
   CxxFLAGS=-Wall -std=c++17 -g -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(SLEPC_DIR)/$(PETSC_ARCH)/lib
else
   CFLAGS=-Wall -O3 
   CxxFLAGS=-Wall -std=c++17 -O3 
endif


CCxx=mpicxx
CxxINCS=-I../../OpenParEMCommon/src -I$(MFEM_DIR) -I$(MFEM_DIR)/linalg -I$(HYPRE_DIR)/src/hypre/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(EIGEN_DIR)
CxxLDIR=-L../../OpenParEMCommon/lib -L$(MFEM_DIR) -L$(HYPRE_DIR)/src/hypre/lib -L$(METIS_DIR) -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/9
CxxLIBS=-lOpenParEMCommon -lpetscmat -lpetscsnes -lpetscdm -lpetscvec -lpetscts -lpetsctao -lpetscsys -lpetscksp -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lflapack -lfblas -lptesmumps -lptscotchparmetisv3 -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -lm -lX11 -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lmfem -lHYPRE -lmetis -lrt -lslepcpep -lslepcsys -lslepceps -lslepclme -lslepcnep -lslepcmfn -lslepcsvd /usr/lib/x86_64-linux-gnu/liblapacke64.a /usr/lib/x86_64-linux-gnu/liblapack64.a -lpetscmat -lpetscsnes -lpetscmat -lpetscsnes -lpetscdm -lpetscvec -lpetscts -lpetsctao -lpetscsys -lpetscksp -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lflapack -lfblas -lslepcpep -lslepcsys -lslepceps -lslepclme -lslepcnep -lslepcmfn -lslepcsvd -lgfortran -lpetscmat -lpetscsnes -lpetscdm -lpetscvec -lpetscts -lpetsctao -lpetscsys -lpetscksp 

CC=mpicc
CINCS=-I../../OpenParEMCommon/src -I$(PETSC_DIR)/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(EIGEN_DIR) -I$(HYPRE_DIR)/src/hypre/include
CLDIR=-L../../OpenParEMCommon/lib -L$(HYPRE_DIR)/src/hypre/lib
CLIBS=-lHYPRE -lm

csr.o: csr.c csr.h
	$(CC) $(CFLAGS) -c csr.c $(CINCS)

fem3D.o: fem3D.cpp fem3D.hpp project.h port.hpp OPEM_L2ZZErrorEstimator.hpp
	$(CCxx) $(CxxFLAGS) -c fem3D.cpp $(CxxINCS)

OpenParEM3D.o: OpenParEM3D.cpp fem3D.hpp port.hpp project.h results.hpp
	$(CCxx) $(CxxFLAGS) -c OpenParEM3D.cpp $(CxxINCS)

port.o: port.cpp port.hpp project.h
	$(CCxx) $(CxxFLAGS) -c port.cpp $(CxxINCS)

process3D.o: process3D.cpp process3D.hpp project.h results.hpp
	$(CCxx) $(CxxFLAGS) -c process3D.cpp $(CxxINCS)

OPEM_L2ZZErrorEstimator.o: OPEM_L2ZZErrorEstimator.cpp OPEM_L2ZZErrorEstimator.hpp
	$(CCxx) $(CxxFLAGS) -c OPEM_L2ZZErrorEstimator.cpp $(CxxINCS)

project.o: project.c project.h 
	$(CC) $(CFLAGS) -c project.c $(CINCS)

results.o: results.cpp results.hpp fem3D.hpp
	$(CCxx) $(CxxFLAGS) -c results.cpp $(CxxINCS)

tline3.o: tline3.cpp
	$(CCxx) $(CxxFLAGS) -c tline3.cpp $(CxxINCS)

solveComplexLinearSystem.o: solveComplexLinearSystem.c solveComplexLinearSystem.h project.h
	$(CC) $(CFLAGS) -c solveComplexLinearSystem.c $(CINCS)

process3D: process3D.o csr.o solveComplexLinearSystem.o port.o fem3D.o project.o results.o OPEM_L2ZZErrorEstimator.o
	$(CCxx) $(CxxFLAGS) csr.o solveComplexLinearSystem.o port.o OPEM_L2ZZErrorEstimator.o fem3D.o process3D.o project.o results.o -o process3D $(CxxLDIR) $(CxxLIBS)

simplify3D.o: simplify3D.cpp simplify3D.hpp project.h
	$(CCxx) $(CxxFLAGS) -c simplify3D.cpp $(CxxINCS)

OpenParEM3D: OpenParEM3D.o csr.o solveComplexLinearSystem.o port.o fem3D.o project.o results.o OPEM_L2ZZErrorEstimator.o
	$(CCxx) $(CxxFLAGS) csr.o solveComplexLinearSystem.o OpenParEM3D.o port.o OPEM_L2ZZErrorEstimator.o fem3D.o project.o results.o -o OpenParEM3D $(CxxLDIR) $(CxxLIBS)

tline3: tline3.o
	$(CCxx) $(CxxFLAGS) tline3.o -o tline3 $(CxxLDIR) $(CxxLIBS)

simplify3D: simplify3D.o project.o
	$(CCxx) $(CxxFLAGS) simplify3D.o project.o -o simplify3D $(CxxLDIR) $(CxxLIBS)


.PHONY: all clean

all: OpenParEM3D process3D tline3 simplify3D

install:
	rm -f ../bin/OpenParEM3D
	cp OpenParEM3D ../bin/OpenParEM3D
	rm -f $(HOME)/bin/OpenParEM3D
	cp OpenParEM3D $(HOME)/bin/OpenParEM3D

	rm -f ../bin/process3D
	cp process3D ../bin/process3D
	rm -f $(HOME)/bin/process3D
	cp process3D $(HOME)/bin/process3D

	rm -f ../bin/tline3
	cp tline3 ../bin/tline3
	rm -f $(HOME)/bin/tline3
	cp tline3 $(HOME)/bin/tline3

	rm -f ../bin/simplify3D
	cp simplify3D ../bin/simplify3D
	rm -f $(HOME)/bin/simplify3D
	cp simplify3D $(HOME)/bin/simplify3D

	rm -f $(HOME)/bin/process3D.sh
	cp ../scripts/process3D.sh $(HOME)/bin/process3D.sh

	rm -f $(HOME)/bin/regression3D.sh
	cp ../scripts/regression3D.sh $(HOME)/bin/regression3D.sh

	rm -f $(HOME)/.FreeCAD/Macro/OpenParEM3D_save.py
	cp ../scripts/OpenParEM3D_save.py $(HOME)/.FreeCAD/Macro/OpenParEM3D_save.py

install_pre-compiled:
	rm -f $(HOME)/bin/OpenParEM3D
	cp ../bin/OpenParEM3D $(HOME)/bin/OpenParEM3D

	rm -f $(HOME)/bin/process3D
	cp ../bin/process3D $(HOME)/bin/process3D

	rm -f $(HOME)/bin/tline3
	cp ../bin/tline3 $(HOME)/bin/tline3

	rm -f $(HOME)/bin/simplify3D
	cp ../bin/simplify3D $(HOME)/bin/simplify3D

	rm -f $(HOME)/bin/process3D.sh
	cp ../scripts/process3D.sh $(HOME)/bin/process3D.sh

	rm -f $(HOME)/bin/regression3D.sh
	cp ../scripts/regression3D.sh $(HOME)/bin/regression3D.sh

	rm -f $(HOME)/.FreeCAD/Macro/OpenParEM3D_save.py
	cp ../scripts/OpenParEM3D_save.py $(HOME)/.FreeCAD/Macro/OpenParEM3D_save.py

clean:
	rm -f OpenParEM3D
	rm -f process3D
	rm -f tline3
	rm -f simplify3D
	rm -f OpenParEM3D.o
	rm -f project.o
	rm -f tline3.o
	rm -f simplify3D.o
	rm -f solveComplexLinearSystem.o
	rm -f csr.o
	rm -f fem3D.o
	rm -f results.o
	rm -f port.o
	rm -f process3D.o
	rm -f OPEM_L2ZZErrorEstimator.o
