nveSim: NVE.f90
	gfortran -o nveSim NVE.f90

testPhaseReader: testPhaseReader.cpp phaseReader.cpp
	clang++ -Wall -Wextra -std=c++11 -O2 testPhaseReader.cpp phaseReader.cpp -o testPhaseReader 

nvtSim: NVT.f90
	gfortran -o nvtSim NVT.f90

argonBox: argonBox.f90
	gfortran -o argonBox argonBox.f90

force1: force1.f90
	gfortran -o force force1.f90

force: force.cpp
	clang++ -Wall -Wextra -std=c++11 -O2 force.cpp
