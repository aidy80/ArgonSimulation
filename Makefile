argonSim: Ft_argonSim.f90
	gfortran -o argonSim Ft_argonSim.f90

argonBox: Ft_argonBox.f90
	gfortran -o argonBox Ft_argonBox.f90

force1: force1.f90
	gfortran -o force force1.f90

force: force.cpp
	clang++ -Wall -Wextra -std=c++11 -O2 force.cpp
