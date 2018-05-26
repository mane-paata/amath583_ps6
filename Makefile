#
# This file is part of the course materials for AMATH483/583 at the University of Washington,
# Spring 2018
#
# Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
# https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# Author: Andrew Lumsdaine, Tommaso Buvoli
#

CXX             = c++

OPTS            = -march=native -DNDEBUG # Also try -Ofast -march=native
LANG            = -std=c++11 -fopenmp
PICKY           = -Wall 

CXXFLAGS        = $(LANG) $(OPTS) $(PICKY) 

.PHONY: clean all

all: pt2n_driver.exe coo_driver.exe csr_driver.exe

# Generic rules
%.exe	     : %.o
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o          : %.cpp
	$(CXX) $(CXXFLAGS) -c $<


# exe dependencies -- consequent is handled by generic rule
pt2n_driver.exe: amath583.o 
coo_driver.exe: amath583.o 
csr_driver.exe: amath583.o 

# Object file (.o) dependencies -- consequent is handled by generic rule
pt2n_driver.o : amath583.hpp Vector.hpp 
amath583.o   : amath583.hpp Vector.hpp 
coo_driver.o : amath583.hpp COOMatrix.hpp
csr_driver.o : amath583.hpp CSRMatrix.hpp

clean:
	/bin/rm -f *.o *.exe *.gch *~ 
