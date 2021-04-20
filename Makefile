# Local Path
SRCPATH	  = ./src/
BINPATH	  = ./bin/
INCPATH	  = ./include/

# Mac Path

MACPATH   = /Library/gurobi901/mac64
INCMAC    = $(MACPATH)/include/
INCSPDLOG = /usr/local/Cellar/spdlog/include/
CPPLIBMAC = -L$(MACPATH)/lib/ -lgurobi_c++ -lgurobi90 $(CPPSTDLIB) -lpthread -lm

# Linux Path
GRBPATH   = /opt/gurobi911/linux64
INCGRB    = $(GRBPATH)/include/
CPPLIBGRB = -L$(GRBPATH)/lib/ -lgurobi_g++5.2 -lgurobi91 $(CPPSTDLIB) -lpthread -lm

# C++ Compiler and Flags
CXX=g++

# MacOS Flags
# CXXFLAGS= -m64 -std=c++17 -O3 -Iinclude/ -I$(INCMAC) 

# Linux Flags
CXXFLAGS=-std=c++17 -O3 -Iinclude/ -I$(INCSPDLOG) -I$(INCGRB) $(CPPLIBGRB)

main: main.o Graph.o Formulation.o MinHeap.o PathSeparation.o TreeSeparation.o
	$(CXX) $(CXXFLAGS) -o $(BINPATH)main main.o Graph.o Formulation.o MinHeap.o PathSeparation.o TreeSeparation.o $(CPPLIBGRB)
main.o: 
	$(CXX) $(CXXFLAGS) -c $(SRCPATH)main.cpp 
Graph.o:
	$(CXX) $(CXXFLAGS) -c $(SRCPATH)Graph.cpp
Formulation.o:
	$(CXX) $(CXXFLAGS) -c $(SRCPATH)Formulation.cpp
MinHeap.o:
	$(CXX) $(CXXFLAGS) -c $(SRCPATH)MinHeap.cpp
PathSeparation.o:
	$(CXX) $(CXXFLAGS) -c $(SRCPATH)PathSeparation.cpp
TreeSeparation.o:
	$(CXX) $(CXXFLAGS) -c $(SRCPATH)TreeSeparation.cpp
clean:
	rm -rf *.o ./bin/main
