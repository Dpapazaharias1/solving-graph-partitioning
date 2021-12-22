#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

SRCPATH	        = ./src/
BINPATH	        = ./bin/
OBJPATH	        = ./obj/
RUNPATH         = ./run/
OUTPATH         = ./out/
GRBPATH         = /opt/gurobi911/linux64
GRBLIB          = $(GRBPATH)/lib/
GRBINC          = $(GRBPATH)/include/
MAINOBJFILES	= $(addprefix $(OBJPATH),$(OBJ))

OBJ             =  main.o \
                   Graph.o \
	           Formulation.o \
		   MinHeap.o \
		   GRBSeparation.o

#-----------------------------------------------------------------------------
# Flagss
#-----------------------------------------------------------------------------
CXX            = g++
CXXFLAGS       = -std=c++17 -O3 -I$(GRBINC)
GRBFLAGS       = -L$(GRBLIB) -lgurobi_g++5.2 -lgurobi91 $(CPPSTDLIB) -lpthread -lm 

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

main: $(BINPATH) $(OBJPATH) $(MAINOBJFILES)
	@echo "-> linking $@"
	$(CXX) $(MAINOBJFILES) -o $(BINPATH)$@ $(GRBFLAGS)


$(OBJPATH)%.o:	$(SRCPATH)%.cpp 
	@echo "-> compiling $@"
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

$(OBJPATH):
	@-mkdir -p $(OBJPATH)

$(BINPATH):
	@-mkdir -p $(BINPATH)

run_all: table2 table3 table4 table5-6 table7-8 table9 table10

table2:
	@echo "-> running instances for table 2"
	@bash $(RUNPATH)table2.sh

table3:
	@echo "-> running instances for table 3"
	@bash $(RUNPATH)table3.sh


table4:
	@echo "-> running instances for table 4"
	@bash $(RUNPATH)table4.sh

table5-6:
	@echo "-> running instances for tables 5 and 6"
	@bash $(RUNPATH)table5-6.sh

table7-8:
	@echo "-> running instances for tables 7 and 8"
	@bash $(RUNPATH)table7-8.sh

table9:
	@echo "-> running instances for table 9"
	@bash $(RUNPATH)table9.sh

table10:
	@echo "-> running instances for table 10"
	@bash $(RUNPATH)table10.sh

clean_out:
	rm $(OUTPATH)*.txt

clean:
	rm -rf  $(BINPATH)main $(OBJPATH)*.o
