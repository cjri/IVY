CC	    = g++
CC_FLAGS	= -g3 -O3 -Wall -pthread  -I /opt/homebrew/Cellar/gsl/2.7.1/include/
LD_FLAGS	= -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas -lm -lstdc++  -g 
LD_FLAGS_TEST   = -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgtest -lgtest_main -pthread -lgsl -lgslcblas -g
BAS	     = io.o utilities.o modelI.o modelV.o modelX.o modelY.o seqmodel.o statespace.o

basic: $(BAS)
	$(CC) $(CC_FLAGS) $(BAS) -o run_ivy  $(LD_FLAGS)
seqmodel.o: seqmodel.cpp
	$(CC) $(CC_FLAGS) -c seqmodel.cpp
modelI.o: modelI.cpp
	$(CC) $(CC_FLAGS) -c modelI.cpp
modelV.o: modelV.cpp
	$(CC) $(CC_FLAGS) -c modelV.cpp
modelY.o: modelY.cpp
	$(CC) $(CC_FLAGS) -c modelY.cpp
modelX.o: modelX.cpp
	$(CC) $(CC_FLAGS) -c modelX.cpp
statespace.o: statespace.cpp
	$(CC) $(CC_FLAGS) -c statespace.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp

