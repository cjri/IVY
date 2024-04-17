CC	    = g++
CC_FLAGS	= -g3 -O3 -Wall -pthread  -I /opt/homebrew/Cellar/gsl/2.7.1/include/
LD_FLAGS	= -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas -lm -lstdc++  -g 
LD_FLAGS_TEST   = -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgtest -lgtest_main -pthread -lgsl 
BAS	     = error.o io.o modelI.o model.o model_tools.o clustering.o likelihood.o process_sequences.o statespace.o ivy.o

basic: $(BAS)
	$(CC) $(CC_FLAGS) $(BAS) -o run_ivy  $(LD_FLAGS)
ivy.o: ivy.cpp
	$(CC) $(CC_FLAGS) -c ivy.cpp
model.o: model.cpp
	$(CC) $(CC_FLAGS) -c model.cpp
modelI.o: modelI.cpp
	$(CC) $(CC_FLAGS) -c modelI.cpp
model_tools.o: model_tools.cpp
	$(CC) $(CC_FLAGS) -c model_tools.cpp
clustering.o: clustering.cpp
	$(CC) $(CC_FLAGS) -c clustering.cpp
likelihood.o: likelihood.cpp
	$(CC) $(CC_FLAGS) -c likelihood.cpp
process_sequences.o: process_sequences.cpp
	$(CC) $(CC_FLAGS) -c process_sequences.cpp
statespace.o: statespace.cpp
	$(CC) $(CC_FLAGS) -c statespace.cpp
error.o: error.cpp
	$(CC) $(CC_FLAGS) -c error.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp

