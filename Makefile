CC := g++
CFLAGS := -g -std=c++11
LDFLAGS := -lm #link editing: math lib
exe = bdd
exe2 = bdd_simple
path = src

all:$(path)/main.cpp $(path)/BDD.hpp $(path)/truth_table.hpp
	@$(CC) $(path)/main.cpp -o $(exe) $(CFLAGS) $(LDFLAGS)
	make clean

simple:$(path)/simple.cpp $(path)/BDD.hpp $(path)/truth_table.hpp
	@$(CC) $(path)/simple.cpp -o $(exe2) $(CFLAGS) $(LDFLAGS)
	make clean

clean:
	@rm -rf *.o *.dSYM
	#@rm -rf *.o *.dSYM $(exe) $(exe2)
