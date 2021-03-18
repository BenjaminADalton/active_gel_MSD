# Makefile for program: Filament Dynamcis
# Benjamin Dalton: 18-09-2016

EX_NAME = execute_dynamics
CPP_FLAGS = -c -std=c++17 -O3

OBJS = main.o filament.o crosslinker.o pair_struct.o pbc.o binding.o utilities.o wall.o wall_calc.o 
$(EX_NAME) : $(OBJS)
	g++ -o $(EX_NAME) $(OBJS)

main.o : main.cpp
	g++ $(CPP_FLAGS) main.cpp
filament.o : filament.cpp filament.h
	g++ $(CPP_FLAGS) filament.cpp
crosslinker.o : crosslinker.cpp crosslinker.h
	g++ $(CPP_FLAGS) crosslinker.cpp
pair_struct.o : pair_struct.cpp pair_struct.h
	g++ $(CPP_FLAGS) pair_struct.cpp
pbc.o : pbc.cpp pbc.h
	g++ $(CPP_FLAGS) pbc.cpp
utilities.o : utilities.cpp utilities.h
	g++ $(CPP_FLAGS) utilities.cpp
binding.o : binding.cpp binding.h
	g++ $(CPP_FLAGS) binding.cpp
wall.o : wall.cpp wall.h
	g++ $(CPP_FLAGS) wall.cpp
wall_calc.o : wall_calc.cpp wall_calc.h
	g++ $(CPP_FLAGS) wall_calc.cpp

clean:
	rm -f core $(EX_NAME) $(OBJS)
