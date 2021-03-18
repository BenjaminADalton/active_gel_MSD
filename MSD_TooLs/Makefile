EX_NAME = execute_MSD_TooLs
CPP_FLAGS = -c -std=c++17 -O3

OBJS = main.o filament.o utilities.o

$(EX_NAME) : $(OBJS)
	g++ -o $(EX_NAME) $(OBJS)

main.o : main.cpp
	g++ $(CPP_FLAGS) main.cpp
filament.o : filament.cpp filament.h
	g++ $(CPP_FLAGS) filament.cpp
utilities.o : utilities.cpp utilities.h
	g++ $(CPP_FLAGS) utilities.cpp

clean:
	rm -f core $(EX_NAME) $(OBJS)
