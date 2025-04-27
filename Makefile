CC=g++
CFLAGS=-O3 -std=c++17 -Wall -Wextra
DEPS = CellCDecon.h CellCDeconIO.h
OBJ = CellCDecon.o CellCDeconIO.o CellCMain.o 

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

CellCDecon: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -f *.o CellCDecon