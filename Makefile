CC=g++
CFLAGS=-O3
DEPS = CellCDecon.h
OBJ = CellCDecon.o CellCMain.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

CellCDecon: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)