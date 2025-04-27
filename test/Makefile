CC=g++
CFLAGS=-g -Wall -Wextra
INCLUDES=-I..
GTEST_DIR?=/usr/src/gtest

all: out

# Build gtest libraries
gtest-all.o: $(GTEST_DIR)/src/gtest-all.cc
	$(CC) -isystem $(GTEST_DIR)/include -I$(GTEST_DIR) -pthread -c $< -o $@

libgtest.a: gtest-all.o
	ar -rv $@ $^

# Build CellCDecon from source
CellCDecon.o: ../CellCDecon.cpp ../CellCDecon.h
	$(CC) $(CFLAGS) -c ../CellCDecon.cpp -o $@

# Build CellCDeconIO from source
CellCDeconIO.o: ../CellCDeconIO.cpp ../CellCDeconIO.h
	$(CC) $(CFLAGS) -c ../CellCDeconIO.cpp -o $@

# Build test file
test_cellc.o: test_cellc.cc test_cellc.h sample15.h
	$(CC) -c $(CFLAGS) $(INCLUDES) -isystem $(GTEST_DIR)/include -pthread $< -o $@

# Build the test executable
out: test_cellc.o libgtest.a CellCDecon.o CellCDeconIO.o
	$(CC) $(CFLAGS) -pthread $^ -o $@

# Add dedicated target for test_cellc
test_cellc: test_cellc.o libgtest.a CellCDecon.o CellCDeconIO.o
	$(CC) $(CFLAGS) -pthread $^ -o $@

clean:
	rm -f *.o *.a out test_cellc

.PHONY: all clean 