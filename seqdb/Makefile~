all: seqdb_test

CC = g++
CFLAGS = -O2

SRC = $(shell ls *.cpp)
OBJ = $(SRC:.cpp=.o)

%.o: %.cpp
	$(CC) -fPIC -c $^

seqdb_test: $(SRC)
	$(CC) $^ -o $@

clean:
	rm *.o seqdb_test