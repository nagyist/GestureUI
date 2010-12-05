all: libseqdb.so

CC = gcc
CFLAGS = -O2

SRC = $(shell ls *.cpp)
OBJ = $(SRC:.cpp=.o)

%.o: %.cpp
	$(CC) -fPIC -c $^

libseqdb.so: $(OBJ)
	$(CC) -O2 -shared -fPIC -Wl,-soname,$@.1 -o $@.1 $^
	ln -s $@.1 $@ 
	ln -s $@ ../bin/$@

clean:
	rm *.o libseqdb.so