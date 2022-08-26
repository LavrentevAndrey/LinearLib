CC=g++
CFLAGS=-c -Wall -Werror -Wextra
SOURCE= main.cpp
SOURCES= matrix.hpp vector.hpp matrix_vector_multiplication.hpp
OBJECT=$(SOURCE:.cpp=.o)
OBJECTS=$(SOURCES:.hpp=.o)
EXECUTABLE=lib

all: $(SOURCE) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECT)
	$(CC) $(OBJECT) -o $@

.hpp.o:
	$(CC) $(CFLAGS) $< -o $@


.PHONY: clean rebuild

clean:
	rm -rf $(OBJECTS) $(OBJECT) $(EXECUTABLE) ./*.o ./*.a ./*.so 

rebuild: clean all