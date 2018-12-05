CC=g++
CFLAGS= -std=c++14 -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -O3 -DNDEBUG



nasteroids-seq: nasteroids-seq.o
	$(CC) nasteroids-seq.o -o nasteroids-seq
	@echo "***************************"
	@echo "Compilation successfully!"
	@echo "***************************"
nasteroids-seq.o: nasteroids-seq.cpp
	$(CC) $(CFLAGS) -c nasteroids-seq.cpp
clean:
	rm -f *.o *.txt nasteroids-seq
	@echo "***************************"
	@echo "Deleted files!"
	@echo "***************************"
