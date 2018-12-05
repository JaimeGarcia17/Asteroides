CC=g++
CFLAGS= -std=c++14 -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -fopenmp -O3 -DNDEBUG -Wno-uninitialized



nasteroids-par: nasteroids-par.o
	$(CC) nasteroids-par.o -o nasteroids-par -fopenmp
	@echo "***************************"
	@echo "Compilation successfully!"
	@echo "***************************"
nasteroids-par.o: nasteroids-par.cpp
	$(CC) $(CFLAGS) -c nasteroids-par.cpp
clean:
	rm -f *.o *.txt nasteroids-par
	@echo "***************************"
	@echo "Deleted files!"
	@echo "***************************"
