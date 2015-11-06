CC = g++

FLAGS = -std=c++11 -g -Wall

MAIN = main

main.out: Vector.h Matrix.h main.cpp r8mat_expm1.o r8lib.o
	$(CC) $(FLAGS) $(MAIN).cpp r8mat_expm1.o r8lib.o -o $(MAIN).out

r8mat_expm1.o: r8mat_expm1.cpp
	$(CC) $(FLAGS) -c r8mat_expm1.cpp -o r8mat_expm1.o

r8lib.o: r8lib.cpp
	$(CC) $(FLAGS) -c r8lib.cpp -o r8lib.o

clean:
	rm -f *.o *.out && clear 

