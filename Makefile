all: LinkedList.o main.cpp
	g++ -std=c++11 LinkedList.o main.cpp -o main

LinkedList.o: LinkedList.cpp
	g++ -c LinkedList.cpp

clean:
	rm *.o
