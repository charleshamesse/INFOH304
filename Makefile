all: LinkedList.o main.cpp
	g++ LinkedList.o main.cpp -o main

LinkedList.o: LinkedList.cpp
	g++ -c LinkedList.cpp

clean:
	rm *.o
