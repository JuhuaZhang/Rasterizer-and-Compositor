.PHONEY: build, run

build: program

run: program
	./program $(file)

program: main.cpp
	g++ -std=c++11 -O3 -I. -g main.cpp lodepng.cpp -Wall -Wextra -pedantic -o program
