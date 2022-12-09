CXX = g++
CXXFLAGS = -g -Wall --pedantic -fno-elide-constructors -std=c++20 -O0 -pg
PHYSIM_PATH = -I ../../../tools

formulas.o: formulas.hpp
	$(CXX) $(CXXFLAGS) $(PHYSIM_PATH) $< -c -o $@

main: main.cpp formulas.o
	$(CXX) $(CXXFLAGS) $(PHYSIM_PATH) $< -o $@

execute: main
	./main
	
clean:
	rm -f formulas main compile 
