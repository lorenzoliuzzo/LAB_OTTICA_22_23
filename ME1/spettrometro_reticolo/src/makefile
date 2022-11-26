PHYSIM_PATH = -I ../../../tools
CXXFLAGS = -g -Wall --pedantic -fno-elide-constructors -std=c++20 -O3 -pg

run: main.cpp
	g++ $(CXXFLAGS) $(PHYSIM_PATH) $< -o $@
	./run

clean:
	rm -f run
