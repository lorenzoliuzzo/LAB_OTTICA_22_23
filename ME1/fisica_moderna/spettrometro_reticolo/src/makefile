PHYSIM_PATH = -I ~/Desktop/codes/physim/include
CXXFLAGS = -g -Wall --pedantic -fno-elide-constructors -std=c++20 -O0 -pg

compile: main.cpp	
	g++ $(CXXFLAGS) $(PHYSIM_PATH) $< -o $@

clean:
	rm -f compile gmon.out
