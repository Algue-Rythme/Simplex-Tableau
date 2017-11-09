all:
	g++ -std=c++14 -O2 -s -pedantic-errors -Wall -fopenmp simplex.cpp -o lpsolve -lgmp

tests:
	g++ -std=c++14 -O2 -s -pedantic-errors -Wall generate_tests.cpp -o gentest -lgmp

assembly:
	g++ -std=c++14 -O2 -s -pedantic-errors -Wall -S simplex.cpp -o asm_simplex
