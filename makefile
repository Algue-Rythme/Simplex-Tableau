all:
	g++ -std=c++14 -O2 -s -pedantic-errors -Wall -fopenmp simplex.cpp -o lpsolve -lgmp -lboost_program_options

tests:
	g++ -std=c++14 -O2 -s -pedantic-errors -Wall -fopenmp generate_tests.cpp -o gentest -lgmp -lboost_program_options
