all:
	g++ -std=c++14 -O2 -s -pedantic-errors -Wall -fopenmp bethune-simplex.cpp -o bethune-exec -lgmp -lboost_program_options

tests:
	g++ -std=c++14 -O2 -s -pedantic-errors -Wall generate_tests.cpp -o gentest -lgmp
