# Simplex-Tableau

## Usage

All the options can be displayed using "./lpsolve --help".  
The usage of .lpsolve is the following :  

./lpsolve --verbose 2 --pivot steepest test_devotail.txt

will use steepest edge rule and print evrything (verbose 2).  
Defaults values are also provided for options, such as maxCoeff or verbose=0.  

## Language and librairies  

The program was written in C++14 (fully standard compliant) in the file simplex.cpp.
generate_tests.cpp is just an utility to creates big files of tests.  

I used Boost for command line parsing.
I also used Boost to get the mpq_rational from libgmp, the GNU Multiple PrecisionArithmetic Library.
Boost provides overloaded operators for this purpose.  
Most of sublibs of Boost are heade-only but these two one needs additional dependencies.  

Finally, I used Eigen, a c++ header-only library of linear algebra, to work with matrix or arrays like in Numpy.  

I provided a Makefile to compile the code, assuming dependencies are on the computer.

## Exemples  

Of course I used the exemples of the homeworks, but also my ones.  

Several exemples are from a book on linear programming written by Vanderbei.  
The first number is the n° of the page, and the second one the n° of the exercise.  
Ex: test_vanderbei_24_22 can be found in page 24 as the 2.2 exemple.
Unfortunaltely I don't provide any solution reference for these inputs. But my program should work on them.  

I made many big files to test speed of my program. Speed is measured and printed on standard output.  

When m*n <= 10*1000 the average computation time is pretty slow (less than few seconds).  
But when mn >= 1000*1000 we can observe many behaviors.  

First, depending of if m or n is the biggest we don't have the same results, depending of the pivot rule used.  

Big number of variables seems to lead us to the worst case.
