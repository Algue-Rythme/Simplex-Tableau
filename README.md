# Simplex-Tableau

Implementation of Simplex Tableau algorithm, for the Approximation and Optimisation course at ENS Lyon 2017-2018.

## Usage

All the options can be displayed using "./bethune-exec --help".
The usage of .bethune-exec is the following :

```
./bethune-exec --verbose 3 --pivot steepest --no-cycles bethune-test_devotail.txt
```

will use steepest edge rule (modified by cycle elimination) and print everything (verbose 3).
Defaults values are also provided for options, such as maxCoeff or verbose=0.

The executable was compiled for Ubuntu 16 on x86-64 bits system, with g++ 5.4.0

## Language and librairies

The program was written in C++14 (fully standard compliant) in the file simplex.cpp.
generate_tests.cpp is just an utility to creates big files of tests.

I used Boost for command line parsing. It requires libboost-all-dev being installed.
I also used Boost to get the mpq_rational from libgmp, the GNU Multiple PrecisionArithmetic Library.
Boost provides overloaded operators for this purpose : libGMP need to be present on the computer.
Finally, I used Eigen, a c++ header-only library of linear algebra, to work with matrix or arrays like in Numpy.

I provided a Makefile to compile the code, assuming dependencies are on the computer.

## Pivots and cycling

I provide 3 differents pivot rules.

1. Bland's rule : the simplest pivot rule. We know it never cycles.
2. Maxcoefff : a simple heuristic, may cycles
3. Steepest : one of the best pivot rule according to studies (in average leads to the fewer number of pivots on real problems), may cycles

To avoid cyles in "maxcoeff" and "steepest edge" I provide an option named `--no-cycle`
This option tries to detect cycles, and when it find one, transform the current pivot rule into a Bland's rule.
When we are sure that we are not in cycle anymore (i.e when the optimum improvement is positive) we come back to the previous rule.

How detect a cycle ? There is many ways !

1. Simplest way (implemented) : there is cycling only if (necessary condition) the improvement (delta) is null.
In this case, we store the list of variables that enters or left the basis, in the order.
This is a word on a alphabet where the symbols are the (x_i, x_j).
If the word is w=uvv we just detect a cycle of size |v| where v is the cycle !
u is just the prefix, because we can make several steps with delta=0 without being in a cycle.
Check this is in O(|w|^2) in worst case (we make an iteration on the all possible |v| such that |v| <= |w|/2).
But, in average, it's O(|w|) because it's really unlikely that each test of a given |v| takes O(|v|).
Avantages and inconvenients :
    + easy to code
    + very fast for small cycles
    + may be recquire a lot of memory and computation time when cycle becomes too long.

In worst case the cycle may be exponential in n, because the number of different basis may exponential in n.

2. Better way (no implemented). I suggest to use the Floyd's Tortoise and Hare.
Hare makes two iteration of simplex at each step. Tortoise only one.
I can update both basis in O(1) (easy). With this, I update a error counter in O(1).
When error counter reaches 0 (in O(1)), the basis are identical. And we successfully detect a cycle !
Avantages and inconvenients :
    + Worst case : factor two of memory or computation time, (because of the second simplex we solve in parallel). That is asymptotically identical.
    + Hard to code.
    + Too heavy for small cycles.

I didn't implement this one.

## Exemples

Of course I used the exemples of the homeworks, but also my ones.

Several exemples are from a book on linear programming written by Vanderbei.
The first number is the n° of the page, and the second one the n° of the exercise.
Ex: test_vanderbei_24_22 can be found in page 24 as the 2.2 exemple.
Unfortunately I don't provide any solution reference for these inputs. But my program should work on them.

I made many big files to test speed of my program. Speed is measured and printed on standard output.
It's HIGHLY recommended to set verbose=0 or 1 on big instances.

When m*n <= 10*1000 the average computation time is pretty slow (less than one second).
But when mn >= 1000*1000 we can observe many behaviors.

First, depending of if m or n is the biggest we don't have the same results, depending of the pivot rule used.
Big number of variables seems to lead us to the worst case.

bethune-test_gentest_n_m are dense test cases with very few zeros, and positive coefficients.
bethune-test_sparse_n_m is corresponding to set cases with many zeros, and positive coefficients.
bethune-test_negpos_n_m is also corresponding to set cases with many zeros, but I authorize negative numbers to test Phase 1/Phase 2.

bethune-test_gentest_1000_1000 =>
- 1000 variables 1000 constraints but null optimal solution,
- only 4 pivots with maxcoeff, ~8.2s on my computer
- 7 pivots with steepest, ~13.3s

bethune-test_sparse_100_10000 =>
- 100 variables 10 000 constraints,
- 4 pivots with maxcoeff, ~3.57s
- 4 pivots with steepest edge, ~3.54s

On big size instances, we can conclude steepest edge is not slower (in average) per iteration.

Big size instances are hard to characterize :

bethune-test_negpos_100_1000 => unfeasible, 1 pivot, less than ~1s
bethune-test_negpos_1000_100 => unbounded, 2 pivots (bland) 3 pivots (with maxcoeff or steepest) and ~0.33s
bethune-test_negpos_1000_1000 => too many pivots steps no results found in reasonable time

On very big sparse positive instances :

bethune-test_sparse_100_10000 =>
- 100 variables 10 000 constraints, null optimal solution
- 4 pivot steps regardless pivot rule used, 3.5s

bethune-test_sparse_10000_100 =>
- 10 000 variables 100 constraints, unbounded
- 202 pivot steps with maxcoeff, 174.2s
- 20 pivot steps with bland, 17.3s
- cycle with steepest, to the end of the time
- 115 pivot steps with steepest and no-cycles option, 99.2s

## Bugs

Actually there is one known (theorically) bug.  
When we aplly Phase 1 / Phase 2 algorithm, at the end of phase 1, I immediatly drop the artificial variables.  
But artificial variables may remains in basis (despite a null-optimum) because the optimal solution may be degenerate.  
In this case the algorithm may fails because the basis looses a variable. I don't know what happen in this case, but I expect at least
inconsistent results, and maybe a segfault.  

I could solve it by one of the following to solutions :

+ Keep artificials variables during phase 2. It's ugly and lead to bad performances, but easy to code.  
+ Replace artificials variable in basis by other variables (we know that a such variable exists) because we're on a degenerate solution. It's the best solution but I didn't implemented it.
