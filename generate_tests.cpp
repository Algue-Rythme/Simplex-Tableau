#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>

#include <boost/multiprecision/gmp.hpp>

using namespace std;
using namespace boost::multiprecision;

typedef mpq_rational Rational;

const int rate = 100;

void randomLine(
    ofstream& out,
    default_random_engine& gen,
    uniform_int_distribution<int>& dis,
    int n) {
    for (int c = 0; c < n; ++c) {
        if (gen()%rate == 0)
            out << dis(gen) << " ";
        else
            out << "0 ";
    }
    out << "\n";
}

int main(int argc, char * argv[])  {
    default_random_engine gen(42);
    uniform_int_distribution<int> pos(0, 1000);
    uniform_int_distribution<int> negpos(-1000, 1000);

    ios::sync_with_stdio(false);
    if (argc < 4) {
        cout << "Please provide n and m\n";
        return EXIT_FAILURE;
    }
    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    string name(argv[3]);
    ofstream out(name);
    out << n << "\n" << m << "\n";
    randomLine(out, gen, negpos, n);
    randomLine(out, gen, negpos, m);
    for (int i = 0; i < m; ++i) {
        randomLine(out, gen, negpos, n);
    }
}
