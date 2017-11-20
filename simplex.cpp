#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <chrono>
#include <tuple>

#include <boost/multiprecision/gmp.hpp>
#include <boost/program_options.hpp>
#include <Eigen/Dense>

using namespace std;
using namespace boost::multiprecision;
using namespace Eigen;

typedef mpq_rational Rational;
typedef Array<Rational, 1, Dynamic> RowVectorXr;
typedef Array<Rational, Dynamic, 1> VectorXr;
typedef Array<Rational, Dynamic, Dynamic> ArrayXr;
typedef Array<bool, Dynamic, 1> VectorXb;
typedef Array<bool, 1, Dynamic> RowVectorXb;
typedef vector<int> Basis;

class SimplexTableau {
public:
    enum class PivotType {
        Steepest,
        Bland,
        MaxCoeff
    };

    enum class Phase {
        Phase1,
        Phase2
    };

    enum class ProblemState {
        NotOptimum,
        Unbounded,
        Unfeasible,
        OptimumFound
    };

    SimplexTableau(int, PivotType, bool);

    void read_from_file(const string&);
    void print_problem() const;
    void print_tableau() const;
    void print_solution() const;
    void print_params() const;

    void solve();

private:

    RowVectorXr::Index pivot_entering_DantzigRule() const;
    RowVectorXr::Index pivot_entering_Bland() const;
    RowVectorXr::Index pivot_entering_SteepestEdgeRule() const;
    VectorXr::Index pivot_leaving_Bland(const VectorXr& A_col) const;

    void do_pivot(RowVectorXr::Index, VectorXr::Index);
    void updateBasis(RowVectorXr::Index, VectorXr::Index);
    void clear_cycling_detection();
    void from_phase1_to_phase2();
    ProblemState simplex_iteration();
    void simplex_method(Phase phase);

    typedef tuple<RowVectorXr::Index,RowVectorXr::Index> PivotOperation;

    int n;
    int m;
    int artificials;
    RowVectorXr c;
    RowVectorXr real_cost;
    VectorXr b;
    ArrayXr A;
    Rational opt;
    Basis basis;
    vector<int> affectation;
    unsigned int nbPivots;
    int verbose;
    PivotType pivot;
    PivotType oldPivot;
    Rational deltaOpt;
    bool detectCycles;
    vector<PivotOperation> history;
};

SimplexTableau::SimplexTableau(int _verboseLevel, PivotType _pivot, bool _detectCycles):
opt(0), nbPivots(0), verbose(_verboseLevel),
pivot(_pivot), oldPivot(_pivot), detectCycles(_detectCycles) {}

void SimplexTableau::read_from_file(const string& name) {
    stringstream in;
    ifstream buffer(name);
    in << buffer.rdbuf();
    in >> n >> m;
    c.resize(n + m);
    c.block(0, n, 1, m) = RowVectorXr::Zero(m);
    for (int variable = 0; variable < n; ++variable) {
        in >> c(variable);
    }
    b.resize(m);
    for (int constraint = 0; constraint < m; ++constraint) {
        in >> b(constraint);
    }
    real_cost = c;
    artificials = (b < 0).count();
    A.resize(m, n + m + artificials);
    A.block(0, n, m, m) = Matrix<Rational,Dynamic,Dynamic>::Identity(m, m);
    for (int constraint = 0; constraint < m; ++constraint) {
        for (int variable = 0; variable < n; ++variable) {
            in >> A(constraint, variable);
        }
    }
    print_problem();
    basis.resize(m);
    if (artificials > 0) {
        c.resize(n + m + artificials);
        c.leftCols(n+m).setZero();
        c.rightCols(artificials).setConstant(-1);
    }
    for (int row = 0, i_artificial = 0; row < m; ++row) {
        if (b(row) < 0) {
            b(row) *= -1;
            A.row(row) *= -1;
            A(row, i_artificial + n + m) = 1;
            opt += b(row);
            c += A.row(row);
            basis[row] = n + m + i_artificial;
            i_artificial += 1;
        } else {
            basis[row] = n + row;
        }
    }
    affectation.resize(n, -1);
}

void SimplexTableau::print_problem() const {
    if (verbose < 2)
        return ;
    cout << "OUTPUT\n";
    cout << "The input linear program is:\n\n";
    cout << "Maximise  ";
    bool first = true;
    for (int x = 0; x < n; ++x) {
        if (real_cost(x) != 0) {
            if (!first)
                cout << " + ";
            first = false;
            cout << real_cost(x) << "x_" << (x+1);
        }
    }
    cout << "\nSuch that ";
    for (int c = 0; c < m; ++c) {
        first = true;
        for (int x = 0; x < n; ++x) {
            if (A(c, x) != 0) {
                if (!first)
                    cout << " + ";
                else if (c != 0)
                    cout << string(10, ' ');
                first = false;
                cout << A(c, x) << "x_" << (x+1);
            }
        }
        cout << " <= " << b(c) << "\n";
    }
    first = true;
    for (int x = 0; x < n; ++x) {
        if (!first)
            cout << ", ";
        else
            cout << string(10, ' ');
        first = false;
        cout << "x_" << (x+1);
    }
    cout << " are non-negative" << endl;
}

void SimplexTableau::print_tableau() const {
    if (verbose < 3)
        return ;
    ArrayXr t;
    t.resize(A.rows() + 1, A.cols() + 1);
    t.topLeftCorner(1,  A.cols()) = c;
    t.topRightCorner(1, 1) = opt;
    t.bottomRightCorner(A.rows(), 1) = b;
    t.bottomLeftCorner(A.rows(), A.cols()) = A;
    cout << t << endl;
}

void SimplexTableau::print_solution() const {
    if (verbose >= 0) {
        cout << "One optimal solution is: ";
        for (int variable = 0; variable < n; ++variable) {
            if (variable != 0)
                cout << ", ";
            cout << "x_" << (variable+1) << " = ";
            if (affectation[variable] != -1)
                cout << b(affectation[variable]);
            else
                cout << "0";
        }
    }
    cout << "\nThe value of the objective for this solution is: " << -opt << "\n";
}

void SimplexTableau::print_params() const {
    cout << "The number of pivots is: " << nbPivots << "\n";
    cout << "Pivot rule used: ";
    switch (oldPivot) {
        case PivotType::Bland:
        cout << "Bland's rule ";
        break;
        case PivotType::Steepest:
        cout << "Steepest Edge rule ";
        break;
        default:
        cout << "Dantzig's rule ";
        break;
    }
    if (detectCycles) {
        cout << "with no-cycles option\n";
    } else {
        cout << "\n";
    }
}

RowVectorXr::Index SimplexTableau::pivot_entering_DantzigRule() const {
    RowVectorXr::Index index;
    c.maxCoeff(&index);
    return index;
}

RowVectorXr::Index SimplexTableau::pivot_entering_Bland() const {
    for (int coeff = 0; coeff < c.cols(); ++coeff) {
        if (c(coeff) > 0)
            return coeff;
    }
    throw runtime_error("Try to perform variable entering on non-positive vector of reduced costs : should never happen\n");
}

RowVectorXr::Index SimplexTableau::pivot_entering_SteepestEdgeRule() const {
    RowVectorXr norms = (c>0).select(A.matrix().colwise().squaredNorm(), 0);
    RowVectorXr steep = (norms!=0).select(c/norms, 0);
    RowVectorXr::Index index;
    steep.maxCoeff(&index);
    return index;
}

VectorXr::Index SimplexTableau::pivot_leaving_Bland(const VectorXr& A_col) const {
    auto candidates = (A_col > 0).eval();
    if (!candidates.any())
        return -1;
    VectorXr isLimitant = candidates.select(b/A_col, 0);
    Rational maxAxis = isLimitant.maxCoeff() + 1;
    VectorXr q = candidates.select(isLimitant, maxAxis);
    VectorXr::Index i_leaving;
    q.minCoeff(&i_leaving);
    return i_leaving;
}

void SimplexTableau::do_pivot(RowVectorXr::Index i_entering, VectorXr::Index i_leaving) {
    VectorXr goal = VectorXr::Zero(b.rows());
    goal(i_leaving) = 1;
    VectorXr coeffPivot = (goal - A.col(i_entering)) / A(i_leaving, i_entering);
    A += (coeffPivot.matrix() * A.row(i_leaving).matrix()).array().eval();
    b += (coeffPivot * b(i_leaving)).eval();
    Rational coeffCost = -c(i_entering);
    c += coeffCost * A.row(i_leaving);
    Rational prev = opt;
    opt += coeffCost * b(i_leaving);
    deltaOpt = opt - prev;
}

void SimplexTableau::clear_cycling_detection() {
    history.clear();
    pivot = oldPivot;
}

void SimplexTableau::updateBasis(RowVectorXr::Index i_entering, VectorXr::Index i_leaving) {
    if (i_entering < n)
        affectation[i_entering] = i_leaving;
    if (basis[i_leaving] < n)
        affectation[basis[i_leaving]] = -1;
    basis[i_leaving] = i_entering;
    if (!detectCycles)
        return ;
    if (deltaOpt > 0) {
        clear_cycling_detection();
    } else {
        if (pivot == PivotType::Bland)
            return ;
        history.emplace_back(i_entering, i_leaving);
        if (history.size() < 2)
            return ;
        for (int s = 1; s*2 <= (int)history.size(); ++s) {
            bool cycle = true;
            for (int i = 0; i < s; ++i) {
                if (history[history.size()-i-1] != history[history.size()-s-i-1]) {
                    cycle = false;
                    break ;
                }
            }
            if (cycle) {
                cout << "CYCLE DETECTED !\n";
                oldPivot = pivot;
                pivot = PivotType::Bland;
                return ;
            }
        }
    }
}

SimplexTableau::ProblemState SimplexTableau::simplex_iteration() {
    if (!(c > 0).any())
        return ProblemState::OptimumFound;
    RowVectorXr::Index i_entering;
    switch (pivot) {
        case PivotType::Bland:
        i_entering = pivot_entering_Bland();
        break;
        case PivotType::Steepest:
        i_entering = pivot_entering_SteepestEdgeRule();
        break;
        default:
        i_entering = pivot_entering_DantzigRule();
        break;
    }
    if (verbose >= 2)
        cout << "Entering x_" << (i_entering+1) << endl;
    VectorXr::Index i_leaving = pivot_leaving_Bland(A.col(i_entering));
    if (i_leaving == -1)
        return ProblemState::Unbounded;
    if (verbose >= 2)
        cout << "Leaving x_" << (basis[i_leaving]+1) << endl;
    do_pivot(i_entering, i_leaving);
    updateBasis(i_entering, i_leaving);
    print_tableau();
    return ProblemState::NotOptimum;
}

void SimplexTableau::from_phase1_to_phase2() {
    /*//
    / WARNING : HERE
    / Because of redundants constraints
    / A artificial variable may remain in basis
    / So this next line could lead to serious failure
    //*/
    A = A.leftCols(n + m).eval();
    c = real_cost;
    for (int col = 0; col < n; ++col) {
        int rowOfBasis = affectation[col];
        if (rowOfBasis != -1) {
            Rational coeff = -c(col);
            c += coeff*A.row(rowOfBasis);
            opt += coeff*b(rowOfBasis);
        }
    }
    artificials = 0;
    clear_cycling_detection();
}

void SimplexTableau::simplex_method(Phase phase) {
    if (phase == Phase::Phase1)
        cout << "Phase 1\n";
    else
        cout << "Phase 2\n";
    if (verbose >= 3) {
        cout << "The initial tableau is:\n";
        print_tableau();
    }
    while (true) {
        if (verbose >= 2)
            cout << "Current optimum : " << -opt << "\n";
        ProblemState problem = simplex_iteration();
        if (problem == ProblemState::Unbounded) {
            cout << "The problem is UNBOUNDED !\n";
            break ;
        } else if (problem == ProblemState::OptimumFound) {
            if (phase == Phase::Phase1) {
                if (opt != 0) {
                    cout << "The problem is UNFEASIBLE !\n";
                    break ;
                } else {
                    cout << "The problem is FEASIBLE !\n";
                    from_phase1_to_phase2();
                    simplex_method(Phase::Phase2);
                    return ;
                }
            }
            print_solution();
            break ;
        }
        nbPivots += 1;
    }
    print_params();
}

void SimplexTableau::solve() {
    if (artificials == 0) {
        cout << "The problem is trivially feasible : no need for phase 1\n";
        simplex_method(SimplexTableau::Phase::Phase2);
    } else {
        cout << "The problem may be unfeasible : let's apply phase 1 / phase 2 method !\n";
        simplex_method(SimplexTableau::Phase::Phase1);
    }
}

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);

    string pivotName = "maxcoeff";
    int verboseLevel = 0;
    string name;

    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("input-file", po::value<string>(&name), "input_file")
        ("verbose", po::value<int>(&verboseLevel)->default_value(0),
        "Set verbose :\n"
        "\t0 : for only pivots (by default)\n"
        "\t1 : for pivots and optimal solution\n"
        "\t2 : for pivots, optimal solution and program\n"
        "\t3 : for pivots, optimal solution, program and tableau\n")
        ("pivot", po::value<string>(&pivotName)->default_value("maxcoeff"),
        "Set pivot type :\n"
        "\tmaxcoeff : for maxcoeff, Dantzig's rule (by default)\n"
        "\tsteepest : for steepest edge rule\n"
        "\tbland : for Bland's rule\n")
        ("no-cycles", "auto-detect cycles and leave them by Bland's rule")
    ;

    po::positional_options_description p;
    p.add("input-file", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    SimplexTableau::PivotType pivot;
    if (pivotName == "steepest")
        pivot = SimplexTableau::PivotType::Steepest;
    else if (pivotName == "bland")
        pivot = SimplexTableau::PivotType::Bland;
    else if (pivotName == "maxcoeff")
        pivot = SimplexTableau::PivotType::MaxCoeff;
    else
        throw runtime_error("Unknown pivot option\n");
    if (verboseLevel < 0 || verboseLevel > 3)
        throw runtime_error("Unknown verbose level\n");

    SimplexTableau tableau(verboseLevel, pivot, vm.count("no-cycles"));
    tableau.read_from_file(name);
    auto start = chrono::steady_clock::now();
    tableau.solve();
    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "Overall computation time : " << chrono::duration <double> (diff).count() << "ms" << endl;
    return EXIT_SUCCESS;
}
