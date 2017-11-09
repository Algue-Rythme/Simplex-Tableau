#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>

#include <boost/multiprecision/gmp.hpp>
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

    SimplexTableau(bool _isVerbose, PivotType _pivot):
    opt(0), nbPivots(0), isVerbose(_isVerbose), pivot(_pivot) {}

    void read_from_file(string name) {
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

    void print_problem(){
        if (!isVerbose)
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

    void print_tableau() {
        if (!isVerbose)
            return ;
        ArrayXr t;
        t.resize(A.rows() + 1, A.cols() + 1);
        t.topLeftCorner(1,  A.cols()) = c;
        t.topRightCorner(1, 1) = opt;
        t.bottomRightCorner(A.rows(), 1) = b;
        t.bottomLeftCorner(A.rows(), A.cols()) = A;
        cout << t << endl;
    }

    RowVectorXr::Index pivot_entering_DantzigRule() const {
        RowVectorXr::Index index;
        c.maxCoeff(&index);
        return index;
    }

    RowVectorXr::Index pivot_entering_Bland() const {
        for (int coeff = 0; coeff < c.cols(); ++coeff) {
            if (c(coeff) > 0)
                return coeff;
        }
        throw exception();
    }

    RowVectorXr::Index pivot_entering_SteepestEdgeRule() const {
        RowVectorXr norms = (c>0).select(A.matrix().colwise().squaredNorm(), 0);
        RowVectorXr steep = (norms!=0).select(c/norms, 0);
        RowVectorXr::Index index;
        steep.maxCoeff(&index);
        return index;
    }

    VectorXr::Index pivot_leaving_Bland(const VectorXr& A_col) const {
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

    enum class ProblemState {
        NotOptimum,
        Unbounded,
        Unfeasible,
        OptimumFound
    };

    void do_pivot(RowVectorXr::Index i_entering, VectorXr::Index i_leaving) {
        VectorXr goal = VectorXr::Zero(b.rows());
        goal(i_leaving) = 1;
        VectorXr coeffPivot = (goal - A.col(i_entering)) / A(i_leaving, i_entering);
        A += (coeffPivot.matrix() * A.row(i_leaving).matrix()).array().eval();
        b += (coeffPivot * b(i_leaving)).eval();
        Rational coeffCost = -c(i_entering);
        c += coeffCost * A.row(i_leaving);
        opt += coeffCost * b(i_leaving);
    }

    ProblemState simplex_iteration() {
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
        cout << "Entering x_" << (i_entering+1) << endl;
        VectorXr::Index i_leaving = pivot_leaving_Bland(A.col(i_entering));
        if (i_leaving == -1)
            return ProblemState::Unbounded;
        cout << "Leaving x_" << (basis[i_leaving]+1) << endl;
        do_pivot(i_entering, i_leaving);
        if (i_entering < n)
            affectation[i_entering] = i_leaving;
        if (basis[i_leaving] < n)
            affectation[basis[i_leaving]] = -1;
        basis[i_leaving] = i_entering;
        print_tableau();
        return ProblemState::NotOptimum;
    }

    enum class Phase {
        Phase1,
        Phase2
    };

    void simplex_method(Phase phase) {
        if (phase == Phase::Phase1)
            cout << "Phase 1\n";
        else
            cout << "Phase 2\n";
        if (isVerbose) {
            cout << "The initial tableau is:\n";
            print_tableau();
        }
        while (true) {
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
                        simplex_method(Phase::Phase2);
                        return ;
                    }
                }
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
                cout << "\nThe value of the objective for this solution is: " << -opt << "\n";
                break ;
            }
            nbPivots += 1;
        }
        cout << "The number of pivots is: " << nbPivots << "\n";
        cout << "Pivot rule used: ";
        switch (pivot) {
            case PivotType::Bland:
            cout << "Bland's rule\n";
            break;
            case PivotType::Steepest:
            cout << "Steepest Edge rule\n";
            break;
            default:
            cout << "Dantzig's rule\n";
            break;
        }
    }

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
    bool isVerbose;
    PivotType pivot;
};

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    if (argc <= 3) {
        cout << "Please provide an input file\n";
        return EXIT_FAILURE;
    }
    string name(argv[3]);
    string pivotName = string(argv[2]);
    SimplexTableau::PivotType pivot;
    if (pivotName == "-steepest")
        pivot = SimplexTableau::PivotType::Steepest;
    else if (pivotName == "-bland")
        pivot = SimplexTableau::PivotType::Bland;
    else if (pivotName == "-maxcoeff")
        pivot = SimplexTableau::PivotType::MaxCoeff;
    else
        throw exception();
    SimplexTableau tableau(string(argv[1]) == "-verbose", pivot);
    tableau.read_from_file(name);
    if (tableau.artificials == 0) {
        cout << "The problem is trivially feasible : no need for phase 1\n";
        tableau.simplex_method(SimplexTableau::Phase::Phase2);
    } else {
        cout << "Th problem may be unfeasible : let's apply phase 1 / phase 2 method !\n";
        tableau.simplex_method(SimplexTableau::Phase::Phase1);
    }
    return EXIT_SUCCESS;
}
