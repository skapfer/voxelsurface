#ifndef GAUSS_SOLVER_H_INCLUDED 
#define GAUSS_SOLVER_H_INCLUDED 
#include <cassert>
#include <algorithm>
#include <iomanip>

template <int NUM>
class GaussSolver {
public:
    void set_coeff (int i, int j, double val) {
        assert (i >= 0);
        assert (i < NUM);
        assert (j >= 0);
        assert (j < NUM);
        ar[i][j] = val;
    }

    void set_right (int i, double val) {
        assert (i >= 0);
        assert (i < NUM);
        ar[i][NUM] = val;
    }

    void solve () {
        for (int i = 0; i != NUM; ++i) {
            pivot (i);
            eliminate (i);
        }
        backsubstitute ();
    }

    double get_solution (int i) {
        assert (i >= 0);
        assert (i < NUM);
        return ar[i][NUM];
    }

    friend std::ostream &operator<< (std::ostream &os, const GaussSolver &this_) {
        for (int i = 0; i != NUM; ++i) {
            for (int j = 0; j != NUM; ++j) {
                os << std::setw (8) << this_.ar[i][j] << " ";
            }
            os << " | " << std::setw (8) << this_.ar[i][NUM] << "\n";
        }
        return os;
    }

private:
    void pivot (int i) {
        double pivv = fabs (ar[i][i]);
        int pivi = i;
        for (int ii = i+1; ii != NUM; ++ii) {
            if (pivv < fabs (ar[ii][i])) {
                pivv = fabs (ar[ii][i]);
                pivi = ii;
            }
        }

        for (int j = i; j != NUM+1; ++j)
            std::swap (ar[pivi][j], ar[i][j]);
        double d = ar[i][i];
        if (d == 0) {
            // system has many solutions, just pick one.
            ar[i][i] = 1.;
        } else {
            assert (fabs (d) > 0.);
            for (int j = i; j != NUM+1; ++j)
                ar[i][j] /= d;
        }
    }

    void eliminate (int i) {
        for (int ii = i+1; ii != NUM; ++ii) {
            for (int j = 0; j != i; ++j)
                assert (ar[ii][j] == 0.);
            double d = ar[ii][i];
            for (int j = i; j != NUM+1; ++j)
                ar[ii][j] -= d * ar[i][j];
            if (ar[ii][i] != 0)
                std::cerr << "\n" << ar[ii][i] << "\n";
            assert (ar[ii][i] == 0.);
        }
    }

    void backsubstitute () {
        for (int ii = NUM-1; ii != -1; --ii) {
            double ac = ar[ii][NUM];
            for (int j = ii+1; j != NUM; ++j) {
                ac -= ar[ii][j] * ar[j][NUM];
                ar[ii][j] = 0.;
            }
            ar[ii][NUM] = ac;
        }
    }

private:
    double ar[NUM][NUM+1];
};

#endif // GAUSS_SOLVER_H_INCLUDED
