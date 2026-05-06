#pragma once

#include <vector>

namespace linalg {

struct SymmetricEigenSystem {
    int n = 0;
    std::vector<double> eigenvalues;
    std::vector<double> eigenvectors;

    double eigenvector(int row, int col) const;
};

bool has_lapacke();

SymmetricEigenSystem diagonalize_symmetric(std::vector<double> matrix, int n);

}
