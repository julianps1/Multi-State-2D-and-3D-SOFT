#include "linalg.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <lapacke.h>

namespace linalg {

double SymmetricEigenSystem::eigenvector(int row, int col) const
{
    return eigenvectors[row * n + col];
}

bool has_lapacke()
{
    return true;
}

SymmetricEigenSystem diagonalize_symmetric(std::vector<double> matrix, int n)
{
    if (n <= 0) {
        throw std::invalid_argument("Matrix dimension must be positive");
    }

    if (matrix.size() != static_cast<std::size_t>(n * n)) {
        throw std::invalid_argument("Matrix size does not match dimension");
    }

    SymmetricEigenSystem result;
    result.n = n;
    result.eigenvalues.resize(n);

    const int info = LAPACKE_dsyev(
        LAPACK_ROW_MAJOR,
        'V',
        'U',
        n,
        matrix.data(),
        n,
        result.eigenvalues.data());

    if (info != 0) {
        throw std::runtime_error(
            "LAPACKE_dsyev failed with info = " + std::to_string(info));
    }

    result.eigenvectors = std::move(matrix);
    return result;
}

}
