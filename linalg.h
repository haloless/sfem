#pragma once

#include "sfem.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>


BEGIN_SFEM_NS;

using Eigen::VectorXd;
using Eigen::MatrixXd;

using Eigen::Vector2d;
using Eigen::Vector3d;


// Eigen Sparse Matrix, column-major
using SparseMatrix = Eigen::SparseMatrix<double>;
// colume-major
using SpMatColMajor = Eigen::SparseMatrix<double, Eigen::ColMajor>;
// row-major
using SpMatRowMajor = Eigen::SparseMatrix<double, Eigen::RowMajor>;

END_SFEM_NS;

