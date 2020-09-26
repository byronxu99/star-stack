/* --------------------------------------------------------------------------
 * File:    matrix.h
 * Author:  Michael Gharbi <gharbi@mit.edu>
 * Created: 2015-10-17
 * --------------------------------------------------------------------------
 *
 *
 *
 * ------------------------------------------------------------------------*/

#pragma once

#include <../Eigen/Dense>
#include <../Eigen/Sparse>

typedef Eigen::MatrixXf Matrix;
typedef Eigen::VectorXf Vector;
typedef Eigen::Vector2f Vec2f;
typedef Eigen::Vector3f Vec3f;

typedef Eigen::SparseMatrix<float> SparseMatrix;
typedef Eigen::Triplet<float> Triplet;

