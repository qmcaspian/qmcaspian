#ifndef SUMC_H_INCLUDED
#define SUMC_H_INCLUDED

#include <Eigen/Dense>

double sumL(const Eigen::MatrixXd &M);
double sumE(const Eigen::MatrixXd &M);
Eigen::MatrixXd sumEMatrices(const Eigen::MatrixXd &M1, const Eigen::MatrixXd &M2);

#endif
