#include <Eigen/Dense>
#include "sumc.h"

double sumL(const Eigen::MatrixXd &M){
    double tmp = 0;
    for (size_t i=0; i < M.size(); i++){
        tmp = tmp + *(M.data() + i);
    }
    return tmp;
}

double sumE(const Eigen::MatrixXd &M){
    return M.sum();
}

Eigen::MatrixXd sumEMatrices(const Eigen::MatrixXd &M1, const Eigen::MatrixXd &M2){
    return M1 + M2;
}
