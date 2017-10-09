#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // Estimation vector size should not be zero
  // Estimation vector size should equal the groudtruth vector size
  if(estimations.size() != ground_truth.size() || estimations.size() == 0){
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  // Accumulate squared residuals
  for(unsigned int i = 0; i < estimations.size(); i++) {
    VectorXd residual = estimations[i] - ground_truth[i];

    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // Take the mean
  rmse = rmse / estimations.size();

  // Squared root of rmse
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3, 4);

  // State params
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // Pre-compute some terms for use in calculating jacobian
  float positionsSquaredAdded = px*px+py*py;
  float sqrtPositionsSquaredAdded = sqrt(positionsSquaredAdded);
  float multipliedPositionsSquaresAndSqrt = positionsSquaredAdded * sqrtPositionsSquaredAdded;

  // Don't divide by zero
  if (fabs(positionsSquaredAdded) < 0.0001) {
    cout << "Jacobian calculation error - Division by zero is not allowed" << endl;
    return Hj;
  }

  // Compute Jacobian Matrix
  Hj << (px/sqrtPositionsSquaredAdded), (py/sqrtPositionsSquaredAdded), 0, 0,
        -(py/positionsSquaredAdded), (px/positionsSquaredAdded), 0, 0,
        py*(vx*py - vy*px)/multipliedPositionsSquaresAndSqrt, px*(px*vy - py*vx)/multipliedPositionsSquaresAndSqrt, px/sqrtPositionsSquaredAdded, py/sqrtPositionsSquaredAdded;

  return Hj;
}
