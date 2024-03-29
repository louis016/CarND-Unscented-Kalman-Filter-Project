#include "tools.h"

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  VectorXd residual;

  // accumulate the residuals
  for (int i = 0; i < estimations.size(); ++i) {
     residual = estimations[i] - ground_truth[i];
     residual = residual.array() * residual.array();
     rmse += residual; 
  }
  
  // caculate the average
  rmse = rmse / estimations.size();
  
  // caculate the rmse
  rmse = rmse.array().sqrt();
  return rmse;
  
}