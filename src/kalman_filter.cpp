#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  Estimate_x_P(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  //for radar use h(x')
  VectorXd y = z - h();

  //Normalize angle
  if(y(1) > PI_)
    y(1) -= 2*PI_;
  if(y(1) < -PI_)
    y(1) += 2*PI_;

  //estimate x and P
  Estimate_x_P(y);
}

void KalmanFilter::Estimate_x_P(const Eigen::VectorXd& y){
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

const VectorXd KalmanFilter::h(){
  VectorXd res = VectorXd(3);
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  //pre-compute x^2 + y^2
  float c1 = px * px + py * py;

  //check division by zero
  if (fabs(c1) < 0.0001) {
    cout << __func__ << "Error - Division by Zero x=" << x_ << endl;
  }else{
    float c2 = sqrt(c1);
    float c3 = px * vx + py * vy;
    res(0) = c2;
    res(1) = atan2(py,px);
    res(2) = c3 / c2;
  }
  return res;
}
