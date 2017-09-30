#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

const float  PI_F=3.14159265358979f;

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
  EstimatexP(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd y = z - h(x_); //for EKF use h(x')
  EstimatexP(y);
}

void KalmanFilter::EstimatexP(const Eigen::VectorXd& y){
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

const VectorXd KalmanFilter::h(const VectorXd& x){
  VectorXd res = VectorXd(3);
  float px = x(0);
  float py = x(1);
  float vx = x(2);
  float vy = x(3);

  //pre-compute a set of terms to avoid repeated calculation
  float c1 = px * px + py * py;

  //check division by zero
  if (fabs(c1) < 0.0001) {
    cout << __func__ << "Error - Division by Zero x=" << x << endl;
  }else{
    float c2 = sqrt(c1);
    float c3 = px * vx + py * vy;
    res(0) = c2;
    res(1) = atan2f(py, px);
    res(2) = c3 / c2;
  }
  return res;
}
