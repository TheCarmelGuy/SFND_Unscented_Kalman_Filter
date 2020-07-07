#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() { 



  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0); 
  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.fill(0);




  // Process noise standard deviation longitudinal acceleration in m/s^2

  std_a_ = 10;
  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
  std_yawdd_ = 3;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  //generate weights
  weights_ = VectorXd(2*n_aug_ + 1);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  weights_.tail(2*n_aug_) = (0.5f/(lambda_ + n_aug_))*VectorXd::Ones(2*n_aug_);


  Xsig_pred_ = Eigen::MatrixXd(n_x_,2*n_aug_ + 1); 


}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  bool useLaserFrame = meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_;
  bool useRadarFrame = meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_;

  if(is_initialized_ == false)
  {
    if(useLaserFrame)
    {
    
      
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);

    //  std::cerr<<"Initial state"<<std::endl;
    //  std::cerr<<"===="<<std::endl;
    //  std::cerr<<x_<<std::endl;
    //  std::cerr<<"===="<<std::endl;
      // set initial cov based on lidar noise... remaining params set arbitarily large
      // since we can't use lidar to measure
      P_(0,0) = std_laspx_*std_laspx_;
      P_(1,1) = std_laspy_*std_laspy_;
      P_(2,2) = 1; 
      P_(3,3) = 1;
      P_(4,4) = 1;

      is_initialized_ = true; 
      time_us_  = meas_package.timestamp_;

    }
    else if(useRadarFrame)
    {
      

      double rho = meas_package.raw_measurements_(0);
      double theta = meas_package.raw_measurements_(1);
      double drho = meas_package.raw_measurements_(2);
      x_(0) = rho*cos(theta);
      x_(1) = rho*sin(theta);
      x_(2) = drho;
      x_(3) = theta;
      x_(4) = drho/rho; 


      P_(0,0) = 3;
      P_(1,1) = 3;
      P_(2,2) = 3; 
      P_(3,3) = 3;
      P_(4,4) = 3;

 
      is_initialized_ = true; 
      time_us_  = meas_package.timestamp_;

     
    }
  }
  

  if(is_initialized_ == true)
  {
   
    double delta_t =  (meas_package.timestamp_ - this->time_us_)/1000000.0f;

    if(useLaserFrame)
    {

      std::cerr<<"INFO: Ukf Processing Lidar Frame"<<std::endl;

      Prediction(delta_t);
     
     UpdateLidar(meas_package);
      this->time_us_ = meas_package.timestamp_;
     
    } 
    else if(useRadarFrame)
    {
      
      std::cerr<<"INFO: Ukf Processing Radar Frame"<<std::endl; 

      Prediction(delta_t);
      //std::cerr<<"State after prediction"<<std::endl;
      //std::cerr<<P_<<std::endl;
       
      UpdateRadar(meas_package);
  
      //std::cerr<<"State after Update"<<std::endl;
      //std::cerr<<P_<<std::endl;
     this->time_us_ = meas_package.timestamp_;

    }
    else 
    {
   ///   std::cerr<<"ERROR: UKF - Sensor type not currently supported"<<std::endl;
    }

  }
}

void UKF::GenerateSigmaPoints(MatrixXd&  _Xsig_aug,VectorXd& _x_aug)
{

   
    _x_aug.head(this->n_x_) = this->x_;
   
    //create augmented motion model cov
    Eigen::MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
    P_aug.topLeftCorner(n_x_,n_x_) = P_;
    P_aug(n_x_,n_x_) =std_a_*std_a_;
    P_aug(n_x_ + 1,n_x_ + 1) =std_yawdd_*std_yawdd_;

  
    // generate Sigma points

    _Xsig_aug.col(0) =  _x_aug;
    Eigen::MatrixXd P_sqrt = P_aug.llt().matrixL();
  

    for(int i = 0; i < n_aug_; i++)
    {

        _Xsig_aug.col(i+1) =  _x_aug + sqrt(lambda_+n_aug_)*P_sqrt.col(i);
        _Xsig_aug.col(i+1 + n_aug_) =  _x_aug - sqrt(lambda_+n_aug_)*P_sqrt.col(i);
         

    }


}


void UKF::PredictSigmaPoints(MatrixXd&  _Xsig_aug,double _delta_t)
{

  // predict sigma points
  for(int i = 0; i< _Xsig_aug.cols(); i++) //iterate through columns
  {
  
    VectorXd sigPoint = _Xsig_aug.col(i);
    double px_k = sigPoint(0); 
    double py_k = sigPoint(1);
    double v_k = sigPoint(2);
    double psi_k = sigPoint(3);
    double dpsi_k = sigPoint(4); 
    double va = sigPoint(5);
    double vpsi = sigPoint(6);
   


    double px = px_k;
    double py = py_k; 
    
    if(fabs(dpsi_k) < 0.001)
    {
        px += v_k*cos(psi_k)*_delta_t;
        py += v_k*sin(psi_k)*_delta_t;
    }
    else
    {
        px += (v_k/dpsi_k)*(sin(psi_k + dpsi_k*_delta_t) - sin(psi_k));
        py += (v_k/dpsi_k)*(-cos(psi_k + dpsi_k*_delta_t) + cos(psi_k));
    }
   

    //adding noise
    px += 0.5*_delta_t*_delta_t*cos(psi_k)*va;
    py += 0.5*_delta_t*_delta_t*sin(psi_k)*va;
    
    double v = v_k + 0 + _delta_t*va;
    double psi = psi_k + dpsi_k*_delta_t + 0.5*_delta_t*_delta_t*vpsi;
    double dpsi = dpsi_k + 0 + _delta_t*vpsi;
  
   
    
    Xsig_pred_.col(i) << px, py, v, psi, dpsi;
     

  }

}

void UKF::UpdateStateAndCovariance()
{


  //update mean state
  VectorXd x(n_x_);
  x.fill(0); 
  for(int i = 0; i < 2*n_aug_ + 1; i++)
  {
    x += weights_(i)*Xsig_pred_.col(i);
  }

  //update covariance of state;
  MatrixXd P(n_x_, n_x_);
  P.fill(0);
  for(int i = 0; i < 2*n_aug_ + 1; i++)
  {
   
    VectorXd x_diff =  Xsig_pred_.col(i) - x;
   

    // Bug here: Angle get's really big 
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;




    P += weights_(i)*x_diff*x_diff.transpose();
  }




  

 

  x_ = x;
  P_ = P;

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  
  Eigen::MatrixXd Xsig_aug(n_aug_, 2*n_aug_ + 1);
 
  Xsig_aug.fill(0);
  Eigen::VectorXd x_aug(n_aug_);
  x_aug.fill(0);


  Xsig_pred_.fill(0); //reset sig prediciton


  //generate sigPoint for X,P k|k
  GenerateSigmaPoints(Xsig_aug, x_aug);



  // predict sigPoints to yield X,P k+1 | k
  PredictSigmaPoints(Xsig_aug, delta_t);

  //Predict mean and covariance: X,P k+1 | k
  UpdateStateAndCovariance();

}

void UKF::PredictLidar(VectorXd& _z_exp, MatrixXd& _S, MatrixXd& _Tc)
{

  MatrixXd Zsig_pred(2,2*n_aug_ + 1);

  //calculate zsig
  for(int i = 0; i<Zsig_pred.cols(); i ++)
  {
    Zsig_pred.col(i) << Xsig_pred_(0,i), Xsig_pred_(1,i);
  }


  VectorXd z_pred(2);
  z_pred.fill(0);
  for(int i = 0; i<Zsig_pred.cols(); i ++)
  {
    z_pred += weights_(i)*Zsig_pred.col(i);

  }

  MatrixXd S(2,2);
  S.fill(0);
  for(int i = 0; i<Zsig_pred.cols(); i ++)
  {

    S += weights_(i)*(Zsig_pred.col(i) - z_pred)*(Zsig_pred.col(i)-z_pred).transpose();
  }

  MatrixXd R(2,2);
  R << std_laspx_*std_laspx_, 0,
    0, std_laspy_*std_laspy_;

  S += R; //add sensor Noise

  //Compute weighted cross correlation bewteen Zsig and x sig

  MatrixXd Tc(n_x_,2); 
  Tc.fill(0);
  for(int i = 0; i< Zsig_pred.cols(); i++)
  {
    Tc += weights_(i)*((Xsig_pred_.col(i) - x_)*(Zsig_pred.col(i) - z_pred).transpose());
  }


  _S = S;
  _z_exp = z_pred;
  _Tc = Tc;

}


void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */


  int z_dim = 2;

  VectorXd z_exp(z_dim);
  MatrixXd S(z_dim,z_dim);
  MatrixXd Tc(n_x_,z_dim);

  //predict Lidar w/ sigma points
  PredictLidar(z_exp, S, Tc);
  

  MatrixXd K = Tc*S.inverse();
  x_ = x_ + K*(meas_package.raw_measurements_ - z_exp);
  P_ = P_ - K*S*K.transpose();


  double NIS = (meas_package.raw_measurements_ - z_exp).transpose()*S.inverse()*(meas_package.raw_measurements_ - z_exp);

  std::cout<<"NIS :"<<NIS<<" "<<meas_package.LASER<<std::endl;

}


void UKF::PredictRadar(VectorXd& _z_exp, MatrixXd& _S, MatrixXd& _Tc)
{

  int z_dim  = 3;
  MatrixXd Zsig_pred(z_dim,2*n_aug_ + 1);
  Zsig_pred.fill(0);

  //calculate zsig
  for(int i = 0; i<Zsig_pred.cols(); i ++)
  {
    VectorXd x_pred = Xsig_pred_.col(i);
    double rho = sqrt(x_pred(0)*x_pred(0) + x_pred(1)*x_pred(1));
 
    double p_x =x_pred(0);
    double p_y= x_pred(1);
    double v = x_pred(2);
    double yaw = x_pred(3);
   
    double v_x =  cos(yaw)*v; 
    double v_y = sin(yaw)*v; 

     
    double theta = atan2(p_y,p_x);
//   double theta = atan(x_pred(1)/x_pred(0));
//    std::cout<<"atan: "<<theta<<std::endl;
    double drho = (p_x*v_x + p_y*v_y)/sqrt(p_x*p_x + p_y*p_y);


    Zsig_pred.col(i) << rho,theta, drho; 
 
  }


  VectorXd z_pred(z_dim);
  z_pred.fill(0);
  for(int i = 0; i<Zsig_pred.cols(); i ++)
  {
    z_pred += weights_(i)*Zsig_pred.col(i);

  }

  MatrixXd S(z_dim,z_dim);
  S.fill(0);
  for(int i = 0; i<Zsig_pred.cols(); i ++)
  {

    VectorXd z_diff = Zsig_pred.col(i) - z_pred;
 
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S += weights_(i)*(z_diff*z_diff.transpose());
  }

  MatrixXd R(z_dim,z_dim);
  R << std_radr_*std_radr_, 0, 0, 
    0, std_radphi_*std_radphi_, 0,
    0,  0, std_radrd_*std_radrd_;

  S += R; //add sensor Noise

  //Compute weighted cross correlation bewteen Zsig and x sig

  MatrixXd Tc(n_x_,z_dim); 
  Tc.fill(0);
  for(int i = 0; i< Zsig_pred.cols(); i++)
  {

    VectorXd z_diff = Zsig_pred.col(i) - z_pred;
 
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
 
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

   
    Tc = Tc +  weights_(i)*(x_diff*z_diff.transpose());
  
  }



  _S = S;
  _z_exp = z_pred;
  _Tc = Tc;
}



void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  int z_dim = 3;

  VectorXd z_exp(z_dim);
  MatrixXd S(z_dim, z_dim);
  MatrixXd Tc(n_x_,z_dim);

  PredictRadar(z_exp, S, Tc); 

  MatrixXd K = Tc*S.inverse();

  //normalize angles
  VectorXd z_diff = meas_package.raw_measurements_ - z_exp;


  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;



  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();


  double NIS = (meas_package.raw_measurements_ - z_exp).transpose()*S.inverse()*(meas_package.raw_measurements_ - z_exp);  

    std::cout<<"NIS :"<<NIS<<" "<<meas_package.RADAR<<std::endl;



}
