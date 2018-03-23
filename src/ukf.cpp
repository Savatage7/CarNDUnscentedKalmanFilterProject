#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
	is_initialized_=false;

	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 2;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.3;

	//DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
	//DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

	n_x_ = 5;
	n_aug_ = 7;
	lambda_ = 3 - n_aug_;

	// the Count for NIS for Radar >7.8
	NIS_radar_G78 = 0.0;
	// the Count for NIS for Lidar >7.8
	NIS_laser_G59 = 0.0;

	NIS_LaserCalulation = 0;
	NIS_RadarCalulation = 0;

	/*
	Complete the initialization. See ukf.h for other member properties.
	Hint: one or more values initialized above might be wildly off...
	 */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	/**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
	 */

	if (!is_initialized_) {

		x_ << 1, 1, 1, 1, 0.1;

		P_ << 0.15,    0, 0, 0, 0,
				0, 0.15, 0, 0, 0,
				0,    0, 1, 0, 0,
				0,    0, 0, 1, 0,
				0,    0, 0, 0, 1;



		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

			float rho = meas_package.raw_measurements_[0]; // range
			float phi = meas_package.raw_measurements_[1]; // bearing
			float rho_dot = meas_package.raw_measurements_[2]; // velocity of rho

			float x = rho * cos(phi);
			float y = rho * sin(phi);
			if ( fabs(x) < 0.0001 ) {
				x = 0.0001;
			}
			if ( fabs(y) < 0.0001 ) {
				y = 0.0001;
			}
			x_ << x, y, rho_dot, phi, 0;
			cout << "INI Radar x_ = " << x_ << endl;

		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
			cout << "INI Laser x_ = " << x_ << endl;

		}
		time_us_ = meas_package.timestamp_;

		// done initializing, no need to predict or update
		is_initialized_ = true;
		cout << "Done initialized" << endl;

		return;

	}
	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds

	Prediction(dt);
	time_us_ = meas_package.timestamp_;

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		// Radar updates

		NIS_RadarCalulation++;
		UpdateRadar(meas_package);
	} else {
		// Laser updates
		NIS_LaserCalulation++;
		UpdateLidar(meas_package);

	}
	std::cout << "NIS_radar_G78 " << NIS_radar_G78<<"/"<<NIS_RadarCalulation<<std::endl;
	std::cout << "NIS_laser_G59 " << NIS_laser_G59<<"/"<<NIS_LaserCalulation<<std::endl;

	std::cout <<std::endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	/**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
	 */

	/*****************************************************************************
	 *  Prediction
	 ****************************************************************************/

	//create augmented mean vector
	VectorXd x_aug = VectorXd(7);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(7, 7);

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	/*****************************************************************************
	 *  Sigma Points
	 ****************************************************************************/
	//create augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5,5) = P_;
	P_aug(5,5) = std_a_*std_a_;
	P_aug(6,6) = std_yawdd_*std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0)  = x_aug;
	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
		Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
	}

	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

	//predict sigma points
	for (int i = 0; i< 2*n_aug_+1; i++)
	{
		//extract values for better readability
		double p_x = Xsig_aug(0,i);
		double p_y = Xsig_aug(1,i);
		double v = Xsig_aug(2,i);
		double yaw = Xsig_aug(3,i);
		double yawd = Xsig_aug(4,i);
		double nu_a = Xsig_aug(5,i);
		double nu_yawdd = Xsig_aug(6,i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd*delta_t;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a*delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd*delta_t;

		//write predicted sigma point into right column
		Xsig_pred_(0,i) = px_p;
		Xsig_pred_(1,i) = py_p;
		Xsig_pred_(2,i) = v_p;
		Xsig_pred_(3,i) = yaw_p;
		Xsig_pred_(4,i) = yawd_p;
	}


	//	VectorXd
	weights_ = VectorXd(2*n_aug_+1);


	double weight_0 = lambda_/(lambda_+n_aug_);
	weights_(0) = weight_0;

	for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights_
		double weight = 0.5/(n_aug_+lambda_);
		weights_(i) = weight;
	}

	//predicted state mean
	x_.fill(0.0);

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		x_ = x_+ weights_(i) * Xsig_pred_.col(i);
	}

	//predicted state covariance matrix
	P_.fill(0.0);

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization

		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

		P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
	}


	//print result
	std::cout << "Predicted state" << std::endl;
	std::cout << x_ << std::endl;
	std::cout << "Predicted covariance matrix" << std::endl;
	std::cout << P_ << std::endl;


}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	/**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
	 */

	MatrixXd H_ = MatrixXd(2, 5);

	H_ << 1, 0, 0, 0, 0,
			0, 1, 0, 0, 0;

	MatrixXd R_ = MatrixXd(2, 2);
	R_ << 0.0225, 0,
			0, 0.0225;

	VectorXd z=meas_package.raw_measurements_;
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;


	if((y.transpose() * Si * y)>5.9){
		NIS_laser_G59++;
	}

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}



/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	/**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
	 */
	//set measurement dimension, radar can measure r, phi, and r_dot
	int n_z_ = 3;
	//transform sigma points into measurement space

	MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		// extract values for better readibility
		double p_x = Xsig_pred_(0,i);
		double p_y = Xsig_pred_(1,i);
		double v  = Xsig_pred_(2,i);
		double yaw = Xsig_pred_(3,i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;


		// measurement model
		Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
		Zsig(1,i) = atan2(p_y,p_x);                                 //phi
		Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_);
	z_pred.fill(0.0);
	for (int i=0; i < 2*n_aug_+1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//innovation covariance matrix S
	MatrixXd S = MatrixXd(n_z_,n_z_);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z_,n_z_);
	R <<    std_radr_*std_radr_, 0, 0,
			0, std_radphi_*std_radphi_, 0,
			0, 0,std_radrd_*std_radrd_;
	S = S + R;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z_);
	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	VectorXd z = VectorXd(n_z_);


	float rho = meas_package.raw_measurements_[0]; // range
	float phi = meas_package.raw_measurements_[1]; // bearing
	float rho_dot = meas_package.raw_measurements_[2]; // velocity of rho

	z << rho, phi, rho_dot;
	//residual
	VectorXd z_diff = z - z_pred;

	//angle normalization
	while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

	if((z_diff.transpose() * S.inverse() * z_diff)>7.8){
		NIS_radar_G78++;
	}

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();
}
