#pragma once
/*
Copyright (c) 2023, Cesare Tonola University of Brescia c.tonola001@unibs.it
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <min_distance_solvers/min_distance_solver.h>

/* The purpose of the class is the definition of a template for safety related velocity scaling (SSM ISO-15066) factor estimator.
 * The goal is to compute an approximation of the scaling factor that the robot will experiment moving from a
 * configuration q1 to a configuration q2, given the obstacles positions (e.g., human's head, arms and torso positions)
*/

namespace ssm15066_estimator
{
class SSM15066Estimator;
typedef std::shared_ptr<SSM15066Estimator> SSM15066EstimatorPtr;

/**
  * @brief The SSM15066Estimator class is a template for safety related velocity scaling factor (SSM ISO-15066) estimator.
  * The goal is to compute an approximation of the scaling factor lambda that the robot will experiment moving from a
  * configuration q1 to a configuration q2, given the obstacles positions (e.g., human's head, arms and torso positions)
  */
class SSM15066Estimator
{
private:
  /**
   * @brief distance_calculator_ is the object in charge of computing the minimum distance between
   * a set of robot's points of interest (poi) and a set of objects.
   */
  MinDistanceSolverPtr distance_calculator_;

//  /**
//   * @brief v_h_ is the cartesian human velocity towards the robot.  NOT CONSIDERED FOR NOW
//   */
//  double v_h_;

  /**
   * @brief t_r_ is the reaction time of the system.
   */
  double t_r_ ;

  /**
   * @brief term1_ stores the constant part of the SSM15066 equation (human velocity not considered).
   */
  double term1_;

  /**
   * @brief a_t_r_ is the term a*Tr of the SSM15066 equation.
   */
  double a_t_r_;

  /**
   * @brief min_distance_ is the minimum human-robot allowed distance.
   */
  double min_distance_ ;

  /**
   * @brief max_cart_acc_ is the maximum cartesian robot acceleration (m/s^2)
   */
  double max_cart_acc_ ;

  /**
   * @brief self_distance_ ??????????????????????????????????????????????????????????????
   */
  double self_distance_;

  /**
   * @brief chain_ is the robot's structure.
   */
  rosdyn::ChainPtr chain_;

  /**
   * @brief inv_max_speed_ is the inverse of the max joints speed.
   */
  Eigen::VectorXd inv_max_speed_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  SSM15066Estimator(const MinDistanceSolverPtr &distance_calculator);

  void setObstaclesPositions(const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions)
  {distance_calculator_->setObstaclesPositions(obstacles_positions);}

  /**
   * @brief computeScalingFactor computes an approximation of the scaling factor the robot will experience travelling from
   * q1 to q2, according to SSM ISO-15066. The maximum robot joints' velocities are considered for this computation.
   * @param q1.
   * @param q2.
   * @return the scaling factor.
   */
  double computeScalingFactor(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2);
};

}