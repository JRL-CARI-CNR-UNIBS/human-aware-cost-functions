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

#include <ros/ros.h>
#include <eigen3/Eigen/Core>
#include <rosdyn_core/primitives.h>

namespace ssm15066_estimator
{
class SSM15066Estimator;
typedef std::shared_ptr<SSM15066Estimator> SSM15066EstimatorPtr;

/**
  * @brief The SSM15066Estimator class is a template for safety related velocity scaling factor (SSM ISO-15066) estimator.
  * The goal is to compute an approximation of the average scaling factor that the robot will experiment moving from a
  * configuration q1 to a configuration q2, given the obstacles positions (e.g., human's head, arms and torso positions)
  */
class SSM15066Estimator
{
protected:

  /**
   * @brief chain_ is the robot's structure.
   */
  rosdyn::ChainPtr chain_;

  /**
   * @brief obstacles_positions_: matrix containing obstacles positions. x,y,z (rows) of obstacles (cols). Number of cols depends on the number of obstacles present in the scene.
   */
  Eigen::Matrix <double,3,Eigen::Dynamic> obstacles_positions_;

  /**
   * @brief inv_max_speed_ is the inverse of the max joints speed.
   */
  Eigen::VectorXd inv_max_speed_;

  /**
  * @brief max_step_size_: max step between consecutive points along a connection for which the distance robot-obstacles is measured.
  */
  double max_step_size_;

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
  double min_distance_;

  /**
   * @brief max_cart_acc_ is the maximum cartesian robot acceleration (m/s^2)
   */
  double max_cart_acc_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  SSM15066Estimator(const rosdyn::ChainPtr &chain, const double& max_step_size=0.05);
  SSM15066Estimator(const rosdyn::ChainPtr &chain, const double& max_step_size,
                    const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions);

  SSM15066Estimator(const urdf::ModelInterfaceSharedPtr &model, const std::string& base_frame, const std::string& tool_frame, const double& max_step_size=0.05);


  void updateMembers();
  void setMaxStepSize(const double& max_step_size);
  void setReactionTime(const double& t_r){t_r_ = t_r;}
  void setMaxCartAcc(const double& max_cart_acc){max_cart_acc_ = max_cart_acc;}
  void setMinDistance(const double& min_distance){min_distance_ = min_distance;}
  void setObstaclesPositions(const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions){obstacles_positions_ = obstacles_positions;}

  /**
   * @brief getObstaclePosition return the obstacles positions matrix
   * @return the matrix
   */
  Eigen::Matrix<double,3,Eigen::Dynamic> getObstaclesPositions(){return obstacles_positions_;}

  /**
   * @brief addObstaclePosition adds an obstacle to the already existing obstacle matrix.
   * @param obstacle_position is the new obstacle position vector (x,y,z)
   */
  void addObstaclePosition(const Eigen::Vector3d& obstacle_position);

  /**
   * @brief computeWorstCaseScalingFactor computes an approximation of the average scaling factor the robot will experience travelling from
   * q1 to q2, according to SSM ISO-15066. The maximum robot joints' velocities are considered for this computation.
   * @param q1.
   * @param q2.
   * @return 0 if a point q in (q1,q2) is associated with 0 scaling factor, the average scaling factor otherwise.
   */
  virtual double computeScalingFactor(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2);

  /**
   * @brief clone creates a copy of the object
   * @return the cloned object
   */
  virtual SSM15066EstimatorPtr clone();
};

}
