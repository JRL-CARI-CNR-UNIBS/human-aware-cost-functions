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
  * @brief The SSM15066Estimator class is a template for safety related velocity scaling factor (SSM ISO/TS-15066) estimator.
  * The goal is to compute an approximation of the average scaling factor that the robot will experiment moving from a
  * configuration q1 to a configuration q2, given the obstacles positions (e.g., human's head, arms and torso positions).
  * The scaling factor is a value from 1.0 to infinite and it is computed as v_r/v_safety.
  */
class SSM15066Estimator
{
protected:

  /**
   * @brief chain_ is the robot's structure.
   */
  rosdyn::ChainPtr chain_;

  /**
   * @brief poi_names_ is a vector containing the names of the points of interest of the robot structure.
   */
  std::vector<std::string> poi_names_;

  /**
   * @brief links_names_ is a vector of the names of all the links of the robot.
   */
  std::vector<std::string> links_names_;

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

  /**
     * @brief human_velocity_ is the cartesian human velocity towards the robot.
     */
  double human_velocity_;

  /**
   * @brief reaction_time_ is the reaction time of the robotic system.
   */
  double reaction_time_;

  /**
   * @brief max_cart_acc_ is the maximum cartesian robot acceleration (m/s^2)
   */
  double max_cart_acc_;

  /**
   * @brief min_distance_ is the minimum human-robot allowed distance.
   */
  double min_distance_;

  /**
   * @brief term1_ stores a constant part of the SSM15066 equation.
   */
  double term1_;

  /**
   * @brief term2_ stores a constant part of the SSM15066 equation.
   */
  double term2_;

  /**
   * @brief verbose_ defines thelevel of  verbosity
   */
  unsigned int verbose_;

  /**
   * @brief safeVelocity applies the SSM equation to compute the maximum robot cartesian velocity given the minimum human-robot distance as input
   * @param distance is the minimum human-robot distance
   * @return the safe robot velocity according to ISO/TS 15066
   */
  double safeVelocity(const double& distance)
  {
    assert(term1_+2.0*max_cart_acc_*distance>=0);
    assert(term1_ == std::pow(human_velocity_,2.0)+std::pow(max_cart_acc_*reaction_time_,2.0)-2.0*max_cart_acc_*min_distance_);
    assert(term2_ == -max_cart_acc_*reaction_time_-human_velocity_);

    return std::max(std::sqrt(term1_+2.0*max_cart_acc_*distance)+term2_,0.0);
  }

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  SSM15066Estimator(const rosdyn::ChainPtr &chain, const double& max_step_size=0.05);
  SSM15066Estimator(const rosdyn::ChainPtr &chain, const double& max_step_size,
                    const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions);

  /**
   * @brief Change a class member using the following functions. By defaults, term1_ and term2_ are updated
   * everytime you modify a class member. If you need to change multiple members and you don't want to call
   * updateMembers() everytime set update = false in each function and call updateMembers() at the end.
   */
  void updateMembers();
  void setMaxCartAcc(const double& max_cart_acc, const bool update = true)
  {
    max_cart_acc_ = max_cart_acc;
    if(update)
      updateMembers();
  }
  void setMinDistance(const double& min_distance, const bool update = true)
  {
    min_distance_ = min_distance;
    if(update)
      updateMembers();
  }
  void setReactionTime(const double& reaction_time, const bool update = true)
  {
    reaction_time_ = reaction_time;
    if(update)
      updateMembers();
  }
  void setHumanVelocity(const double& human_velocity, const bool update = true)
  {
    human_velocity_ = human_velocity;
    if(update)
      updateMembers();
  }
  void setPoiNames(const std::vector<std::string> poi_names)
  {
    if(poi_names.empty())
    {
      ROS_ERROR("Poi names void");
      return;
    }
    poi_names_ = poi_names;
  }

  void setMaxStepSize(const double& max_step_size);
  void setVerbose(const unsigned int& verbose){verbose_ = verbose;}

  /**
    Getters
   */
  rosdyn::ChainPtr         getChain        (){return chain_         ;}
  std::vector<std::string> getPoiNames     (){return poi_names_     ;}
  double                   getMaxCartAcc   (){return max_cart_acc_  ;}
  double                   getMinDistance  (){return min_distance_  ;}
  double                   getMaxStepSize  (){return max_step_size_ ;}
  double                   getReactionTime (){return reaction_time_ ;}
  double                   getHumanVelocity(){return human_velocity_;}

  /**
   * @brief getObstaclePosition return the obstacles positions matrix
   * @return the matrix
   */
  Eigen::Matrix<double,3,Eigen::Dynamic> getObstaclesPositions(){return obstacles_positions_;}

  /**
   * @brief setObstaclesPositions sets the matrix of obstacles locations
   * @param obstacles_positions is the matrix containing in the columns the location of each obstacle as x,y,z
   */
  virtual void setObstaclesPositions(const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions){obstacles_positions_ = obstacles_positions;}

  /**
   * @brief addObstaclePosition adds an obstacle to the already existing obstacles matrix.
   * @param obstacle_position is the new obstacle location vector (x,y,z)
   */
  virtual void addObstaclePosition(const Eigen::Vector3d& obstacle_position);

  /**
   * @brief clearObstaclePosition clears the matrix of obstacles locartions
   */
  void clearObstaclesPositions(){obstacles_positions_.resize(3,0);}

  /**
   * @brief computeWorstCaseScalingFactor computes an approximation of the average scaling factor the robot will experience travelling from
   * q1 to q2, according to SSM ISO-15066. The maximum robot joints' velocities are considered for this computation.
   * The scaling factor is computed as the average of the scaling factor at each qi between q1 and q2. The scaling factor of qi
   * is measurd as v_r(i)/v_safety(i).
   * @param q1.
   * @param q2.
   * @return the average scaling factor.
   */
  virtual double computeScalingFactor(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2) = 0;

  /**
   * @brief clone creates a copy of the object
   * @return the cloned object
   */
  virtual SSM15066EstimatorPtr clone() = 0;
};

}
