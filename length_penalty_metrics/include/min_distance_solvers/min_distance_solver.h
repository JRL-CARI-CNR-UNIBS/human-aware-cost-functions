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

#include <rosdyn_core/primitives.h>
#include <min_distance_solvers/util.h>

/* The purpose of the class is the computation of the minimum distance between the human and the robot given a connection as input */

namespace ssm15066_estimator
{
class MinDistanceSolver;
typedef std::shared_ptr<MinDistanceSolver> MinDistanceSolverPtr;

/**
 * @brief The MinDistanceSolver class computes the minimum distance between a set of objects and a set of robot's points of interests (poi)
 */
class MinDistanceSolver
{
protected:
  /**
  * @brief max_step_size_: max step between consecutive points along a connection for which the distance robot-obstacles is measured.
  */
  double max_step_size_;

  /**
   * @brief chain_: robot chain.
   */
  rosdyn::ChainPtr chain_;

  /**
   * @brief obstacles_positions_: x,y,z (rows) of obstacles (cols). Number of cols depends on the number of obstacles present in the scene.
   */
  Eigen::Matrix<double,3,Eigen::Dynamic> obstacles_positions_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  MinDistanceSolver(const rosdyn::ChainPtr& chain, const double& max_step_size);
  MinDistanceSolver(const rosdyn::ChainPtr &chain, const double& max_step_size,
                    const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions);

  rosdyn::ChainPtr getChain(){return chain_;}

  void setMaxStepSize(const double& max_step_size);
  void setChain(const rosdyn::ChainPtr& chain){chain_ = chain;}
  void setObstaclesPositions(const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions)
  {obstacles_positions_ = obstacles_positions;}

  /**
   * @brief computeMinDistance: computes the minimum distance between the robot's points of interests (poi) and the obstacles present in the scene.
   * It discretizes the connection (q1,q2) into more points and the distance between robot's poi and obstacles is computed for each one.
   * The maximum distance between consecutive points is equal to max_step_size_.
   * @param q1: first configuration of the connection.
   * @param q2: last configuration of the connection.
   * @return a DistancePtr object containing the minimum distance information.
   */
  virtual DistancePtr computeMinDistance(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2);
};

}
