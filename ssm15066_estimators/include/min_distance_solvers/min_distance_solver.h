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
   * @brief chain_: robot chain.
   */
  rosdyn::ChainPtr chain_;

  /**
   * @brief obstacles_positions_: x,y,z (rows) of obstacles (cols). Number of cols depends on the number of obstacles present in the scene.
   */
  Eigen::Matrix<double,3,Eigen::Dynamic> obstacles_positions_;

  /**
   * @brief poi_names_ is a vector containing the names of the points of interest of the robot structure.
   */
  std::vector<std::string> poi_names_;

  /**
   * @brief links_names_ is a vector of the names of all the links of the robot.
   */
  std::vector<std::string> links_names_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  MinDistanceSolver(const rosdyn::ChainPtr& chain);
  MinDistanceSolver(const rosdyn::ChainPtr &chain, const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions);

  rosdyn::ChainPtr getChain(){return chain_;}

  void setChain(const rosdyn::ChainPtr& chain){chain_ = chain;}
  void setPoiNames(const std::vector<std::string> poi_names){poi_names_ = poi_names;}
  void setObstaclesPositions(const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions){obstacles_positions_ = obstacles_positions;}

  Eigen::Matrix<double,3,Eigen::Dynamic> getObstaclesPositions(){return obstacles_positions_;}

  /**
   * @brief addObstaclePosition adds an obstacle to the already existing obstacle matrix.
   * @param obstacle_position is the new obstacle position vector (x,y,z)
   */
  void addObstaclePosition(const Eigen::Vector3d& obstacle_position);

  /**
   * @brief clearObstaclePosition clears the matrix of obstacles locations
   */
  void clearObstaclesPositions(){obstacles_positions_.resize(3,0);}

  /**
   * @brief computeMinDistance: computes the minimum distance between the robot's points of interests (poi) and the obstacles present in the scene.
   * @param q: robot configuration
   * @return a DistancePtr object containing the minimum distance information.
   */
  virtual DistancePtr computeMinDistance(const Eigen::VectorXd& q);

  /**
   * @brief clone gives a cloned and indipendent copy of the object
   * @return the cloned object
   */
  virtual MinDistanceSolverPtr clone();
};

}
