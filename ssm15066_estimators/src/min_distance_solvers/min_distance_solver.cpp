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

namespace ssm15066_estimator
{

MinDistanceSolver::MinDistanceSolver(const rosdyn::ChainPtr &chain):
  chain_(chain)
{
  obstacles_positions_.resize(3,0);

  links_names_ = chain_->getLinksName();
  poi_names_ = links_names_;
}

MinDistanceSolver:: MinDistanceSolver(const rosdyn::ChainPtr &chain,
                                      const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions):
  chain_(chain)
{
  setObstaclesPositions(obstacles_positions);

  links_names_ = chain_->getLinksName();
  poi_names_ = links_names_;
}

void MinDistanceSolver::addObstaclePosition(const Eigen::Vector3d& obstacle_position)
{
  obstacles_positions_.conservativeResize(Eigen::NoChange, obstacles_positions_.cols()+1);

  if(obstacle_position.rows()>1) // column vector
    obstacles_positions_.col(obstacles_positions_.cols()-1) = obstacle_position;
  else
    obstacles_positions_.col(obstacles_positions_.cols()-1) = obstacle_position.transpose();  // make it a column vector
}

DistancePtr MinDistanceSolver::computeMinDistance(const Eigen::VectorXd& q)
{
  DistancePtr res = std::make_shared<Distance>();
  if(obstacles_positions_.cols() == 0)
  {
    res->distance_ = std::numeric_limits<double>::infinity(); //set infinity when there are no obstacles
    return res;
  }

  double distance;
  unsigned int poi_idx, obs_idx;
  Eigen::Vector3d distance_vector, min_distance_vector, i_poi_fk, poi_fk;

  double min_distance = std::numeric_limits<double>::infinity();
  std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> poi_poses_in_base = chain_->getTransformations(q);

  for (Eigen::Index i_obs=0;i_obs<obstacles_positions_.cols();i_obs++)
  {
    for (size_t i_poi=0;i_poi<poi_poses_in_base.size();i_poi++)
    {
      //consider only links inside the poi_names_ list
      if(std::find(poi_names_.begin(),poi_names_.end(),links_names_[i_poi])>=poi_names_.end())
        continue;

      i_poi_fk = poi_poses_in_base.at(i_poi).translation();
      distance_vector = obstacles_positions_.col(i_obs)-i_poi_fk; //in base
      distance = distance_vector.norm();

      if(distance<min_distance)
      {
        min_distance = distance;
        min_distance_vector = distance_vector;
        poi_fk = i_poi_fk;

        poi_idx = i_poi;
        obs_idx = i_obs;
      }
    }
  }

  res->poi_fk_              = poi_fk             ;
  res->obstacle_            = obs_idx            ;
  res->distance_            = min_distance       ;
  res->robot_poi_           = poi_idx            ;
  res->distance_vector_     = min_distance_vector;
  res->robot_configuration_ = q                  ;

  return res;
}

MinDistanceSolverPtr MinDistanceSolver::clone()
{
  MinDistanceSolverPtr clone = std::make_shared<MinDistanceSolver>(chain_->clone(),obstacles_positions_);
  return clone;
}


}
