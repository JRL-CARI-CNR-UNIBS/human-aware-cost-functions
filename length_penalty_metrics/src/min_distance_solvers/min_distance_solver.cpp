/*
Copyright (c) 2019, Cesare Tonola University of Brescia c.tonola001@unibs.it
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

MinDistanceSolver::MinDistanceSolver(const rosdyn::ChainPtr &chain, const double& max_step_size):
  chain_(chain)
{
  setMaxStepSize(max_step_size);
  obstacles_positions_.resize(0,0);
}

MinDistanceSolver:: MinDistanceSolver(const rosdyn::ChainPtr &chain, const double& max_step_size,
                                      const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions):
  chain_(chain)
{
  setMaxStepSize(max_step_size);
  setObstaclesPositions(obstacles_positions);
}

void MinDistanceSolver::setMaxStepSize(const double& max_step_size)
{
  max_step_size_ = max_step_size;
  if(max_step_size_<=0)
  {
    ROS_ERROR("max_step_size must be positive, set equal to 0.05");
    max_step_size_ = 0.05;
  }
}

DistancePtr MinDistanceSolver::computeMinDistance(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2)
{
  DistancePtr res = std::make_shared<Distance>();
  if(obstacles_positions_.cols() == 0)
  {
    res->distance_ = std::numeric_limits<double>::infinity(); //set infinity when there are no obstacles
    return res;
  }

  unsigned int iter = std::ceil((q2-q1).norm()/max_step_size_);

  Eigen::VectorXd q = q1;
  Eigen::VectorXd delta_q = (q2-q1)/iter;

  std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> poi_poses_in_base;

  double distance, min_distance;
  unsigned int poi_idx, obj_idx;
  Eigen::VectorXd robot_configuration;
  Eigen::Vector3d distance_vector, min_distance_vector;

  min_distance = std::numeric_limits<double>::infinity();

  for(unsigned int i=0;i<iter+1;i++)  //try it in parallel --> TO DO
  {
    poi_poses_in_base = chain_->getTransformations(q);

    for (Eigen::Index i_obj=0;i_obj<obstacles_positions_.cols();i_obj++)
    {
      for (size_t i_poi=0;i_poi<poi_poses_in_base.size();i_poi++)
      {
        distance_vector = obstacles_positions_.col(i_obj)-poi_poses_in_base.at(i_poi).translation(); //in base
        distance = distance_vector.norm();

        if(distance<min_distance)
        {
          min_distance = distance;
          min_distance_vector = distance_vector;

          robot_configuration = q;

          poi_idx = i_poi;
          obj_idx = i_obj;
        }
      }
    }

    q = q+delta_q;
  }

  res->distance_            = min_distance       ;
  res->distance_vector_     = min_distance_vector;
  res->robot_configuration_ = robot_configuration;
  res->robot_poi_           = poi_idx            ;
  res->object_              = obj_idx            ;

  return res;
}

Eigen::VectorXd MinDistanceSolver::getMaxJointsVelocity(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2)
{

}

}
