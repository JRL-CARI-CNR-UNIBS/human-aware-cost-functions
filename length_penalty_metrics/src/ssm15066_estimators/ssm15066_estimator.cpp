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

#include <ssm15066_estimators/ssm15066_estimator.h>

namespace ssm15066_estimator
{

SSM15066Estimator::SSM15066Estimator(const rosdyn::ChainPtr &chain, const double& max_step_size):
  chain_(chain)
{
  setMaxStepSize(max_step_size);

  inv_max_speed_ = chain_->getDQMax().cwiseInverse();

  min_distance_   = 0.15;
  max_cart_acc_   = 2.50;
  reaction_time_  = 0.15;
  human_velocity_ = 1.60; // from safety standards

  updateMembers();

  links_names_ = chain_->getLinksName();
  poi_names_ = links_names_;

  verbose_ = 0;
}

SSM15066Estimator::SSM15066Estimator(const rosdyn::ChainPtr &chain, const double &max_step_size, const Eigen::Matrix<double,3,Eigen::Dynamic> &obstacles_positions):
  chain_(chain), obstacles_positions_(obstacles_positions)
{
  setMaxStepSize(max_step_size);

  inv_max_speed_ = chain_->getDQMax().cwiseInverse();

  min_distance_   = 0.15;
  max_cart_acc_   = 2.50;
  reaction_time_  = 0.15;
  human_velocity_ = 1.60; // from safety standards

  updateMembers();

  links_names_ = chain_->getLinksName();
  poi_names_ = links_names_;

  verbose_ = 0;
}

void SSM15066Estimator::updateMembers()
{
  double a_tr = max_cart_acc_*reaction_time_;
  term1_ = std::pow(human_velocity_,2.0)+std::pow(a_tr,2.0)-2.0*max_cart_acc_*min_distance_;
  term2_ = -a_tr-human_velocity_;
}

void SSM15066Estimator::setMaxStepSize(const double& max_step_size)
{
  max_step_size_ = max_step_size;
  if(max_step_size_<=0)
  {
    ROS_ERROR("max_step_size must be positive, set equal to 0.05");
    max_step_size_ = 0.05;
  }
}

void SSM15066Estimator::addObstaclePosition(const Eigen::Vector3d& obstacle_position)
{
  obstacles_positions_.conservativeResize(Eigen::NoChange, obstacles_positions_.cols()+1);

  if(obstacle_position.rows()>1) // column vector
    obstacles_positions_.col(obstacles_positions_.cols()-1) = obstacle_position;
  else
    obstacles_positions_.col(obstacles_positions_.cols()-1) = obstacle_position.transpose();  // make it a column vector
}
}
