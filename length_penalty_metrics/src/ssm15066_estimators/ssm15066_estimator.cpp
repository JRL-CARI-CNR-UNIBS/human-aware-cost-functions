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

#include <ssm15066_estimators/ssm15066_estimator.h>

namespace ssm15066_estimator
{

SSM15066Estimator::SSM15066Estimator(const rosdyn::ChainPtr &chain, const double& max_step_size):
  chain_(chain)
{
  obstacles_positions_.resize(0,0);

  setMaxStepSize(max_step_size);

  inv_max_speed_ = chain_->getDQMax().cwiseInverse();

  self_distance_ = 0.00;
  min_distance_  = 0.30;
  max_cart_acc_  = 0.10;
  t_r_           = 0.15;

  updateMembers();
}

SSM15066Estimator::SSM15066Estimator(const rosdyn::ChainPtr &chain, const double &max_step_size, const Eigen::Matrix<double,3,Eigen::Dynamic> &obstacles_positions):
  chain_(chain), obstacles_positions_(obstacles_positions)
{
  setMaxStepSize(max_step_size);

  inv_max_speed_ = chain_->getDQMax().cwiseInverse();

  self_distance_  = 0.00;
  min_distance_   = 0.30;
  max_cart_acc_   = 0.10;
  t_r_            = 0.15;

  updateMembers();
}

void SSM15066Estimator::updateMembers()
{
  a_t_r_ = max_cart_acc_*t_r_;
  term1_=std::pow(a_t_r_,2)-2*max_cart_acc_*min_distance_;
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

double SSM15066Estimator::computeScalingFactor(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2)
{
  double scaling_factor = 1.0;

  if (obstacles_positions_.cols()==0)  //no obstacles in the scene
    return scaling_factor;

  Eigen::Vector3d distance_vector;
  double distance, tangential_speed, tmp_scaling_factor, v_safety;
  std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>>  poi_poses_in_base;
  std::vector< Eigen::Vector6d, Eigen::aligned_allocator<Eigen::Vector6d>> poi_twist_in_base;

  /* Compute the time of each joint to move from q1 to q2 at its maximum speed and consider the longest time */
  double slowest_joint_time = (inv_max_speed_.cwiseProduct(q2 - q1)).cwiseAbs().maxCoeff();

  /* The "slowest" joint will move at its highest speed while the other ones will
   * move at (t_i/slowest_joint_time)*max_speed_i, where slowest_joint_time >= t_i */
  Eigen::VectorXd dq = (q2-q1)/slowest_joint_time;

  unsigned int iter = std::ceil((q2-q1).norm()/max_step_size_);

  Eigen::VectorXd q = q1;
  Eigen::VectorXd delta_q = (q2-q1)/iter;

  for(unsigned int i=0;i<iter+1;i++)
  {
    poi_poses_in_base = chain_->getTransformations(q);
    poi_twist_in_base = chain_->getTwist(q,dq);

    for (Eigen::Index i_obs=0;i_obs<obstacles_positions_.cols();i_obs++)
    {
      for (size_t i_poi=0;i_poi<poi_poses_in_base.size();i_poi++)
      {
        distance_vector=obstacles_positions_.col(i_obs)-poi_poses_in_base.at(i_poi).translation();
        distance = distance_vector.norm();

        if(distance<self_distance_)
          continue;

        tangential_speed = ((poi_twist_in_base.at(i_poi).block(0,0,3,1)).dot(distance_vector))/distance;

        if(tangential_speed<=0)  // robot is going away
        {
          tmp_scaling_factor = 1.0;
        }
        else if(distance>min_distance_)
        {
          v_safety=std::sqrt(term1_+2.0*max_cart_acc_*distance)-a_t_r_;  //NB: human velocity not considered for now
          tmp_scaling_factor=v_safety/tangential_speed;                  // no division by 0
        }
        else  // distance<=min_distance -> you have found the minimum scaling factor, return
        {
          scaling_factor = 0.0;
          return scaling_factor;
        }

        if(tmp_scaling_factor<scaling_factor)
          scaling_factor = tmp_scaling_factor;
      }
    }

    q = q+delta_q;
  }

  return scaling_factor;
}

SSM15066EstimatorPtr SSM15066Estimator::clone()
{
  SSM15066EstimatorPtr clone = std::make_shared<SSM15066Estimator>(rosdyn::createChain(chain_),max_step_size_,obstacles_positions_);

  clone->setMaxCartAcc(max_cart_acc_);
  clone->setMinDistance(min_distance_);
  clone->setMaxStepSize(max_step_size_);
  clone->setSelfDistance(self_distance_);

  clone->updateMembers();

  return clone;
}


}
