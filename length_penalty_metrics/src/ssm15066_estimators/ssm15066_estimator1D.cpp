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

#include <ssm15066_estimators/ssm15066_estimator1D.h>

namespace ssm15066_estimator
{
SSM15066Estimator1D::SSM15066Estimator1D(const rosdyn::ChainPtr &chain, const double& max_step_size):
  SSM15066Estimator(chain,max_step_size)
{
  min_distance_solver_ = std::make_shared<MinDistanceSolver>(chain);
}

SSM15066Estimator1D::SSM15066Estimator1D(const rosdyn::ChainPtr &chain, const double &max_step_size, const Eigen::Matrix<double,3,Eigen::Dynamic> &obstacles_positions):
  SSM15066Estimator(chain,max_step_size,obstacles_positions)
{
  min_distance_solver_ = std::make_shared<MinDistanceSolver>(chain,obstacles_positions);
}

void SSM15066Estimator1D::addObstaclePosition(const Eigen::Vector3d& obstacle_position)
{
  SSM15066Estimator::addObstaclePosition(obstacle_position);
  min_distance_solver_->addObstaclePosition(obstacle_position);
}

void SSM15066Estimator1D::setObstaclesPositions(const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions)
{
  SSM15066Estimator::setObstaclesPositions(obstacles_positions);
  min_distance_solver_->setObstaclesPositions(obstacles_positions);
}

double SSM15066Estimator1D::computeScalingFactor(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2)
{
  assert(obstacles_positions_ == min_distance_solver_->getObstaclesPositions());
  if(obstacles_positions_.cols()==0)  //no obstacles in the scene
    return 1.0;

  double sum_scaling_factor = 0.0;

  double min_distance, velocity, scaling_factor, min_scaling_factor_of_q, v_safety;
  std::vector< Eigen::Vector6d, Eigen::aligned_allocator<Eigen::Vector6d>> poi_twist_in_base;

  /* Compute the time of each joint to move from q1 to q2 at its maximum speed and consider the longest time */
  double slowest_joint_time = (inv_max_speed_.cwiseProduct(q2 - q1)).cwiseAbs().maxCoeff();

  /* The "slowest" joint will move at its highest speed while the other ones will
   * move at (t_i/slowest_joint_time)*max_speed_i, where slowest_joint_time >= t_i */
  Eigen::VectorXd dq = (q2-q1)/slowest_joint_time;

  unsigned int iter = std::ceil((q2-q1).norm()/max_step_size_);

  Eigen::VectorXd q;
  Eigen::VectorXd delta_q = (q2-q1)/iter;

  for(unsigned int i=0;i<iter+1;i++)
  {
    q = q1+i*delta_q;
    min_scaling_factor_of_q = 1.0;

    poi_twist_in_base = chain_->getTwist(q,dq);
    min_distance = min_distance_solver_->computeMinDistance(q)->distance_;

    v_safety = safeVelocity(min_distance);

    for(size_t i_poi=0;i_poi<poi_twist_in_base.size();i_poi++)
    {
      //consider only links inside the poi_names_ list
      if(std::find(poi_names_.begin(),poi_names_.end(),links_names_[i_poi])>=poi_names_.end())
        continue;

      velocity = poi_twist_in_base[i_poi].norm();

      if(velocity<1e-02)
      {
        scaling_factor = 1.0;
      }
      else if(min_distance>min_distance_)
      {
        scaling_factor = v_safety/velocity; // no division by 0

        if(scaling_factor < 1e-02)
          return 0.0;
      }
      else  // distance<=min_distance -> you have found the minimum scaling factor, return
      {
        return 0.0;  //if one point q has 0.0 scaling factor, return it
      }

      if(scaling_factor<min_scaling_factor_of_q)
      {
        min_scaling_factor_of_q = scaling_factor;
      }
    } // end robot poi for loop

    sum_scaling_factor += min_scaling_factor_of_q;
  }
  assert((q2-q).norm()<1e-08);

  // return the average scaling factor (if no q have zero scaling factor)
  return sum_scaling_factor/((double) iter+1);
}

SSM15066EstimatorPtr SSM15066Estimator1D::clone()
{
  SSM15066Estimator1DPtr clone = std::make_shared<SSM15066Estimator1D>(chain_->clone(),max_step_size_,obstacles_positions_);

  clone->setMaxCartAcc(max_cart_acc_);
  clone->setMinDistance(min_distance_);
  clone->setMaxStepSize(max_step_size_);

  clone->updateMembers();

  return clone;
}


}
