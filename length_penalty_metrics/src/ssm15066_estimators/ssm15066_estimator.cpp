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

SSM15066Estimator::SSM15066Estimator(const MinDistanceSolverPtr& distance_calculator):
  distance_calculator_(distance_calculator)
{
  chain_ = distance_calculator_->getChain();
  inv_max_speed_ = chain_->getDQMax().cwiseInverse();

  self_distance_ = 0.00;
  min_distance_  = 0.30;
  max_cart_acc_  = 0.10;
  t_r_           = 0.15;

  a_t_r_ = max_cart_acc_*t_r_;
  term1_=std::pow(a_t_r_,2)-2*max_cart_acc_*min_distance_;
}

double SSM15066Estimator::computeScalingFactor(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2)
{
  DistancePtr distance_solver_res;
  distance_solver_res = distance_calculator_->computeMinDistance(q1,q2);

  double scaling_factor = 1.0;

  if(distance_solver_res->distance_ == std::numeric_limits<double>::infinity()) //no obstacles in the scene
    return scaling_factor;
  else
  {
    assert(abs(distance_solver_res->distance_-distance_solver_res->distance_vector_.norm())<1e-06);

    /* Compute the cartesian robot velocity doing these two approximations:
     * 1) Consider the maximum velocity of each joint
     * 2) Consider the Jacobian at the configuration which is the closest one the objects */

    /* Compute the time of each joint to move from q1 to q2
     * at its maximum speed and consider the longest time */
    double slowest_joint_time = (inv_max_speed_.cwiseProduct(q2 - q1)).cwiseAbs().maxCoeff();

    /* The "slowest" joint will move at its highest speed while the other ones will
     * move at (t_i/slowest_joint_time)*max_speed_i, where slowest_joint_time >= t_i */
    Eigen::VectorXd dq = (q2-q1)/slowest_joint_time;

    Eigen::VectorXd q = distance_solver_res->robot_configuration_;
    std::vector< Eigen::Vector6d, Eigen::aligned_allocator<Eigen::Vector6d> > poi_twist_in_base = chain_->getTwist(q,dq);

    double distance = distance_solver_res->distance_;

    if(distance<self_distance_)
      return scaling_factor;

    Eigen::Vector3d poi_velocity_vector = poi_twist_in_base.at(distance_solver_res->robot_poi_).block(0,0,3,1);
    double tangential_speed = (poi_velocity_vector).dot(distance_solver_res->distance_vector_)/distance;

    if(tangential_speed<=0)
      return scaling_factor;
    else if(distance>min_distance_)
    {
      double v_max = std::sqrt(term1_+2.0*max_cart_acc_*distance)-a_t_r_;  //NB: human velocity not considered for now
      scaling_factor = std::min(v_max/tangential_speed,1.0);

      assert(scaling_factor>0);
    }
    else  //distance<=min_distance
    {
      scaling_factor = 0.0;
    }

    return scaling_factor;
  }
}
}
