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

#include <ssm15066_estimators/ssm15066_estimator2D.h>

namespace ssm15066_estimator
{

SSM15066Estimator2D::SSM15066Estimator2D(const rosdyn::ChainPtr &chain, const double& max_step_size):
  SSM15066Estimator(chain,max_step_size){}

SSM15066Estimator2D::SSM15066Estimator2D(const rosdyn::ChainPtr &chain, const double &max_step_size, const Eigen::Matrix<double,3,Eigen::Dynamic> &obstacles_positions):
  SSM15066Estimator(chain,max_step_size,obstacles_positions){}

double SSM15066Estimator2D::computeScalingFactor(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2)
{
  if(obstacles_positions_.cols()==0)  //no obstacles in the scene
    return 1.0;

  double sum_scaling_factor = 0.0;

  bool tmp_saved = false;
  Eigen::Vector3d distance_vector, tmp_distance_vector;
  double distance, tangential_speed, scaling_factor, min_scaling_factor_of_q, v_safety, tmp_speed;
  std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>>  poi_poses_in_base, tmp_pose;
  std::vector< Eigen::Vector6d, Eigen::aligned_allocator<Eigen::Vector6d>> poi_twist_in_base, tmp_twist;

  /* Compute the time of each joint to move from q1 to q2 at its maximum speed and consider the longest time */
  double slowest_joint_time = (inv_max_speed_.cwiseProduct(q2 - q1)).cwiseAbs().maxCoeff();

  /* The "slowest" joint will move at its highest speed while the other ones will
   * move at (t_i/slowest_joint_time)*max_speed_i, where slowest_joint_time >= t_i */
  Eigen::VectorXd dq = (q2-q1)/slowest_joint_time;

  unsigned int iter = std::max(std::ceil((q2-q1).norm()/max_step_size_),1.0);

  Eigen::VectorXd q;
  Eigen::VectorXd delta_q = (q2-q1)/iter;

  if(verbose_>0)
    ROS_ERROR_STREAM("delta q "<<delta_q.transpose()<<" iter "<<iter);

  for(unsigned int i=0;i<iter+1;i++)
  {
    q = q1+i*delta_q;

    min_scaling_factor_of_q = 1.0;

    poi_twist_in_base = chain_->getTwist(q,dq);
    poi_poses_in_base = chain_->getTransformations(q);

    for(Eigen::Index i_obs=0;i_obs<obstacles_positions_.cols();i_obs++)
    {
      for(size_t i_poi=0;i_poi<poi_poses_in_base.size();i_poi++)
      {
        //consider only links inside the poi_names_ list
        if(std::find(poi_names_.begin(),poi_names_.end(),links_names_[i_poi])>=poi_names_.end())
          continue;

        if(verbose_)
        {
          Eigen::Vector3d diff = obstacles_positions_.col(i_obs)-poi_poses_in_base.at(i_poi).translation();
          ROS_ERROR_STREAM("poi name "<<links_names_[i_poi]);
          ROS_ERROR_STREAM("poi pos "<<poi_poses_in_base.at(i_poi).translation().transpose()<<" obs pos "<<obstacles_positions_.col(i_obs).transpose()<<" diff "<<diff.transpose());
        }

        distance_vector = obstacles_positions_.col(i_obs)-poi_poses_in_base.at(i_poi).translation();
        distance = distance_vector.norm();

        tangential_speed = ((poi_twist_in_base.at(i_poi).block(0,0,3,1)).dot(distance_vector))/distance;

        if(tangential_speed<=0)  // robot is going away
        {
          scaling_factor = 1.0;
        }
        else if(distance>min_distance_)
        {
          v_safety = safeVelocity(distance);
          scaling_factor = v_safety/tangential_speed; // no division by 0

          assert(v_safety>=0);

          if(scaling_factor < 1e-02)
          {
            if(verbose_)
              ROS_ERROR("scaling factor close to 0");

            return 0.0;
          }
        }
        else  // distance<=min_distance -> you have found the minimum scaling factor, return
        {
          if(verbose_)
            ROS_ERROR("scaling factor 0");

          return 0.0;  //if one point q has 0.0 scaling factor, return it
        }

        if(scaling_factor<min_scaling_factor_of_q)
        {
          min_scaling_factor_of_q = scaling_factor;

          if(verbose_>0)
          {
            tmp_distance_vector = distance_vector;
            tmp_speed = tangential_speed;
            tmp_pose = poi_poses_in_base;
            tmp_twist = poi_twist_in_base;
            tmp_saved = true;
          }
        }
      } // end robot poi for loop
    } // end obstacles for loop

    if(verbose_>0)
    {
      if(tmp_saved)
        ROS_ERROR_STREAM("q "<<q.transpose()<<" || dist vector: "<<tmp_distance_vector.transpose()<<" || tangential speed: "<<tmp_speed<<
                         " || scaling factor at q: "<<min_scaling_factor_of_q);

      if(verbose_ ==2 )
      {
        for(Eigen::Affine3d p:tmp_pose)
          ROS_ERROR_STREAM("POI poses \n"<<p.matrix());

        for(Eigen::Vector6d t:tmp_twist)
          ROS_ERROR_STREAM("POI twists "<<t.transpose());
      }

      ROS_ERROR("------------");
    }

    sum_scaling_factor += min_scaling_factor_of_q;
  } //end q for loop

  assert([&]() ->bool{
           double err = (q2-q).norm();
           if(err<1e-08)
           {
             return true;
           }
           else
           {
             ROS_INFO_STREAM("error "<<err<<" q "<<q.transpose()<<" q2 "<<q2.transpose());
             return false;
           }
         }());

  // return the average scaling factor (if no q have zero scaling factor)
  double res = sum_scaling_factor/((double) iter+1);
  assert([&]() ->bool{
           if(res>=0 && res<=1)
           {
             return true;
           }
           else
           {
             ROS_INFO_STREAM("Scaling factor "<<res<<" sum "<<sum_scaling_factor<<" denominator "<<iter+1);
             return false;
           }
         }());

  return res;
}

SSM15066EstimatorPtr SSM15066Estimator2D::clone()
{
  SSM15066Estimator2DPtr clone = std::make_shared<SSM15066Estimator2D>(chain_->clone(),max_step_size_,obstacles_positions_);

  clone->setMaxCartAcc(max_cart_acc_);
  clone->setMinDistance(min_distance_);
  clone->setMaxStepSize(max_step_size_);

  clone->updateMembers();

  return clone;
}


}
