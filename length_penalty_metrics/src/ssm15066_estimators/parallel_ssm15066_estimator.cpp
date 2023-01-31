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

#include <ssm15066_estimators/parallel_ssm15066_estimator.h>

namespace ssm15066_estimator
{

ParallelSSM15066Estimator::ParallelSSM15066Estimator(const rosdyn::ChainPtr &chain, const unsigned int& n_threads):
  SSM15066Estimator(chain),n_threads_(n_threads){init();}


ParallelSSM15066Estimator::ParallelSSM15066Estimator(const rosdyn::ChainPtr &chain, const double& max_step_size,
                                                     const unsigned int& n_threads):
  SSM15066Estimator(chain,max_step_size),n_threads_(n_threads){init();}

ParallelSSM15066Estimator::ParallelSSM15066Estimator(const rosdyn::ChainPtr &chain, const double& max_step_size,
                                                     const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions,
                                                     const unsigned int& n_threads):
  SSM15066Estimator(chain,max_step_size,obstacles_positions),n_threads_(n_threads){init();}

void ParallelSSM15066Estimator::init()
{
  stop_ = false;
  running_threads_ = 0;

  if(n_threads_<=0)
    n_threads_ = DEFAULT_THREADS_NUMBER;

  queues_ .clear();
  chains_ .clear();
  threads_.clear();

  queues_ .resize(n_threads_);
  threads_.resize(n_threads_);

  for(unsigned int i=0;i<n_threads_;i++)
  {
    chains_.push_back(rosdyn::createChain(chain_));
    queues_.push_back(Queue());
  }
}

void ParallelSSM15066Estimator::resetQueues()
{
  stop_ = true;

  for(unsigned int i=0;i<running_threads_;i++)
  {
    threads_[i].wait();
    queues_[i].reset();
  }

  running_threads_ = 0;
}

unsigned int ParallelSSM15066Estimator::fillQueues(const Eigen::VectorXd& q1, const Eigen::VectorXd q2)
{
  bool all_threads = false;
  unsigned int n_addends = 0;
  unsigned int thread_iter = 0;
  unsigned int iter = std::ceil((q2-q1).norm()/max_step_size_);

  Eigen::VectorXd q = q1;
  Eigen::VectorXd delta_q = (q2-q1)/iter;

  for(unsigned int i=0;i<iter+1;i++)
  {
    if(thread_iter>=n_threads_)
    {
      thread_iter = 0;
      all_threads = true;
    }

    queues_[thread_iter].insert(q);

    n_addends++;
    thread_iter++;

    q = q+delta_q;
  }

  assert((q2-q).norm()<1e-10);
  assert(n_addends == iter+1);

  if(all_threads)
    running_threads_ = n_threads_;
  else
    running_threads_ = thread_iter;

  return n_addends;
}

double ParallelSSM15066Estimator::computeScalingFactor(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2)
{
  ros::WallTime tic, tic_init, toc;
  double time_tot, time_reset, time_fill, time_thread, time_join;

  if(obstacles_positions_.cols()==0)  //no obstacles in the scene
    return 1.0;

  tic_init = ros::WallTime::now();
  tic = ros::WallTime::now();
  resetQueues();
  toc = ros::WallTime::now();
  time_reset = (toc-tic).toSec();

  stop_ = false;                  //must be after resetQueues()

  /* Compute the time of each joint to move from q1 to q2 at its maximum speed and consider the longest timeà
   * The "slowest" joint will move at its highest speed while the other ones will
   * move at (t_i/slowest_joint_time)*max_speed_i, where slowest_joint_time >= t_i */
  double slowest_joint_time = (inv_max_speed_.cwiseProduct(q2 - q1)).cwiseAbs().maxCoeff();
  dq_max_ = (q2-q1)/slowest_joint_time;

  tic = ros::WallTime::now();
  unsigned int n_addends = fillQueues(q1,q2);
  toc = ros::WallTime::now();
  time_fill = (toc-tic).toSec();

  if(n_addends == 0)
    return 1.0;

  tic = ros::WallTime::now();
  for(unsigned int i=0;i<running_threads_;i++)
    threads_[i] = std::async(std::launch::async,&ParallelSSM15066Estimator::computeScalingFactorAsync,this,i);
  toc = ros::WallTime::now();
  time_thread = (toc-tic).toSec();

  double queue_sum;
  double sum_scaling_factors = 0.0;
  tic = ros::WallTime::now();
  for(unsigned int i=0;i<running_threads_;i++)
  {
    queue_sum = threads_[i].get();

    if(queue_sum == 0.0)
      return 0.0;
    else
      sum_scaling_factors += queue_sum;
  }
  toc = ros::WallTime::now();
  time_join = (toc-tic).toSec();
  time_tot = (toc-tic_init).toSec();

  ROS_INFO_STREAM("time reset queues "<<(time_reset/time_tot)*100<<"% || time fill queues "<<(time_fill/time_tot)*100<<"% || time threads creation "<<(time_thread/time_tot)*100<<"% || time executions "<<(time_join/time_tot)*100<<"%");

  return (sum_scaling_factors/((double) n_addends));
}

double ParallelSSM15066Estimator::computeScalingFactorAsync(const unsigned int& idx_queue)
{
  rosdyn::ChainPtr chain = chains_[idx_queue];

  Eigen::Vector3d distance_vector;
  std::vector<Eigen::Vector6d, Eigen::aligned_allocator<Eigen::Vector6d>> poi_twist_in_base;
  std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> poi_poses_in_base;
  double distance, tangential_speed, v_safety, scaling_factor, min_scaling_factor_of_q, sum_scaling_factor;

  sum_scaling_factor = 0.0;

  for(const Eigen::VectorXd& q: queues_[idx_queue].queue_)
  {
    min_scaling_factor_of_q = 1.0;

    poi_twist_in_base = chain->getTwist(q,dq_max_);
    poi_poses_in_base = chain->getTransformations(q);

    for(Eigen::Index i_obs=0;i_obs<obstacles_positions_.cols();i_obs++)
    {
      for(size_t i_poi=0;i_poi<poi_poses_in_base.size();i_poi++)
      {
        distance_vector = obstacles_positions_.col(i_obs)-poi_poses_in_base.at(i_poi).translation();
        distance = distance_vector.norm();

        if(distance<self_distance_)
          continue;

        tangential_speed = ((poi_twist_in_base.at(i_poi).block(0,0,3,1)).dot(distance_vector))/distance;

        if(tangential_speed<=0)  // robot is going away
        {
          scaling_factor = 1.0;
        }
        else if(distance>min_distance_)
        {
          v_safety=std::sqrt(term1_+2.0*max_cart_acc_*distance)-a_t_r_;  //NB: human velocity not considered for now
          scaling_factor=v_safety/tangential_speed;                      // no division by 0

          if(scaling_factor < 1e-02)
          {
            sum_scaling_factor = 0.0;
            stop_ = true;

            return sum_scaling_factor;
          }
        }
        else  // distance<=min_distance -> you have found the minimum scaling factor, return
        {
          sum_scaling_factor = 0.0;
          stop_ = true;

          return sum_scaling_factor;
        }

        if(scaling_factor<min_scaling_factor_of_q)
          min_scaling_factor_of_q = scaling_factor;

        if(stop_)
          break;
      } // end for robot links
      if(stop_)
        break;
    } // end for obstacles
    if(stop_)
      break;

    sum_scaling_factor += min_scaling_factor_of_q;
  } //end for q in queue

  return sum_scaling_factor;
}

SSM15066EstimatorPtr ParallelSSM15066Estimator::clone()
{
  ParallelSSM15066EstimatorPtr clone = std::make_shared<ParallelSSM15066Estimator>(rosdyn::createChain(chain_),max_step_size_,obstacles_positions_,n_threads_);

  clone->setMaxCartAcc(max_cart_acc_);
  clone->setMinDistance(min_distance_);
  clone->setMaxStepSize(max_step_size_);
  clone->setSelfDistance(self_distance_);

  clone->updateMembers();

  return clone;
}

}
