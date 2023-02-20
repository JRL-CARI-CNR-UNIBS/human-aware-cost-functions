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

#include <ssm15066_estimators/parallel_ssm15066_estimator2D.h>

namespace ssm15066_estimator
{

ParallelSSM15066Estimator2D::ParallelSSM15066Estimator2D(const rosdyn::ChainPtr &chain, const unsigned int& n_threads):
  SSM15066Estimator2D(chain),n_threads_(n_threads){init();}


ParallelSSM15066Estimator2D::ParallelSSM15066Estimator2D(const rosdyn::ChainPtr &chain, const double& max_step_size,
                                                         const unsigned int& n_threads):
  SSM15066Estimator2D(chain,max_step_size),n_threads_(n_threads){init();}

ParallelSSM15066Estimator2D::ParallelSSM15066Estimator2D(const rosdyn::ChainPtr &chain, const double& max_step_size,
                                                         const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions,
                                                         const unsigned int& n_threads):
  SSM15066Estimator2D(chain,max_step_size,obstacles_positions),n_threads_(n_threads){init();}


void ParallelSSM15066Estimator2D::init()
{
  verbose_ = 0;
  stop_ = true;
  running_threads_ = 0;

  if(n_threads_<=0)
    n_threads_ = std::thread::hardware_concurrency();
  else if(n_threads_>std::thread::hardware_concurrency())
  {
    ROS_ERROR_STREAM("number of threads ("<<n_threads_<<") should not be higher than hardware max concurrency ("<<std::thread::hardware_concurrency()<<")");
    n_threads_ = std::thread::hardware_concurrency();
  }

  queues_ .clear();
  chains_ .clear();
  futures_.clear();

  chains_ .resize(n_threads_);
  queues_ .resize(n_threads_);
  futures_.resize(n_threads_);

  for(unsigned int i=0;i<n_threads_;i++)
  {
    chains_[i] = chain_->clone();
    queues_[i] = Queue();
  }

  pool_ = std::make_shared<BS::thread_pool>(n_threads_);
}

void ParallelSSM15066Estimator2D::resetQueues()
{
  stop_ = true;

  pool_->wait_for_tasks();
  std::for_each(queues_.begin(),queues_.end(),[&](Queue& queue){
    queue.reset();
  });

  running_threads_ = 0;

  assert([&]() ->bool{
           for(const Queue& queue:queues_)
           {
             if(queue.queue_.size() != 0)
             {
               return false;
             }
           }
           return true;
         }());

  assert(pool_->get_tasks_total  () == 0);
  assert(pool_->get_tasks_queued () == 0);
  assert(pool_->get_tasks_running() == 0);
}

unsigned int ParallelSSM15066Estimator2D::fillQueues(const Eigen::VectorXd& q1, const Eigen::VectorXd q2)
{
  bool all_threads = false;
  unsigned int n_addends = 0;
  unsigned int thread_iter = 0;
  unsigned int iter = std::ceil((q2-q1).norm()/max_step_size_);

  Eigen::VectorXd q;
  Eigen::VectorXd delta_q = (q2-q1)/iter;

  for(unsigned int i=0;i<iter+1;i++)
  {
    q = q1+i*delta_q;

    if(thread_iter>=n_threads_)
    {
      thread_iter = 0;
      all_threads = true;
    }

    queues_[thread_iter].insert(q);

    n_addends++;
    thread_iter++;
  }

  assert((q2-q).norm()<1e-08);
  assert(n_addends == iter+1);

  if(all_threads)
    running_threads_ = n_threads_;
  else
    running_threads_ = thread_iter;

  return n_addends;
}

double ParallelSSM15066Estimator2D::computeScalingFactor(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2)
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

  tic = ros::WallTime::now();
  unsigned int n_addends = fillQueues(q1,q2);
  toc = ros::WallTime::now();
  time_fill = (toc-tic).toSec();

  if(n_addends == 0)
    return 1.0;

  /* Compute the time of each joint to move from q1 to q2 at its maximum speed and consider the longest time
   * The "slowest" joint will move at its highest speed while the other ones will
   * move at (t_i/slowest_joint_time)*max_speed_i, where slowest_joint_time >= t_i */
  double slowest_joint_time = (inv_max_speed_.cwiseProduct(q2 - q1)).cwiseAbs().maxCoeff();
  dq_max_ = (q2-q1)/slowest_joint_time;

  stop_ = false; //must be after resetQueues()

  tic = ros::WallTime::now();
  for(unsigned int i=0;i<running_threads_;i++)
    futures_[i]= pool_->submit(&ParallelSSM15066Estimator2D::computeScalingFactorAsync,this,i);
  toc = ros::WallTime::now();
  time_thread = (toc-tic).toSec();

  double queue_sum;
  double sum_scaling_factors = 0.0;
  tic = ros::WallTime::now();
  for(unsigned int i=0;i<running_threads_;i++)
  {
    queue_sum = futures_[i].get();

    if(queue_sum == 0.0)
      return 0.0;
    else
      sum_scaling_factors += queue_sum;
  }
  toc = ros::WallTime::now();
  time_join = (toc-tic).toSec();
  time_tot = (toc-tic_init).toSec();

  if(verbose_>0)
    ROS_INFO_STREAM("time reset queues "<<(time_reset/time_tot)*100<<"% || time fill queues "<<(time_fill/time_tot)*100<<"% || time threads creation "<<(time_thread/time_tot)*100<<"% || time executions "<<(time_join/time_tot)*100<<"%");

  return (sum_scaling_factors/((double) n_addends));
}

double ParallelSSM15066Estimator2D::computeScalingFactorAsync(const unsigned int& idx_queue)
{
  Eigen::Vector3d distance_vector, tmp_distance_vector;
  std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> poi_poses_in_base, tmp_pose;
  std::vector<Eigen::Vector6d, Eigen::aligned_allocator<Eigen::Vector6d>> poi_twist_in_base, tmp_twist;
  double distance, tangential_speed, v_safety, scaling_factor, min_scaling_factor_of_q, sum_scaling_factor, tmp_speed;

  rosdyn::ChainPtr chain = chains_[idx_queue];

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
        //consider only links inside the poi_names_ list
        if(std::find(poi_names_.begin(),poi_names_.end(),links_names_[i_poi])>=poi_names_.end())
          continue;

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

          if(scaling_factor<1e-02)
          {
            stop_ = true;
            return 0.0;
          }
        }
        else  // distance<=min_distance -> you have found the minimum scaling factor, return
        {
          stop_ = true;
          return 0.0;
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
          }
        }

        if(stop_)
          break;
      } // end robot poi for loop
      if(stop_)
        break;
    } // end obstacles for loop

    if(verbose_>0)
    {
      mtx_.lock();
      ROS_INFO_STREAM("q "<<q.transpose()<<" || dist vector: "<<tmp_distance_vector.transpose()<<" || tangential speed: "<<tmp_speed<<
                      " || scaling factor at q: "<<min_scaling_factor_of_q);

      if(verbose_ == 2)
      {
        for(Eigen::Affine3d p:poi_poses_in_base)
          ROS_INFO_STREAM("POI poses \n"<<p.matrix());

        for(Eigen::Vector6d t:tmp_twist)
          ROS_INFO_STREAM("POI twists "<<t.transpose());
      }
      mtx_.unlock();
    }

    sum_scaling_factor += min_scaling_factor_of_q;

    if(stop_)
      break;
  } //end q in queue for loop

  return sum_scaling_factor;
}

SSM15066EstimatorPtr ParallelSSM15066Estimator2D::clone()
{
  ParallelSSM15066Estimator2DPtr clone = std::make_shared<ParallelSSM15066Estimator2D>(chain_->clone(),max_step_size_,obstacles_positions_,n_threads_);

  clone->setMaxCartAcc(max_cart_acc_);
  clone->setMinDistance(min_distance_);
  clone->setMaxStepSize(max_step_size_);

  clone->updateMembers();

  return clone;
}

}
