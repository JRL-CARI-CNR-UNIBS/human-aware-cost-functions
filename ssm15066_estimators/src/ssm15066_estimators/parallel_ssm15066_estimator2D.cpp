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
    queues_[i] = std::make_shared<Queue>();
  }

  pool_ = std::make_shared<BS::thread_pool>(n_threads_);
}

void ParallelSSM15066Estimator2D::resetQueues()
{
  stop_ = true;

  pool_->wait_for_tasks();
  assert(pool_->get_tasks_total  () == 0);
  assert(pool_->get_tasks_queued () == 0);
  assert(pool_->get_tasks_running() == 0);

  for(const QueuePtr& queue:queues_)
    queue->reset();

  running_threads_ = 0;

  assert([&]() ->bool{
           for(const QueuePtr& queue:queues_)
           {
             if(not queue->queue_.empty())
             return false;
           }
           return true;
         }());
}

unsigned int ParallelSSM15066Estimator2D::fillQueues(const Eigen::VectorXd& q1, const Eigen::VectorXd q2)
{
  unsigned int n_addends = 0;
  unsigned int thread_iter = 0;

  Eigen::VectorXd connection_vector = q2-q1;
  unsigned int iter = std::max(std::ceil(connection_vector.norm()/max_step_size_),1.0);

  Eigen::VectorXd q;
  Eigen::VectorXd delta_q = connection_vector/iter;

  assert([&]() ->bool{
           for(const QueuePtr& queue: queues_)
           {
             if(not queue->queue_.empty())
             return false;
           }

           return true;
         }());

  for(unsigned int i=0;i<iter+1;i++)
  {
    q = q1+i*delta_q;

    if(thread_iter>=n_threads_)
      thread_iter = 0;

    queues_[thread_iter]->insert(q);

    n_addends++;
    thread_iter++;
  }

  assert([&]() ->bool{
           double diff = (q2-q).norm();
           if(diff<1e-08)
           {
             return true;
           }
           else
           {
             ROS_INFO_STREAM("Diff "<<diff);
             ROS_INFO_STREAM("q "<<q.transpose()<<"\nq2 "<<q2.transpose());
             ROS_INFO_STREAM("iter "<<iter);
             ROS_INFO_STREAM("iter+1 "<<(iter+1)<<" delta_q "<<delta_q.transpose()<<" product "<<((iter+1)*delta_q).transpose());
             return false;
           }
         }());

  if(n_addends>n_threads_)
    running_threads_ = n_threads_;
  else
    running_threads_ = n_addends;

  if(verbose_>0)
  {
    for(unsigned int i=0;i<queues_.size();i++)
    {
      ROS_WARN_STREAM("QUEUE "<<i);
      for(const Eigen::VectorXd& q:queues_[i]->queue_)
        ROS_INFO_STREAM(q.transpose());
    }
  }

  return n_addends;
}

double ParallelSSM15066Estimator2D::computeScalingFactor(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2)
{
  if(verbose_>0)
    ROS_WARN("--------");

  ros::WallTime tic, tic_init, toc;
  double time_tot, time_reset, time_fill, time_thread, time_join;

  if(obstacles_positions_.cols()==0)  //no obstacles in the scene
  {
    if(verbose_>0)
      ROS_ERROR("--------");

    return 1.0;
  }

  if(verbose_>0)
  {
    ROS_ERROR_STREAM("number of obstacles: "<<obstacles_positions_.cols()<<", number of poi: "<<poi_names_.size());
    for(unsigned int i=0;i<obstacles_positions_.cols();i++)
      ROS_ERROR_STREAM("obs location -> "<<obstacles_positions_.col(i).transpose());
  }

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
  {
    if(verbose_>0)
      ROS_ERROR("--------");

    return 1.0;
  }

  /* Compute the time of each joint to move from q1 to q2 at its maximum speed and consider the longest time
   * The "slowest" joint will move at its highest speed while the other ones will
   * move at (t_i/slowest_joint_time)*max_speed_i, where slowest_joint_time >= t_i */
  Eigen::VectorXd connection_vector = q2-q1;
  double slowest_joint_time = (inv_max_speed_.cwiseProduct(q2 - q1)).cwiseAbs().maxCoeff();
  dq_max_ = connection_vector/slowest_joint_time;

  assert([&]() ->bool{
           Eigen::VectorXd q_v  = connection_vector/connection_vector.norm();
           Eigen::VectorXd dq_v = dq_max_/dq_max_.norm();

           double err = (q_v-dq_v).norm();
           if(err<1e-08)
           {
             return true;
           }
           else
           {
             ROS_ERROR_STREAM("q_v "<<q_v.transpose()<<" dq_v "<<dq_v.transpose()<<" err "<<err<<" slowest time "<<slowest_joint_time);
             ROS_ERROR_STREAM("q1 "<<q1.transpose()<<" q2 "<<q2.transpose()<<" dq_inv "<<inv_max_speed_.transpose());

             return false;
           }
         }());

  if(verbose_>0)
    ROS_ERROR_STREAM("joint velocity "<<dq_max_.norm());

  stop_ = false; //must be after resetQueues()

  tic = ros::WallTime::now();
  for(unsigned int i=0;i<running_threads_;i++)
    futures_[i]= pool_->submit(&ParallelSSM15066Estimator2D::computeScalingFactorAsync,this,i);
  toc = ros::WallTime::now();
  time_thread = (toc-tic).toSec();

  assert(pool_->get_tasks_total() <= running_threads_);

  double queue_sum;
  double sum_scaling_factors = 0.0;
  tic = ros::WallTime::now();

  pool_->wait_for_tasks();
  for(unsigned int i=0;i<running_threads_;i++)
  {
    queue_sum = futures_[i].get();
    sum_scaling_factors += queue_sum;

    if(verbose_>0)
      ROS_INFO_STREAM("thread "<<i<<" finished");
  }
  stop_ = true;

  toc = ros::WallTime::now();
  time_join = (toc-tic).toSec();
  time_tot = (toc-tic_init).toSec();

  assert([&]() ->bool{
           if(pool_->get_tasks_running() != 0)
           {
             ROS_INFO_STREAM("running threads "<<pool_->get_tasks_running());
             return false;
           }
           return true;
         }());

  double scaling_factor = (sum_scaling_factors/((double) n_addends));

  assert([&]() ->bool{
           if(scaling_factor>=1.0)
           {
             return true;
           }
           else
           {
             ROS_INFO_STREAM("Scaling factor "<<scaling_factor);
             ROS_INFO_STREAM("sum "<<sum_scaling_factors<<" addends "<<(double) n_addends<<" iter "<< std::max(std::ceil(connection_vector.norm()/max_step_size_),1.0));
             return false;
           }
         }());

  if(verbose_>0)
  {
    if(verbose_>1)
      ROS_INFO_STREAM("time reset queues "<<(time_reset/time_tot)*100<<"% || time fill queues "<<(time_fill/time_tot)*100<<"% || time threads creation "<<(time_thread/time_tot)*100<<"% || time executions "<<(time_join/time_tot)*100<<"%");

    ROS_ERROR("--------");
  }

  return scaling_factor;
}

double ParallelSSM15066Estimator2D::computeScalingFactorAsync(const unsigned int& idx_queue)
{
  Eigen::Vector3d distance_vector;
  std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> poi_poses_in_base;
  std::vector<Eigen::Vector6d, Eigen::aligned_allocator<Eigen::Vector6d>> poi_twist_in_base;
  double distance, tangential_speed, v_safety, scaling_factor, max_scaling_factor_of_q, sum_scaling_factor;

  rosdyn::ChainPtr chain = chains_[idx_queue];

  sum_scaling_factor = 0.0;
  for(const Eigen::VectorXd& q: queues_[idx_queue]->queue_)
  {
    max_scaling_factor_of_q = 1.0;

    poi_twist_in_base = chain->getTwist(q,dq_max_);
    poi_poses_in_base = chain->getTransformations(q);

    if(verbose_>0)
      ROS_INFO_STREAM("q -> "<<q.transpose()<<" from queue "<<idx_queue);

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

          if(v_safety == 0.0)
          {
            if(verbose_>0)
              ROS_INFO("stop -> v_safety = 0");

            stop_ = true;
            return std::numeric_limits<double>::infinity();
          }
          else
            scaling_factor = tangential_speed/v_safety; // no division by 0

          assert(v_safety>=0.0);
        }
        else  // distance<=min_distance -> you have found the maximum scaling factor, return
        {
          if(verbose_>0)
            ROS_INFO("stop -> distance < min_distance");

          stop_ = true;
          return std::numeric_limits<double>::infinity();
        }

        if(scaling_factor>max_scaling_factor_of_q)
          max_scaling_factor_of_q = scaling_factor;

        if(stop_)
          break;
      } // end robot poi for-loop
      if(stop_)
        break;
    } // end obstacles for-loop

    if(verbose_>0)
      ROS_INFO_STREAM("q "<<q.transpose()<<" -> scaling factor: "<<max_scaling_factor_of_q);

    sum_scaling_factor += max_scaling_factor_of_q;

    if(stop_)
      break;
  } //end q in queue for-loop

  return sum_scaling_factor;
}

pathplan::CostPenaltyPtr ParallelSSM15066Estimator2D::clone()
{
  ParallelSSM15066Estimator2DPtr cloned_ssm = std::make_shared<ParallelSSM15066Estimator2D>(chain_->clone(),max_step_size_,obstacles_positions_,n_threads_);

  cloned_ssm->setPoiNames(poi_names_);
  cloned_ssm->setMaxStepSize(max_step_size_);
  cloned_ssm->setObstaclesPositions(obstacles_positions_);

  cloned_ssm->setMaxCartAcc(max_cart_acc_,false);
  cloned_ssm->setMinDistance(min_distance_,false);
  cloned_ssm->setReactionTime(reaction_time_,false);
  cloned_ssm->setHumanVelocity(human_velocity_,false);

  cloned_ssm->updateMembers();

  pathplan::CostPenaltyPtr clone = cloned_ssm;

  return clone;
}

}
