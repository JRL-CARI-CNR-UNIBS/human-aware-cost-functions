#pragma once
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

#include <mutex>
#include <thread>
#include <Eigen/StdVector>
#include <ssm15066_estimators/ssm15066_estimator.h>

namespace ssm15066_estimator
{
class ParallelSSM15066Estimator;
typedef std::shared_ptr<ParallelSSM15066Estimator> ParallelSSM15066EstimatorPtr;

/**
  * @brief The ParallelSSM15066Estimator class is a template for safety related velocity scaling factor (SSM ISO-15066) estimator.
  * The goal is to compute an approximation of the scaling factor that the robot will experiment moving from a
  * configuration q1 to a configuration q2, given the obstacles positions (e.g., human's head, arms and torso positions)
  */
class ParallelSSM15066Estimator: public SSM15066Estimator
{
protected:
#define DEFAULT_THREADS_NUMBER 5

  /**
   * @brief The Queue struct defines an makes more readable the queue struct used by each thread
   */
  struct Queue
  {
    std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > queue_; //don't remove spaces
    void reset(){queue_.clear();}
    void insert(const Eigen::VectorXd& q){queue_.push_back(q);}
  };

  /**
   * @brief These are class members related to threads management
   */
  bool stop_;
  std::mutex mtx_;
  unsigned int n_threads_;
  std::vector<Queue> queues_;
  unsigned int running_threads_;
  std::vector<std::thread> threads_;
  std::vector<rosdyn::ChainPtr> chains_;

  /**
   * @brief dq_max_ is the maximum joint velocity to move from q1 to q2.
   */
  Eigen::VectorXd dq_max_;

  /**
   * @brief scaling_factor_ is the reference scaling factor used to compare the scaling factors computed by the threads
   */
  double ref_scaling_factor_;

  /**
   * @brief init intiializes threds-related class members
   */
  void init();

  /**
   * @brief resetQueues reset queues_
   */
  void resetQueues();

  /**
   * @brief fillQueues fill queue_ with the robot configurations q belonging to (q1,q2) to evaluate
   * @param q1
   * @param q2
   */
  void fillQueues(const Eigen::VectorXd& q1, const Eigen::VectorXd q2);

  /**
   * @brief computeScalingFactorThread computes an approximation of the scaling factor the robot will experience at configuration
   * q, according to SSM ISO-15066. The maximum robot joints' velocities are considered for this computation. It takes q from queue_
   * @param idx_thread.
   * @return the scaling factor.
   */
  void computeScalingFactorThread(const unsigned int& idx_thread);

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  ParallelSSM15066Estimator(const rosdyn::ChainPtr &chain, const unsigned int& n_threads=DEFAULT_THREADS_NUMBER);
  ParallelSSM15066Estimator(const rosdyn::ChainPtr &chain, const double& max_step_size, const unsigned int& n_threads=DEFAULT_THREADS_NUMBER);
  ParallelSSM15066Estimator(const rosdyn::ChainPtr &chain, const double& max_step_size,
                            const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions,
                            const unsigned int& n_threads=DEFAULT_THREADS_NUMBER);

  /**
   * @brief computeWorstCaseScalingFactor computes an approximation of the scaling factor the robot will experience travelling from
   * q1 to q2, according to SSM ISO-15066. The maximum robot joints' velocities are considered for this computation.
   * Ita speeds up the computation using multiple threads
   * @param q1.
   * @param q2.
   * @return the scaling factor.
   */
  double computeScalingFactor(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2);

  inline virtual SSM15066EstimatorPtr clone();

};

}
