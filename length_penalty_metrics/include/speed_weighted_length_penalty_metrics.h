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

#include <length_penalty_metrics.h>

namespace pathplan
{
/**
 * @brief The SpeedWeightedLengthPenaltyMetrics class computes the Euclidean distance between two nodes, weighted on the maximum speed of each joint, and penalizes it based on a penalty.
 *
 *                            c(q1,q2) = ||(q2-q1)./q'_max||*lambda,     with lambda >= 1.0
 *
 *                    The penalty lambda is computed by the CostPenalty class.
*/

class SpeedWeightedLengthPenaltyMetrics: public LengthPenaltyMetrics
{
protected:
  Eigen::VectorXd qp_max_; //max joints velocities
  Eigen::VectorXd inv_qp_max_; //inverse of max joints velocities

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  SpeedWeightedLengthPenaltyMetrics(const CostPenaltyPtr& penalizer, const Eigen::VectorXd& qp_max);

  void setMaxJointsSpeed(const Eigen::VectorXd& qp_max)
  {
    qp_max_ = qp_max;
  }

  virtual double cost(const NodePtr& node1,
                      const NodePtr& node2);

  virtual double cost(const Eigen::VectorXd& configuration1,
                      const Eigen::VectorXd& configuration2);


  virtual double utopia(const NodePtr& node1,
                        const NodePtr& node2);

  virtual double utopia(const Eigen::VectorXd& configuration1,
                        const Eigen::VectorXd& configuration2);

  virtual MetricsPtr clone();
};
}
