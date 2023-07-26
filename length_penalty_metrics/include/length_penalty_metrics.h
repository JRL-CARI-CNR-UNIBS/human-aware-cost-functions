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

#include <graph_core/metrics.h>
#include <graph_core/graph/node.h>

namespace pathplan
{
class CostPenalty;
typedef std::shared_ptr<CostPenalty> CostPenaltyPtr;

class LengthPenaltyMetrics;
typedef std::shared_ptr<LengthPenaltyMetrics> LengthPenaltyMetricsPtr;

/**
 * @brief The LengthPenaltyMetrics class computes the Euclidean distance between two nodes, and penalizes it based on a penalty.
 * Each joint can be weighted using a scale (default set to 1)
 *
 *                            c(q1,q2) = ||(q2-q1).*scale||*lambda,     with lambda >= 1.0
 *
 *                    The penalty lambda is computed by the CostPenalty class.
*/


class LengthPenaltyMetrics: public Metrics
{
protected:
  CostPenaltyPtr penalizer_; // computes lambda
  Eigen::VectorXd scale_; // scales the distance vector

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  static constexpr double lambda_penalty_ = 1.0e12;

  LengthPenaltyMetrics(const CostPenaltyPtr& penalizer, const Eigen::VectorXd& scale);

  CostPenaltyPtr getPenalizer()
  {
    return penalizer_;
  }

  void setPenalizer(const CostPenaltyPtr& penalizer)
  {
    penalizer_ = penalizer;
  }

  void setScale(const Eigen::VectorXd& scale)
  {
    scale_ = scale;
  }

  Eigen::VectorXd getScale()
  {
    return scale_;
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

/**
 * @brief The CostPenalty class is a template class for the computation of the penalty lambda
 */
class CostPenalty
{
protected:
  unsigned int verbose_;

  /**
   * @brief computePenalty computes the penalty
   * @param q1 parent configuration
   * @param q2 child configuration
   * @return the penalty computed
   */
  virtual double computePenalty(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2) = 0;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  CostPenalty()
  {
    verbose_ = 0;
  }

  virtual void setVerbose(const unsigned int& verbose)
  {
    verbose_ = verbose;
  }

  /**
   * @brief getPenalty computes and return the penalty computed
   * @param q1 parent configuration
   * @param q2 child configuration
   * @return the penalty computed
   */
  virtual double getPenalty(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2)
  {
    double penalty = computePenalty(q1,q2);
    assert(penalty >= 1.0);

    return penalty;
  }

  /**
   * @brief clone clones the object
   * @return a cloned object
   */
  virtual CostPenaltyPtr clone() = 0;
};

}
