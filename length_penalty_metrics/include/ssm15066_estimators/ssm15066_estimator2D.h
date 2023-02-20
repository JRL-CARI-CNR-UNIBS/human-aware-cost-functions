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

#include <ssm15066_estimators/ssm15066_estimator.h>

namespace ssm15066_estimator
{
class SSM15066Estimator2D;
typedef std::shared_ptr<SSM15066Estimator2D> SSM15066Estimator2DPtr;

/**
  * @brief The SSM15066Estimator2D class is a 2D dSSM estimator, that means that the robot velocity vector towards the human is considered.
  * It computes the scaling factor for each configuration xi along a connection (xs,xg) and then the mean value.
  * The computation is sequential.
  */
class SSM15066Estimator2D: public SSM15066Estimator
{
protected:

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  SSM15066Estimator2D(const rosdyn::ChainPtr &chain, const double& max_step_size=0.05);
  SSM15066Estimator2D(const rosdyn::ChainPtr &chain, const double& max_step_size,
                    const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions);

  virtual double computeScalingFactor(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2) override;
  virtual SSM15066EstimatorPtr clone() override;
};

}
