#pragma once
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

#include <ros/ros.h>
#include <eigen3/Eigen/Core>

/* The purpose of the class is the definition of a template for safety related velocity scaling factor estimator.
 * The goal is to compute an approximation of the scaling factor lambda that the robot will experiment moving from a
 * configuration q1 to a configuration q2, given the obstacles positions (e.g., human's head, arms and torso positions)
*/

 namespace ssm15066
 {
 class VelocityScalingEstimator;
 typedef std::shared_ptr<VelocityScalingEstimator> VelocityScalingEstimatorPtr;


 class VelocityScalingEstimator
 {
 private:
   Eigen::Matrix<double,3,Eigen::Dynamic> obstacles_positions_; //x,y,z (raws) of obstacles (cols). Number of cols depends on the number of obstacles of the scene

 public:
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW
   VelocityScalingEstimator();
   VelocityScalingEstimator(const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions);

   void setObstaclesPositions(const Eigen::Matrix<double,3,Eigen::Dynamic>& obstacles_positions);
   double computeScalingFactor(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2) = 0;
 };

 }
