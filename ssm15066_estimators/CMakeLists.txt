cmake_minimum_required(VERSION 3.0.2)
project(ssm15066_estimators)
add_compile_options(-pthread -std=c++17 -funroll-loops -Wall -Ofast)
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

find_package(catkin REQUIRED COMPONENTS
roscpp
graph_core
rosdyn_core
thread-pool
length_penalty_metrics
)
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS roscpp graph_core rosdyn_core thread-pool length_penalty_metrics
  DEPENDS
  )
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  )
add_library(${PROJECT_NAME}
src/ssm15066_estimators/ssm15066_estimator.cpp
src/ssm15066_estimators/ssm15066_estimator1D.cpp
src/ssm15066_estimators/ssm15066_estimator2D.cpp
src/ssm15066_estimators/parallel_ssm15066_estimator2D.cpp
src/min_distance_solvers/min_distance_solver.cpp
)
add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})
