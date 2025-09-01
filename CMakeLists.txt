cmake_minimum_required(VERSION 3.22)
project(ecg2o LANGUAGES CXX)



set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

#----------------------------------------------------------------------------------------
# Find required packages
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(fmt REQUIRED)
find_package(g2o 1.0.0 REQUIRED)
 
# OpenMP support
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    add_definitions(-DG2O_USE_OPENMP)
endif()

 

 

# Specify the path to g2o include directory and libraries
set(G2O_INCLUDE_DIRS "/usr/include/include" CACHE PATH "g2o include directory")
set(G2O_LIBRARY_DIRS "/usr/include/lib" CACHE PATH "g2o library directory")
set(G2O_LIBRARIES 
    g2o::core 
    g2o::stuff
    g2o::solver_cholmod
)




# Manually specify g2o include directories and libraries
include_directories(${G2O_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include)  
link_directories(${G2O_LIBRARY_DIRS})


# add executable example_gn
add_executable(example_gn examples/example_gn.cpp)
target_link_libraries(example_gn ${G2O_LIBRARIES} Eigen3::Eigen  fmt::fmt) 


# add executable example_al
add_executable(example_al examples/example_al.cpp)
target_link_libraries(example_al ${G2O_LIBRARIES} Eigen3::Eigen  fmt::fmt) 



# add executable example_al
add_executable(example_soft examples/example_soft.cpp)
target_link_libraries(example_soft ${G2O_LIBRARIES} Eigen3::Eigen  fmt::fmt) 



# add executable example_gtsam_gn
add_executable(example_gtsam_gn examples/comparison_gtsam/example_gtsam_gn.cpp)
target_link_libraries(example_gtsam_gn ${G2O_LIBRARIES} Eigen3::Eigen  fmt::fmt) 



# add executable example_gtsam_al
add_executable(example_gtsam_al examples/comparison_gtsam/example_gtsam_al.cpp)
target_link_libraries(example_gtsam_al ${G2O_LIBRARIES} Eigen3::Eigen  fmt::fmt) 



# add executable oc_test
add_executable(oc_test oc_tracking/oc_test.cpp)
target_link_libraries(oc_test ${G2O_LIBRARIES} Eigen3::Eigen  fmt::fmt)
target_include_directories(oc_test PRIVATE ${CMAKE_SOURCE_DIR}/oc_tracking/include)

