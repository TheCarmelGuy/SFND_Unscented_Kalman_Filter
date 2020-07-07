cmake_minimum_required(VERSION 2.8 FATAL_ERROR)


set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

project(playback)


set(BUILD_MODERN OFF) #Toggle OFF to run w/ original Deps 
set(BUILD_TEST OFF) #Toggle on to Build Unit tests 




if(BUILD_MODERN) 
  find_package(PCL 1.8 REQUIRED)
  add_definitions(-std=c++14)
else()
  find_package(PCL 1.2 REQUIRED) 
  add_definitions(-std=c++11)
endif()



if(BUILD_TEST)
  find_package(GTest REQUIRED)
endif()




include_directories(${GTEST_INCLUDE_DIRS})

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
list(REMOVE_ITEM PCL_LIBRARIES "vtkproj4")


add_executable (ukf_highway src/main.cpp src/ukf.cpp src/tools.cpp src/render/render.cpp)
target_link_libraries (ukf_highway ${PCL_LIBRARIES})

if(BUILD_TEST)
  add_executable(test_ukf src/test_ukf.cpp src/ukf.cpp) 
  target_link_libraries (test_ukf ${GTEST_LIBRARIES} pthread)
endif()
