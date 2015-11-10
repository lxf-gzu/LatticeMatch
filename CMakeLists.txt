project(LatticeMatch)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

SET(CMAKE_CXX_FLAGS "-Wall ${CMAKE_CXX_FLAGS} -std=c++11 -DNDEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG "-Wall ${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -O0")

IF(UNIX)
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} m)
ENDIF(UNIX)
