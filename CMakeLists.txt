cmake_minimum_required(VERSION 3.24)

project(LLas CXX)

# #################################################
# #### Compile Options ############################
# #################################################
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_POSITION_INDEPENDENT_CODE on)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DLLAS_MEASURE_TIME")
  set(CMAKE_CXX_FLAGS_DEBUG "/Z7 /Wall /DLLAS_PRINT_BYTES /DLLAS_LOG_DEBUG /DLLAS_CHECK_BYTE_SIZE")
  set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DLLAS_LOG_INFO")
  elseif(UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLLAS_MEASURE_TIME")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -DLLAS_PRINT_BYTES -DLLAS_LOG_DEBUG -DLLAS_CHECK_BYTE_SIZE")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DLLAS_LOG_INFO")
endif()


# #################################################
# #### Envirionment Variables #####################
# #################################################
# #### Set environment variables in project
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(PROJECT_TEST_DIR ${PROJECT_SOURCE_DIR}/test)

# #################################################
# #### External Libraries #########################
# #################################################
# Nothing

# #################################################
# #### Test Projects ##############################
# #################################################
set(PROJECT_NAME_TEST_LLAS_READ test_llas_read)

project(${PROJECT_NAME_TEST_LLAS_READ} CXX)

add_executable(
  ${PROJECT_NAME_TEST_LLAS_READ}
  "${PROJECT_SOURCE_DIR}/test/read.cpp"
)

target_include_directories(
  ${PROJECT_NAME_TEST_LLAS_READ}
  PRIVATE
  ${PROJECT_INCLUDE_DIR}
)

# #################################################
# #### Install ####################################
# #################################################
install(
	DIRECTORY ${PROJECT_INCLUDE_DIR}
	DESTINATION ${CMAKE_INSTALL_PREFIX}
	FILES_MATCHING PATTERN "*.hpp"
)

# #################################################
# #### Message ####################################
# #################################################
message("# =======================================================================================================")
message("#  '${PROJECT_NAME}' build configuration ")
message("# =======================================================================================================")
message("#  [CMake]")
message("#    CMake version                        : ${CMAKE_VERSION}")
message("#    CMake build type                     : ${CMAKE_BUILD_TYPE}")
message("#    CMAKE_CXX_FLAGS                      : ${CMAKE_CXX_FLAGS}")
message("#    CMAKE_CXX_FLAGS_DEBUG                : ${CMAKE_CXX_FLAGS_DEBUG}")
message("#    CMAKE_CXX_FLAGS_RELEASE              : ${CMAKE_CXX_FLAGS_RELEASE}")
message("# ")
message("#  [C/C++]")
message("#    C++ Compiler                         : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} ${CMAKE_CXX_COMPILER_WRAPPER} ${CMAKE_CXX_COMPILER}")
message("#    CXX STANDARD                         : ${CMAKE_CXX_STANDARD}")
message("# ")
message("#  [Environment Variables]")
message("#    PROJECT_INCLUDE_DIR                  : ${PROJECT_INCLUDE_DIR}")
message("#    PROJECT_TEST_DIR                     : ${PROJECT_TEST_DIR}")
message("# =======================================================================================================")
