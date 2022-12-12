############################################################################
# CMakeLists.txt file for building analysis code for ICAL
# @author
# @mail
############################################################################

cmake_minimum_required(VERSION 3.9)

# set the project name
project(icalanalysis VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_BUILD_TYPE Debug)

# Set the output folder where your program will be created
#set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})
#set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# sub directorires
#add_subdirectory(src)

#including root
find_package(ROOT CONFIG REQUIRED)
find_package(ROOT COMPONENTS Imt Minuit)

#----------------------------------------------------------------------------
# Setup Root include directories and compile definitions
#
include("${ROOT_USE_FILE}")
include("${ROOT_DIR}/modules/RootMacros.cmake")
#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include)
link_directories(${PROJECT_SOURCE_DIR}/src)

add_compile_options(-Wall -Wextra -pedantic)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

separate_arguments(ROOT_EXE_LINKER_FLAGS)
add_executable(icalanalysis multiplicity.cc ${sources} ${headers})

target_link_libraries(icalanalysis PUBLIC ${ROOT_LIBRARIES} ${ROOT_EXE_LINKER_FLAGS})

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# Testing
enable_testing()
