cmake_minimum_required(VERSION 3.1)
project(dog_editor)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)

set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -ffast-math -O2")

set(LIBIGL_ROOT "${CMAKE_CURRENT_LIST_DIR}/ext/libigl/")
set(LIBIGL_SOURCE_DIR "${LIBIGL_ROOT}/include")
set(LIBIGL_EXTERNAL "${LIBIGL_ROOT}/external")

list(APPEND CMAKE_MODULE_PATH "${LIBIGL_ROOT}/cmake")
include(libigl)

include_directories(
  ${CMAKE_CURRENT_LIST_DIR}/ext/libigl/include
  ${CMAKE_CURRENT_LIST_DIR}/ext/libigl/external/eigen
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -ffast-math -O2")

# Check to see if PARDISO is available
find_library(PARDISO
  pardiso500-INTEL1301-X86-64
  pardiso500-GNU461-X86-64
  pardiso500-GNU472-X86-64
  pardiso500-GNU481-X86-64
  pardiso600-GNU720-X86-64
  pardiso600-GNU800-X86-64
  pardiso500-MPI-INTEL1301-X86-64
  pardiso500-MPI-GNU450-X86-64
  pardiso500-MPI-GNU461-X86-64
  pardiso500-MPI-GNU463-X86-64
  pardiso500-MPI-GNU472-X86-64
  pardiso500-WIN-X86-64
  pardiso600-WIN-X86-64
  pardiso500-MACOS-X86-64
  pardiso600-MACOS-X86-64
  PATHS ${CMAKE_CURRENT_SOURCE_DIR}/ext/pardiso
)

if (${PARDISO} MATCHES "NOTFOUND")
  message("Could not find Pardiso solver.\n")
else()
  message("Found PARDISO library at ${PARDISO}")
  add_library(libpardiso SHARED IMPORTED)
  set_target_properties(libpardiso PROPERTIES IMPORTED_LOCATION ${PARDISO})
  set(ADDITIONAL_LIBRARIES libpardiso gomp gfortran pthread openblas)
  set(ADDITIONAL_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/Optimization/Solvers/Pardiso/PardisoSolver.cpp)
endif()

# Add project files
file(GLOB SRCFILES src/*.cpp)
file(GLOB SRC_OPTIMIZATION 
      src/Optimization/*.cpp
      src/Optimization/Solvers/*.cpp
      src/Optimization/Solvers/lbfgs/*.cpp)
file(GLOB SRC_QUAD src/QuadMesh/*.cpp)
file(GLOB SRC_DOG 
    src/Dog/*.cpp
    src/Dog/Objectives/*.cpp
    src/Dog/Solvers/*.cpp)
file(GLOB SRC_GUI src/Gui/*.cpp)
file(GLOB SRC_GEOMETRY_PRIMITIVES src/GeometryPrimitives/*.cpp)

add_executable(${PROJECT_NAME}_bin ${SRCFILES} ${SRC_OPTIMIZATION} ${SRC_QUAD} ${SRC_DOG} ${SRC_GUI} ${SRC_GEOMETRY_PRIMITIVES} ${ADDITIONAL_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin PUBLIC igl::core igl::opengl_glfw igl::opengl_glfw_imgui gomp gfortran pthread ${ADDITIONAL_LIBRARIES})