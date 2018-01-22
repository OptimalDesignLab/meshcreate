#This file shows how to link to a PUMI
#installation using CMake
#it represents a simple 'CMakeLists.txt'
#file for a new project

cmake_minimum_required(VERSION 3.0.0)

project(MeshCreate VERSION 1.0.0 LANGUAGES CXX)

# Starting here are the critical lines:

# Allow the user to indicate where they installed SCOREC
# via "-DSCOREC_PREFIX=/home/somewhere" when calling `cmake`
set(SCOREC_PREFIX "" CACHE STRING "Directory where SCOREC is installed")

# If SCOREC_PREFIX was specified, only link to that directory,
# i.e. don't link to another installation in /usr/lib by mistake
if (SCOREC_PREFIX)
  find_package(SCOREC 2.1.0 REQUIRED CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)
else()
# IF SCOREC_PREFIX was not specified, look in typical system directories,
# and also in CMAKE_PREFIX_PATH (environment variable)
  find_package(
      SCOREC #package name, has to be SCOREC
      2.1.0  #version. can be omitted, and will match any installed version
             #greater than or equal to this one, as long as the major number
             #is the same
      REQUIRED #indicate that SCOREC is really needed to compile
      CONFIG   #skip the 'MODULE' search system, save some time and confusion
      )
endif()


MESSAGE( STATUS "CMAKE_BINARY_DIR:     " ${CMAKE_BINARY_DIR})
MESSAGE( STATUS "Pumi configuration file location: " ${SCOREC_CONFIG})

#this is just example code, do your own thing
#add_executable(mylibrary mylibrary.cpp)
add_executable(a1tri2 a1tri2.cc)
add_executable(a1tri2_su2 a1tri2_su2.cc)
add_executable(a1tri2_su2tri a1tri2_su2tri.cc)
add_executable(a1tri2_polar a1tri2_polar.cc)
add_executable(a1tri3_polar a1tri3_polar.cc)
add_executable(a1tet a1tet.cc)
add_executable(renderit renderit.cc)
add_executable(print_geo_coords print_geo_coords.cc)
add_executable(verify verify.cc)
#add_executable(zsplit zsplit.cc)

#for any targets that use PUMI, just use this command
#to it to include the right directories and link to all
#the scorec libraries automatically.
#we recommend PUBLIC if the target is a library and
#PRIVATE if the target is an executable
#target_link_libraries(mylibrary PUBLIC SCOREC::core)
target_link_libraries(a1tri2 PRIVATE SCOREC::core)
target_link_libraries(a1tri2_su2 PRIVATE SCOREC::core)
target_link_libraries(a1tri2_su2tri PRIVATE SCOREC::core)
target_link_libraries(a1tri2_polar PRIVATE SCOREC::core)
target_link_libraries(a1tri3_polar PRIVATE SCOREC::core)
target_link_libraries(a1tet PRIVATE SCOREC::core)
target_link_libraries(renderit PRIVATE SCOREC::core)
target_link_libraries(print_geo_coords PRIVATE SCOREC::core)
target_link_libraries(verify PRIVATE SCOREC::core)
#target_link_libraries(zsplit PRIVATE SCOREC::core)
install (TARGETS a1tri2 a1tri2_su2 a1tri2_su2tri a1tri2_polar a1tri3_polar a1tet renderit print_geo_coords verify DESTINATION bin)
