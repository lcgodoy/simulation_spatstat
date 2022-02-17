# File : ---- CMakeLists.txt
# Date : 23/12/2011
# Author : Tristan Lecorgne
# INRIA Serpico team
#
# Description : This file can be used in order to prepare the build of a cimg project
#
# Usage :
#    - Read all the comments
#    - Put your preferences in the line that ends with a "# EDIT"
#

#-----------------------------------------------------------------------------
### Paths

#
# If your system use atypical developpement folders, uncomment and set the following variables.
# This is useful on Windows systems in particular.
# Note: \ needs to be escaped (example C:\\include)
#

#set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "your include path") # EDIT 
#set (CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "your library path") # EDIT

#-----------------------------------------------------------------------------
### Check CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0 FATAL_ERROR)
 
#-----------------------------------------------------------------------------
### Project name
project(GcoPS) # EDIT : put your project name

#-----------------------------------------------------------------------------
### Add project directories to include path.
include_directories (
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

##-----------------------------------------------------------------------------
### Add our own CMake modules to the find_package search path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/CMakeModules")


#-----------------------------------------------------------------------------
### Compilation

# default compilation is release
if (NOT CMAKE_BUILD_TYPE MATCHES DEBUG)
    set (CMAKE_BUILD_TYPE "RELEASE")
endif (NOT CMAKE_BUILD_TYPE MATCHES DEBUG)

message ("")
message ("Configuring ${PROJECT_NAME} in ${CMAKE_BUILD_TYPE} mode")

# warning
set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -W -ansi -pedantic")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -W -ansi -pedantic")
set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -w")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -w")

# optimisation
if (MINGW OR CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
endif (MINGW OR CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
if (MSVC OR MSVC_IDE OR MSVC60 OR  MSVC70 OR MSVC71 OR MSVC80 OR CMAKE_COMPILER_2005 OR MSVC90 OR MSVC10)
    set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")        
endif (MSVC OR MSVC_IDE OR MSVC60 OR  MSVC70 OR MSVC71 OR MSVC80 OR CMAKE_COMPILER_2005 OR MSVC90 OR MSVC10)

# debug
set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Dcimg_debug=3")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Dcimg_debug=3")


#-----------------------------------------------------------------------------
### CIMG Support

#
# CIMG is very flexible and makes use of optional third party libraries.
# The following lines allow you to select which libraries you need to compile your project.
#
# When set to ON, cmake will search for the lib on the system and crash if it's absent. 
#   If present, the program will be linked to this library
#
# When set to OFF, cmake will ignore this library and the program will not use its functionnalities.
#   In consequence, the part of the code in CImg.h within the macro "#ifdef cimg_use_xxx" will not be compiled.
#

option (USE_DISPLAY_LIBS "Set to off if your application does not use any display (=> cimg_diplay=0)" OFF) # EDIT
option (USE_JPEG_LIBS "Set to off if you don't need libjpeg (=> cimg_use_jpeg undefined)" OFF) # EDIT
option (USE_PNG_LIBS "Set to off if you don't need libpng (=> cimg_use_png undefined)" OFF) # EDIT
option (USE_TIFF_LIBS "Set to off if you don't need libtiff (=> cimg_use_tiff undefined)" ON) # EDIT
option (USE_LAPACK_LIBS "Set to off if you don't need lapack libraries (=> cimg_use_lapack undefined)" OFF) # EDIT
option (USE_BOARD_LIBS "Set to off if you don't need libboard (=> cimg_use_board undefined)" OFF) # EDIT
option (USE_FFMPEG_LIBS "Set to off if you don't need libffmpeg (=> cimg_use_ffmpeg undefined)" OFF) # EDIT
option (USE_FFTW3_LIBS "Set to off if you don't need libfftw3 (=> cimg_use_fftw3 undefined)" ON) # EDIT
option (USE_MAGICK_LIBS "Set to off if you don't need libmagick++ (=> cimg_use_magick undefined)" OFF) # EDIT

#
# Call findCIMG.cmake and use the variables defined by it
#

# FFTW3
find_package (Fftw3ForCImgOnMac)

# CIMG
find_package (CIMG)
add_definitions (${CIMG_DEFINITIONS})
add_definitions (-Dcimg_use_vt100)
add_definitions (-Dcimg_use_fftw3)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CIMG_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CIMG_CXX_FLAGS}")

#-----------------------------------------------------------------------------
### Add project directories to include path.
include_directories (
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CIMG_INCLUDE_DIRS}
  ${FFTW3_INCLUDE_DIR}
)


#-----------------------------------------------------------------------------
### Executables

#
# Here you can declare what executable you want to compile with cmake commands. 
# For example :
# - add_executable (myExecutable myExecutable.cpp myExecutable.h otherRequiredFile.h ...)
# - add_library (myLibrary myLibrary.cpp ...)
# - target_link_libraries (myExecutable myLibrary  ...)
#
# Then, link ${CIMG_LIBRARIES} with your executables :
#
# - target_link_libraries (myExecutable ${CIMG_LIBRARIES})
#
# For more information on these commands and on other possible commands, go to cmake.org
#

add_executable(GcoPS colocalizationTest.cpp vesicleDisplaying2D.cpp vesicleDisplaying3D.cpp)
target_link_libraries (GcoPS ${CIMG_LIBRARIES} ${FFTW3_LIBRARY} ${FFTW3_THREADS_LIBRARY})

#-----------------------------------------------------------------------------
### Installation

#
# Here you can specify which programs and libraries should be copied on the system when user run "make install"
# Uncomment the following INSTALL command and replace myExecutable and myLibrary by the names of your application binaries
#

INSTALL(TARGETS GcoPS
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)


