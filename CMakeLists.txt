cmake_minimum_required(VERSION 3.13)
project(ImageCompletion)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenCV)

include_directories(${OpenCV_INCLUDE_DIRS})

add_executable(ImageCompletion main.cpp OpenCvUtility.cpp OpenCvUtility.h PointManager.cpp PointManager.h
        StructurePropagation.cpp StructurePropagation.h Photometric.cpp Photometric.h TextureCompletion.cpp
        TextureCompletion.h)
target_link_libraries(ImageCompletion ${OpenCV_LIBS})