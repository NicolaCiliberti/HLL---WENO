cmake_minimum_required(VERSION 3.26)
project(Projet)

set(CMAKE_CXX_STANDARD 17)

add_executable(Projet main.cpp
        Grid.h
        Grid.cpp
        Flux.h
        Flux.cpp
        Time.h
        helpers.cpp
        Time.cpp
        Print.cpp
        Grid1D.cpp
)
