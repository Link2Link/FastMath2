cmake_minimum_required(VERSION 3.16)

project(FastMath2 VERSION 0.2.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#if (CMAKE_COMPILER_IS_GNUCXX)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")
#elseif (MSVC)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  /std:c++17")
#else ()
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++17 -O3")
#endif (CMAKE_COMPILER_IS_GNUCXX)


set(PROJECT_TOP_DIR ${PROJECT_SOURCE_DIR})
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/bin)


#add_library(FastMath INTERFACE)
#target_include_directories(FastMath PUBLIC src)

option(ENABLE_UNIT_TESTS "Enable unit tests" ON)
option(ENABLE_DEVELOPMENT "Enable development" ON)
option(ENABLE_EXAMPLE "Enable example" ON)

if(${ENABLE_UNIT_TESTS})
    add_definitions(-DENABLE_UNIT_TESTS)
    add_subdirectory(test)
endif()

if(${ENABLE_DEVELOPMENT})
    add_executable(develop develop.cpp)
    target_include_directories(develop PUBLIC ${PROJECT_TOP_DIR}/FastMath)
endif()

if (${ENABLE_EXAMPLE})
    SET(TEST_INSTALL_DIR ${PROJECT_TOP_DIR}/bin/examples)
    FILE(GLOB EXAMPLE_SRC "${PROJECT_SOURCE_DIR}/examples/*.cpp")
    foreach(item IN LISTS EXAMPLE_SRC)
        get_filename_component(SRCNAME ${item} NAME_WE) #得到文件的名字，不带扩展名
        add_executable(${SRCNAME} ${item})
        target_include_directories(${SRCNAME} PUBLIC ${PROJECT_TOP_DIR}/FastMath)
        install(TARGETS ${SRCNAME} DESTINATION ${TEST_INSTALL_DIR})
        message(STATUS "Example : ${SRCNAME} <- ${item}")

    endforeach ()


endif ()

