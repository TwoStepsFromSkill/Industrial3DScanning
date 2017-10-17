cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(Industrial3DScanning)

set(CMAKE_CXX_STANDARD 11)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(warnings "-Wall -Wextra -pedantic")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(warnings "/W4 /EHsc")
endif()

if (NOT CONFIGURED_ONCE)
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${warnings}"
    CACHE STRING "Flags used by the compiler for all build types" FORCE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${warnings}"
    CACHE STRING "Flags used by the compiler for all build types" FORCE)
endif()

find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

add_definitions( ${QT_DEFINITIONS} )
include_directories( ${CMAKE_BINARY_DIR} )

include_directories(${Qt5Widgets_INCLUDES})
add_definitions(${Qt5Widgets_DEFINITIONS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")

# OpenGL & GLEW library
find_package(OpenGL REQUIRED)

set(HEADER)
set(SOURCE)
add_subdirectory(include)
add_subdirectory(source)

include_directories(include)

add_executable(Industrial3DScanning ${SOURCE} ${HEADER} ${MOC_SRCS})
qt5_use_modules(Industrial3DScanning Widgets OpenGL)
target_link_libraries(Industrial3DScanning ${QT_LIBRARIES} ${OPENGL_LIBRARIES})

set(CONFIGURED_ONCE TRUE CACHE INTERNAL "A flag showing that CMake has been configured at least once.")
