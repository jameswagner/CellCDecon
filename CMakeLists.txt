cmake_minimum_required(VERSION 3.10)
project(CellCDecon VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Library sources
set(LIB_SOURCES
    src/CellCDecon.cpp
    src/CellCDeconIO.cpp
)

# Include directories
include_directories(include)

# Create library
add_library(cellcdecon ${LIB_SOURCES})

# Create executable
add_executable(cellcdecon_app src/main.cpp)
target_link_libraries(cellcdecon_app PRIVATE cellcdecon)

# Options
option(BUILD_TESTS "Build the test programs" ON)

# Add tests
if(BUILD_TESTS)
    add_subdirectory(test)
endif()

# Install rules
install(TARGETS cellcdecon
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

install(TARGETS cellcdecon_app
    RUNTIME DESTINATION bin
)

install(FILES include/CellCDecon.h include/CellCDeconIO.h
    DESTINATION include/CellCDecon
) 