# Find GTest package
find_package(GTest REQUIRED)

# Include directories
include_directories(${GTEST_INCLUDE_DIRS})
include_directories(../include)
include_directories(unit)

# Define test sources
set(TEST_SOURCES
    unit/test_cellc.cpp
    unit/test_cellc_utils.cpp
)

# Create the test executable
add_executable(cellcdecon_tests ${TEST_SOURCES})

# Link libraries
target_link_libraries(cellcdecon_tests
    PRIVATE
    cellcdecon
    ${GTEST_LIBRARIES}
    pthread
)

# Add tests
include(GoogleTest)
gtest_discover_tests(cellcdecon_tests)

# Copy test data for tests
file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}) 