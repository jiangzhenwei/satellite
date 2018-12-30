# tests/CMakeLists.txt

# PCE
add_executable(gpsClock gpsClock.cpp)
target_link_libraries(gpsClock rocket)

# POD
add_executable(gpsOrbit gpsOrbit.cpp)
target_link_libraries(gpsOrbit rocket)

# others

