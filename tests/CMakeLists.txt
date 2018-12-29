# tests/CMakeLists.txt

# OI
add_executable(gpsOrbitFit gpsOrbitFit.cpp)
target_link_libraries(gpsOrbitFit rocket)

add_executable(gpsOI gpsOI.cpp)
target_link_libraries(gpsOI rocket)

# RINEX
add_executable(recClock recClock.cpp)
target_link_libraries(recClock rocket)

add_executable(satClock satClock.cpp)
target_link_libraries(satClock rocket)

# PCE
add_executable(gpsClockDualOld gpsClockDualOld.cpp)
target_link_libraries(gpsClockDualOld rocket)

add_executable(gpsClockDualNew gpsClockDualNew.cpp)
target_link_libraries(gpsClockDualNew rocket)

# POD
add_executable(gpsOrbitDualOld gpsOrbitDualOld.cpp)
target_link_libraries(gpsOrbitDualOld rocket)

add_executable(gpsOrbitDualNew gpsOrbitDualNew.cpp)
target_link_libraries(gpsOrbitDualNew rocket)

# others

