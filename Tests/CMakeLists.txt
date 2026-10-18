# Native build of the engine-independent code (calculator, geometry, switch and sound models) and
# its tests. Unreal compiles the same .cpp files inside the game module; this build exists so the
# math can be checked without the engine:  cmake -S Tests -B build && cmake --build build && ctest --test-dir build
cmake_minimum_required(VERSION 3.16)
project(MechanicalCalculatorTests CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(MODULE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../Source/MechanicalCalculator)
file(GLOB CORE_SOURCES ${MODULE_DIR}/Calc/*.cpp ${MODULE_DIR}/Geometry/*.cpp ${MODULE_DIR}/Keyboard/*.cpp)
file(GLOB TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# -DMCALC_ONLY_LIB=Orbits builds just that library file (the others register nothing), for working
# on one library while the rest are being changed. Pair it with MCALC_CASES=orbits.tsv at run time.
set(MCALC_ONLY_LIB "" CACHE STRING "Build only Calc/Lib<name>.cpp")
if(MCALC_ONLY_LIB)
  list(FILTER CORE_SOURCES EXCLUDE REGEX ".*/Calc/Lib[A-Za-z]+\\.cpp$")
  if(EXISTS ${MODULE_DIR}/Calc/Lib${MCALC_ONLY_LIB}.cpp)   # "None" builds the core with no library
    list(APPEND CORE_SOURCES ${MODULE_DIR}/Calc/Lib${MCALC_ONLY_LIB}.cpp)
  endif()
  list(FILTER TEST_SOURCES EXCLUDE REGEX ".*/TestLib[A-Za-z]+\\.cpp$")
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/TestLib${MCALC_ONLY_LIB}.cpp)
    list(APPEND TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/TestLib${MCALC_ONLY_LIB}.cpp)
  endif()
  set(STUBS "#include \"CalcLibrary.h\"\nnamespace mcalc\n{\n")
  foreach(LIB Control Flight Radar Orbits Ops Atlas Fields)
    if(NOT LIB STREQUAL MCALC_ONLY_LIB)
      string(APPEND STUBS "void Register${LIB}Library(Registry& R) { (void)R; }\n")
    endif()
  endforeach()
  string(APPEND STUBS "}\n")
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lib_stubs.cpp "${STUBS}")
  list(APPEND CORE_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/lib_stubs.cpp)
endif()

# Unreal's constraints: no exceptions, no RTTI, shadowing is an error; plus the macro poison header.
set(STRICT_FLAGS -Wall -Wextra -Wshadow -Werror -fno-exceptions -fno-rtti
    -include ${CMAKE_CURRENT_SOURCE_DIR}/ue_macro_poison.h)

add_library(mcalc_core STATIC ${CORE_SOURCES})
target_include_directories(mcalc_core PUBLIC ${MODULE_DIR}/Calc ${MODULE_DIR}/Geometry ${MODULE_DIR}/Keyboard)
target_compile_options(mcalc_core PRIVATE ${STRICT_FLAGS})

add_executable(mcalc_tests ${TEST_SOURCES})
target_link_libraries(mcalc_tests PRIVATE mcalc_core)
target_compile_options(mcalc_tests PRIVATE -Wall -Wextra -Werror)
target_compile_definitions(mcalc_tests PRIVATE MCALC_CASES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/reference")

# Unity build: every core .cpp in one translation unit, as Unreal Build Tool does by default,
# so helper names that collide across files fail here.
set(UNITY_FILE ${CMAKE_CURRENT_BINARY_DIR}/unity_all.cpp)
set(UNITY_CONTENT "")
foreach(SRC ${CORE_SOURCES})
  string(APPEND UNITY_CONTENT "#include \"${SRC}\"\n")
endforeach()
file(WRITE ${UNITY_FILE}.in "${UNITY_CONTENT}")
configure_file(${UNITY_FILE}.in ${UNITY_FILE} COPYONLY)
add_library(mcalc_unity OBJECT ${UNITY_FILE})
target_include_directories(mcalc_unity PRIVATE ${MODULE_DIR}/Calc ${MODULE_DIR}/Geometry ${MODULE_DIR}/Keyboard)
target_compile_options(mcalc_unity PRIVATE ${STRICT_FLAGS})

# The calculator on the command line.
add_executable(mcalc_cli ${CMAKE_CURRENT_SOURCE_DIR}/tools/mcalc_cli.cpp)
target_link_libraries(mcalc_cli PRIVATE mcalc_core)

enable_testing()
add_test(NAME mcalc_tests COMMAND mcalc_tests)

# Preview: writes the assembled board as OBJ for Tests/preview/render.mjs (not part of the tests).
add_executable(mcalc_export_obj ${CMAKE_CURRENT_SOURCE_DIR}/preview/export_obj.cpp)
target_link_libraries(mcalc_export_obj PRIVATE mcalc_core)
target_compile_options(mcalc_export_obj PRIVATE -Wall -Wextra -Werror)
