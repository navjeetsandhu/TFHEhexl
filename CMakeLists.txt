cmake_minimum_required(VERSION 3.16)
project(tfhe++ CXX)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS
      "-march=native -O3 -g -funroll-loops -Wall -Wextra -pedantic -Wno-sign-compare"
  )
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

option(USE_HEXL "Use Intel HEXL" ON)
option(USE_RANDEN "Use randen as CSPRNG" ON)
option(ENABLE_TEST "Build tests" ON)


set(TFHEpp_DEFINITIONS
    ""
    PARENT_SCOPE)

if(USE_RANDEN)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_RANDEN"
      PARENT_SCOPE)
  add_compile_definitions(USE_RANDEN)
  add_subdirectory(thirdparties/randen)
endif()

if(USE_HEXL)
  set(hexl-ntt_DEFINITIONS
          "${hexl-ntt_DEFINITIONS};USE_HEXL"
          PARENT_SCOPE)
  add_compile_definitions(USE_HEXL)
  add_subdirectory(thirdparties/hexl)
  # set(CMAKE_CXX_FLAGS "-march=native -O3 -g -funroll-loops -Wall -Wextra
  # -pedantic -Wno-sign-compare -mprefer-vector-width=512" )
endif()


add_subdirectory(thirdparties/spqlios)

# For OpenMP
find_package(OpenMP)
if(OpenMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


add_subdirectory(src)

if(ENABLE_TEST)
  add_subdirectory(unit_test/other)
  add_subdirectory(unit_test/ntt)
  add_subdirectory(unit_test/tutorial)
  add_subdirectory(unit_test/encrypt)
  add_subdirectory(unit_test/gate)
  add_subdirectory(unit_test/externalproduct)
  add_subdirectory(unit_test/decomposition)
  add_subdirectory(unit_test/bootstrapping)
endif()

install(TARGETS tfhe++ LIBRARY DESTINATION lib)
