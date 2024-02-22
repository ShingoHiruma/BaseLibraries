project(BaseTest CXX)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${Basedir}/build)

add_executable(testBasic.out main.cpp)
if(CMAKE_CXX_COMPILER MATCHES "icpc")
 message("Intel Classical Compiler use in Example !") 
 target_compile_options(testBasic.out PUBLIC -std=c++14 -O3 -I${Basedir} -qopenmp -ipo)
 target_link_libraries(testBasic.out PUBLIC -qmkl -ipo BasicFuncs)
#
elseif(CMAKE_CXX_COMPILER MATCHES "icpx")
 message("Intel New Compiler use in Example !") 
 target_compile_options(testBasic.out PUBLIC -std=c++14 -O3 -I${Basedir} -qopenmp -ipo)
 target_link_libraries(testBasic.out PUBLIC -qmkl -ipo BasicFuncs)
#
else()
 target_compile_options(testBasic.out PUBLIC -std=c++14 -O3 -I${Basedir} -fopenmp)
 target_link_libraries(testBasic.out PUBLIC -lgomp BasicFuncs)
endif()
