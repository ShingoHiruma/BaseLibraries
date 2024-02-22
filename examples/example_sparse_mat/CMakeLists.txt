project(MatTest CXX)
#出力先はベース
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${Basedir}/build)

#コンパイル設定
add_executable(testMat.out main.cpp)
#オプションとライブラリリンクの設定
if(CMAKE_CXX_COMPILER MATCHES "icpc")
 message("Intel Classical Compiler use in Example !") 
 target_compile_options(testMat.out PUBLIC -std=c++14 -O3 -I${Basedir} -qopenmp -ipo)
 target_link_libraries(testMat.out PUBLIC -qmkl -ipo SparseMat)
#
elseif(CMAKE_CXX_COMPILER MATCHES "icpx")
 message("Intel New Compiler use in Example !") 
 target_compile_options(testMat.out PUBLIC -std=c++14 -O3 -I${Basedir} -qopenmp -ipo)
 target_link_libraries(testMat.out PUBLIC -qmkl -ipo SparseMat)
#
else()
 target_compile_options(testMat.out PUBLIC -std=c++14 -O3 -I${Basedir} -fopenmp)
 target_link_libraries(testMat.out PUBLIC -lgomp SparseMat)
endif()
