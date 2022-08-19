project(SparseMat CXX)
#出力先はベース
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${Basedir})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${Basedir})

#コンパイル設定
add_library(SparseMat
  STATIC
    SparseMat.cpp
    SparseMatC.cpp
    SparseMatOperators.cpp
    MatSolvers_ICCG.cpp
    MatSolvers_ICCOCG
    MatSolvers_Eigen.cpp
    MatSolvers_MKL.cpp
)
#オプションとライブラリリンクの設定
if(CMAKE_CXX_COMPILER MATCHES "icpc")
 message("Intel Compiler use in SparseMat !")
 target_compile_options(SparseMat PUBLIC -std=c++14 -O3 -xhost -I${Basedir} -qopenmp -ipo)
 target_link_libraries(SparseMat PUBLIC -qmkl)
else()
 target_compile_options(SparseMat PUBLIC -std=c++14 -O3 -I${Basedir} -fopenmp)
 target_link_libraries(SparseMat PUBLIC -lgomp)
endif()

