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
    MatSolvers_ICCOCG.cpp
    MatSolvers_Eigen.cpp
    MatSolvers_MKL.cpp
    MatSolvers_ABMCICCG.cpp
)
#オプションとライブラリリンクの設定
if(CMAKE_CXX_COMPILER MATCHES "icpc")
 message("Intel Classical Compiler use in SparseMat !") 
 target_compile_options(SparseMat PUBLIC -std=c++14 -O3 -I${Basedir} -qopenmp -ipo)
 target_link_libraries(SparseMat PUBLIC -qmkl)
#
elseif(CMAKE_CXX_COMPILER MATCHES "icpx")
 message("Intel New Compiler use in SparseMat !") 
 set(CMAKE_AR "/opt/intel/oneapi/compiler/latest/linux/bin-llvm/llvm-ar")
 set(CMAKE_NM "/opt/intel/oneapi/compiler/latest/linux/bin-llvm/llvm-nm")
 set(CMAKE_RANLIB "/opt/intel/oneapi/compiler/latest/linux/bin-llvm/llvm-ranlib")
 target_compile_options(SparseMat PUBLIC -std=c++14 -O3 -I${Basedir} -qopenmp -ipo)
 target_link_libraries(SparseMat PUBLIC -qmkl)
#
else()
  find_package(OpenMP REQUIRED)
  target_include_directories(SparseMat PUBLIC ${Basedir})
  target_compile_options(SparseMat PUBLIC -std=c++14 -O3 -fopenmp -fPIC)
  target_link_libraries(SparseMat PUBLIC OpenMP::OpenMP_CXX)
endif()

if(PYBIND_EXPORT)
  find_package(pybind11 REQUIRED)
  pybind11_add_module(SparseMatPy MODULE SparseMatPy.cpp)
  target_link_libraries(SparseMatPy PUBLIC SparseMat)
endif()
