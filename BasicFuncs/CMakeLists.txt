project(BasicFuncs CXX)
#出力先はベース
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${Basedir})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${Basedir})

#コンパイル設定
add_library(BasicFuncs
  STATIC
    CommonFuncs.cpp
    Mt.cpp
    WaveFormModifier.cpp
)

#オプションとライブラリリンクの設定
if(CMAKE_CXX_COMPILER MATCHES "icpc")
 message("Intel Classical Compiler use in SparseMat !") 
 target_compile_options(BasicFuncs PUBLIC -std=c++14 -O3 -I${Basedir} -qopenmp -ipo)
 target_link_libraries(BasicFuncs PUBLIC -qmkl)
#
elseif(CMAKE_CXX_COMPILER MATCHES "icpx")
 message("Intel New Compiler use in SparseMat !") 
 set(CMAKE_AR "/opt/intel/oneapi/compiler/latest/linux/bin-llvm/llvm-ar")
 set(CMAKE_NM "/opt/intel/oneapi/compiler/latest/linux/bin-llvm/llvm-nm")
 set(CMAKE_RANLIB "/opt/intel/oneapi/compiler/latest/linux/bin-llvm/llvm-ranlib")
 target_compile_options(BasicFuncs PUBLIC -std=c++14 -O3 -I${Basedir} -qopenmp -ipo)
 target_link_libraries(BasicFuncs PUBLIC -qmkl)
#
else()
 target_compile_options(BasicFuncs PUBLIC -std=c++14 -O3 -I${Basedir} -fopenmp)
 target_link_libraries(BasicFuncs PUBLIC -lgomp)
endif()

