project(BasicFuncs CXX)
#蜃ｺ蜉帛・縺ｯ繝吶・繧ｹ
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${Basedir})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${Basedir})

#繧ｳ繝ｳ繝代う繝ｫ險ｭ螳・add_library(BasicFuncs
  STATIC
    CommonFuncs.cpp
    Mt.cpp
    WaveFormModifier.cpp
)

#繧ｪ繝励す繝ｧ繝ｳ縺ｨ繝ｩ繧､繝悶Λ繝ｪ繝ｪ繝ｳ繧ｯ縺ｮ險ｭ螳・if(CMAKE_CXX_COMPILER MATCHES "icpc")
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

