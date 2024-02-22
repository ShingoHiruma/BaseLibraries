# CMakeのバージョンを設定
cmake_minimum_required(VERSION 3.10)
#コンパイラ変更
set(CMAKE_CXX_COMPILER icpc)
if(CMAKE_CXX_COMPILER MATCHES "icpc")
 message("Intel Compiler use !")
endif()
#
# プロジェクト名と使用する言語を設定
project(BaseLibs CXX)
set(Basedir ${CMAKE_CURRENT_SOURCE_DIR})
# 行列ライブラリを登録
add_subdirectory(SparseMat)
# サンプルを登録
add_subdirectory(examples/example_sparse_mat)
add_subdirectory(examples/example_abmc)
