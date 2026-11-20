cmake_minimum_required(VERSION 3.20)
project(nfft3_consumer LANGUAGES C)

find_package(NFFT3 REQUIRED)

add_executable(use_nfft3 use_nfft3.c)
target_link_libraries(use_nfft3 PRIVATE NFFT3::nfft3)

add_executable(use_nfft3_nocomplex use_nfft3_nocomplex.c)
target_compile_definitions(use_nfft3_nocomplex PRIVATE NFFT_PRECISION_DOUBLE)
target_link_libraries(use_nfft3_nocomplex PRIVATE NFFT3::nfft3)
