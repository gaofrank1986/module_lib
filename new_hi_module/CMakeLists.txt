set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
add_library(hi2_mod SHARED geo2D.f90 param.f90 SIEPPEM.f90 matrix.f90 gaussian.f90 ./sgb/func.f90)

target_link_libraries(hi2_mod data)
#set_target_properties(hello PROPERTIES OUTPUT_NAME "hello")

