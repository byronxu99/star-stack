include_directories(.)
file(GLOB src_files "*.cpp")
add_library(image_utils OBJECT ${src_files})