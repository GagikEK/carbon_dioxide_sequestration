cmake_minimum_required(VERSION 3.10)
project(Exercice3 LANGUAGES CXX)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# collect sources
file(GLOB_RECURSE SRC_FILES
	"${CMAKE_SOURCE_DIR}/scripts/*.cpp"
)

add_executable(mon_exe ${SRC_FILES})

target_include_directories(mon_exe PRIVATE "${CMAKE_SOURCE_DIR}/scripts")
target_link_libraries(mon_exe PRIVATE Eigen3::Eigen)
