cmake_minimum_required (VERSION 2.8)
project (FINDCLUSTER)

include_directories (${FINDCLUSTER_SOURCE_DIR})
link_directories (${FINDCLUSTER_BINARY_DIR})

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "-O3 -Wall -DNDEBUG")        ## Optimize
# set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary

add_executable (findcluster.x findcluster.cpp findcluster_cluster.cpp findcluster_path.cpp findcluster_radius.cpp findcluster_box.cpp findcluster_writemeasure.cpp findcluster_helper.cpp findcluster_write.cpp)

# Link the executable to the Hello library.
# target_link_libraries (calcalljack calcint calcjack) 
