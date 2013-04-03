cmake_minimum_required (VERSION 2.8)
project (POLLHISTO)

# Create a library called "Hello" which includes the source file "hello.cxx".
# The extension is already found. Any number of sources could be listed here.
# add_library (neib neib.cpp)
# add_library (chargeconfig chargeconfig.cpp

include_directories (${POLLHISTO_SOURCE_DIR})

# Make sure the linker can find the Hello library once it is built.
link_directories (${POLLHISTO_BINARY_DIR})

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "-O3")        ## Optimize
# set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary

add_executable (bin/pollhisto.x pollhisto.cpp)
add_executable (bin/genpollhisto.x genpollhisto.cpp)

# Link the executable to the Hello library.
# target_link_libraries (calcalljack calcint calcjack) 
