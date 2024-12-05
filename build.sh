rm -rf build
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
cmake --build .
ctest --rerun-failed --output-on-failure .