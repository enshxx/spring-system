rm -rf build
mkdir plots/csv_data
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .