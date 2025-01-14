mkdir build
cd build
cmake ..
make -j16
cd ..
cd run
cp ../src/salt/base/flute/POST9.dat ./
cp ../src/salt/base/flute/POWV9.dat ./