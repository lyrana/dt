# Make the Python module called dnt_cpp

# dnt_cpp: dnt_cpp.cpython-36m-x86_64-linux-gnu.so

dnt_cpp.cpython-36m-x86_64-linux-gnu.so: prec.o pseg.o
	c++ -O3 -Wall -shared -std=c++11 -fPIC -I/usr/include/eigen3 `python3 -m pybind11 --includes` dnt_cpp.cpp -o dnt_cpp`python3-config --extension-suffix` prec.o pseg.o

# Compile the source files

prec.o: prec.cpp
	c++ -c -O3 -Wall -std=c++11 -fPIC -I/usr/include/eigen3 `python3 -m pybind11 --includes` prec.cpp -o prec.o

pseg.o: pseg.cpp
	c++ -c -O3 -Wall -std=c++11 -fPIC -I/usr/include/eigen3 `python3 -m pybind11 --includes` pseg.cpp -o pseg.o
