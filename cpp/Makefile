# Make the Python module called dnt_cpp

# dnt_cpp: dnt_cpp.cpython-36m-x86_64-linux-gnu.so

dnt_cpp.cpython-36m-x86_64-linux-gnu.so: dnt_cpp.cpp predicates.o fstruct.o pstruct.o dolfin_functions.o pseg.o SegmentedArrayPair_Cpp.h
	c++ -O3 -Wall -shared -std=c++11 -fPIC -fvisibility=hidden -I/usr/include/eigen3 `python3 -m pybind11 --includes` dnt_cpp.cpp -o dnt_cpp`python3-config --extension-suffix` predicates.o fstruct.o pstruct.o dolfin_functions.o pseg.o SegmentedArrayPair_Cpp.h

# Compile the source files

dolfin_functions.o: dolfin_functions.cpp dolfin_functions.h fstruct.h pstruct.h predicates.h
	c++ -c -O3 -Wall -std=c++17 -fPIC -fvisibility=hidden -I/usr/include/eigen3 `python3 -m pybind11 --includes` dolfin_functions.cpp -o dolfin_functions.o

fstruct.o: fstruct.cpp fstruct.h
	c++ -c -O3 -Wall -std=c++11 -fPIC -fvisibility=hidden -I/usr/include/eigen3 `python3 -m pybind11 --includes` fstruct.cpp -o fstruct.o

predicates.o: predicates.cpp predicates.h
	c++ -c -O3 -Wall -std=c++11 -fPIC -fvisibility=hidden -I/usr/include/eigen3 `python3 -m pybind11 --includes` predicates.cpp -o predicates.o

pseg.o: pseg.cpp pseg.h
	c++ -c -O3 -Wall -std=c++11 -fPIC -fvisibility=hidden -I/usr/include/eigen3 `python3 -m pybind11 --includes` pseg.cpp -o pseg.o

pstruct.o: pstruct.cpp pstruct.h
	c++ -c -O3 -Wall -std=c++11 -fPIC -fvisibility=hidden -I/usr/include/eigen3 `python3 -m pybind11 --includes` pstruct.cpp -o pstruct.o

#SegmentedArrayPair_Cpp.o: SegmentedArrayPair_Cpp.cpp SegmentedArrayPair_Cpp.h pstruct.h
#	c++ -c -O3 -Wall -std=c++11 -fPIC -fvisibility=hidden `python3 -m pybind11 --includes` SegmentedArrayPair_Cpp.cpp -o SegmentedArrayPair_Cpp.o

#test_segmented_arrays: test_segmented_arrays.cpp SegmentedArrayPair_Cpp.cpp
#	c++ -O0 -g -Wall -std=c++11 -fPIC -fvisibility=hidden `python3-config --cflags` test_segmented_arrays.cpp -o test_segmented_arrays `python3-config --ldflags`

clean:
	/bin/rm *.so *.o
