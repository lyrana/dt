# Make the Python module called dnt_cpp

# dnt_cpp: dnt_cpp.cpython-36m-x86_64-linux-gnu.so

dnt_cpp.cpython-36m-x86_64-linux-gnu.so: dnt_cpp.cpp predicates.o fstruct.o pstruct.o dolfin.o pseg.o
	c++ -O3 -Wall -shared -std=c++11 -fPIC -fvisibility=hidden -I/usr/include/eigen3 `python3 -m pybind11 --includes` dnt_cpp.cpp -o dnt_cpp`python3-config --extension-suffix` predicates.o fstruct.o pstruct.o dolfin.o pseg.o 

# Compile the source files

dolfin.o: dolfin.cpp dolfin.h fstruct.h pstruct.h predicates.h
	c++ -c -O3 -Wall -std=c++11 -fPIC -fvisibility=hidden -I/usr/include/eigen3 `python3 -m pybind11 --includes` dolfin.cpp -o dolfin.o

fstruct.o: fstruct.cpp fstruct.h
	c++ -c -O3 -Wall -std=c++11 -fPIC -fvisibility=hidden -I/usr/include/eigen3 `python3 -m pybind11 --includes` fstruct.cpp -o fstruct.o

predicates.o: predicates.cpp predicates.h
	c++ -c -O3 -Wall -std=c++11 -fPIC -fvisibility=hidden -I/usr/include/eigen3 `python3 -m pybind11 --includes` predicates.cpp -o predicates.o

pseg.o: pseg.cpp pseg.h
	c++ -c -O3 -Wall -std=c++11 -fPIC -fvisibility=hidden -I/usr/include/eigen3 `python3 -m pybind11 --includes` pseg.cpp -o pseg.o

pstruct.o: pstruct.cpp pstruct.h
	c++ -c -O3 -Wall -std=c++11 -fPIC -fvisibility=hidden -I/usr/include/eigen3 `python3 -m pybind11 --includes` pstruct.cpp -o pstruct.o

clean:
	/bin/rm *.so *.o
