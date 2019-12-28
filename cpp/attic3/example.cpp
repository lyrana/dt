#include <pybind11/pybind11.h>

namespace py = pybind11;

int add(int i, int j) {
  return i + j;
}

// Create a variable m of type py::module. This will make a file called example.so which
// can be imported using 'import example':
PYBIND11_MODULE(example, m) {
  
  m.doc() = "This is my first pybind11 example plugin"; // optional module docstring

// Create binding code to call add() from Python using module::def()
  m.def("add", &add, "A function which adds two numbers");
}

