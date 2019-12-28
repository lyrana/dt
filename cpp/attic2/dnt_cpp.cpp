#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_ex1(py::module &);
void init_ex2(py::module &);
/* ... */
PYBIND11_MODULE(dnt_cpp, m) {
init_ex1(m);
init_ex2(m);
/* ... */
}
