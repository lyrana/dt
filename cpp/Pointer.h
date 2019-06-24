#include <pybind11/pybind11.h>
#include <memory>
#include <iostream>

namespace py = pybind11;

//! PointerWrapper is a class that contains just a pointer to an object of type T
/*!

*/

template <typename T> class PointerWrapper
{
 public:
 // The constructors
  // Null argument
  PointerWrapper() : ptr(nullptr) {}
  // A pointer-to-T argument
  PointerWrapper(T* ptr) : ptr(ptr) {}
  // A PointerWrapper object as argument
  PointerWrapper(const PointerWrapper& other) : ptr(other.ptr) {}

 // The member functions
  // The * op
  T& operator* () const { return *ptr; }
  // The -> op
  T* operator->() const { return  ptr; }
  // The get() function
  T* get() const { return ptr; }
  // The destroy() function
  void destroy() { delete ptr; }
  // The array item [] op
  T& operator[](std::size_t idx) const { return ptr[idx]; }

  
 private:
  // The pointer-to-T
  T* ptr;
};

float array[3] = { 3.14, 2.18, -1 };

//! Define get_ptr()
/*!
    This returns an object of type 

    It calls a constructor with argument "array".
 */
PointerWrapper<float> get_ptr(void) { return array; }

//! Define use_ptr()
void use_ptr(PointerWrapper<float> ptr) {
  for (int i = 0; i < 3; ++i)
    std::cout << ptr[i] << " ";
  std::cout << "\n";
}

PYBIND11_MODULE(Ptr,m)
{
    py::class_<PointerWrapper<float>>(m,"pfloat");
    m.def("get_ptr", &get_ptr);
    m.def("use_ptr", &use_ptr);
}
