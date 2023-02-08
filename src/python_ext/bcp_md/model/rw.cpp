#include "bcp_md/model/cpp_rw.h"
#include <pybind11/pybind11.h>

PYBIND11_MODULE(rw, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("fn_rw", &rw_test);
}
