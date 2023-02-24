#include "bcp_md/model/cpp_rw.h"
#include <pybind11/pybind11.h>

PYBIND11_MODULE(rw, m) {
    m.doc() = R"pbdoc(
        random walk module
        ------------------
        .. currentmodule:: rw
        
        .. autosummary::
           :toctree: _generate

           fn_rw
        
        )pbdoc"; // optional module docstring

    m.def("fn_rw", &rw_test, R"pbdoc(
        self-avoiding random walk

        read in parameters from in.parameters, which is written by a python wrapper model.Data.generate
    )pbdoc");
}
