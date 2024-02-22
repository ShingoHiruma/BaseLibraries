
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <memory>
#include <iostream>

#include "SparseMat.hpp"
#include "MatSolvers.hpp"

using namespace std;
using namespace SRLfem;
namespace py = pybind11;

PYBIND11_MODULE(SparseMatPy, m)
{
  py::class_<SparseMat, shared_ptr<SparseMat>>(m, "SparseMat")
      .def(py::init<>([](const py::array_t<int> &rows,
                         const py::array_t<int> &cols,
                         const py::array_t<double> &vals,
                         int N)
                      {
      auto mat = make_shared<SparseMat>(N);
      for(int i=0; i<rows.size(); i++)
      {
        mat->add(rows.at(i), cols.at(i), vals.at(i));
      }
      mat->fix();
      return mat; }));

  m.def("ICCG", [](const SparseMat &A,
                   const py::array_t<double> &b,
                   py::array_t<double> &x,
                   int max_iter, double tol, double acc)
        {
    MatSolvers::solveICCG(b.size(), tol, max_iter, acc, A, b.data(), const_cast<double*>(x.data()));
    return x; });
}
