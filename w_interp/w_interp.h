#ifndef _W_INTERP_H
#define _W_INTERP_H
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "pArray.h"
using std::vector;

namespace WMATH {

// line-like interpolation support double and float
#define INTERP_1D(type)\
  template <typename T>\
    bool interp_1D_##type(const T *, const T *, int,\
                          const T *,       T *, int);\
  bool interp_1D_##type(const pArray &, const pArray &,\
                        const pArray &,       pArray &);\
  template <typename T>\
    bool interp_1D_##type(const vector<T> &, const vector<T> &,\
                          const vector<T> &,       vector<T> &);

INTERP_1D(line)
INTERP_1D(log)
INTERP_1D(loglog)
#undef INTERP_1D

// gsl 1D interpolation only support double
bool interp_1D_gsl(const double *, const double *, int,
                   const double *,       double *, int,
                   const gsl_interp_type * type = gsl_interp_steffen);
bool interp_1D_gsl(const pArray &, const pArray &,
                   const pArray &,       pArray &,
                   const gsl_interp_type * type = gsl_interp_steffen);
bool interp_1D_gsl(const vector<double> &, const vector<double> &,
                   const vector<double> &,       vector<double> &,
                   const gsl_interp_type * type = gsl_interp_steffen);
}
#endif
