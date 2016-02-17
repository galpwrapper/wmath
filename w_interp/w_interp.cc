#include <iostream>
#include <cstring>
#include <vector>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "pArray.h"
#include "w_interp.h"
using std::vector;                                                                                                                                             
using std::cout;                                                                
using std::endl;
using std::copy;

// in (line, line) space
#define line_do(x0, y0, x1, y1, x, y)\
  double sl = (y1 - y0)/(x1 - x0);\
  y = sl * (x - x0) + y0;

// in (line, log) space
#define log_do(x0, y0, x1, y1, x, y)\
  double sl = (x - x0)/(x1 - x0);\
  y = y0 * pow(y1/y0, sl);

// in (log, log) space
#define loglog_do(x0, y0, x1, y1, x, y)\
  if (x0 && y0 && x1 && y1) {\
    double sl = log(y1/y0)/log(x1/x0);\
    y = y0 * pow(x/x0, sl);\
  } else {\
    y = 0;\
  }


#define INTERP_1D_FUN(type)\
  template <typename T>\
    bool WMATH::interp_1D_##type(const T * Xorg, const T * Yorg, int Norg,\
                                 const T * Xout,       T * Yout, int Nout) {\
      if(Norg <2) {\
        cout<<"interp_1D_"<<#type<<" ERROR: input points < 2"<<endl;\
        return false;\
      }\
      int small_pre = 0;\
      for(int i=0; i<Nout; i++) {\
        int small, big;\
        if(Xout[i] <= Xorg[0]) {\
          small = 0; big = 1;\
        } else if (Xout[i] >= Xorg[Norg-1]) {\
          small = Norg-2; big = Norg-1;\
        } else {\
          for(int ii=small_pre; ii<Norg; ii++) {\
            if(Xout[i] < Xorg[ii]) {\
              big = ii; small = big - 1;\
              small_pre = small;\
              break;\
            }\
          }\
        }\
        if(Xout[i] == Xorg[small]) {\
          Yout[i] = Yorg[small];\
        } else {\
          type##_do(Xorg[small], Yorg[small], Xorg[big], Yorg[big], Xout[i], Yout[i])\
        }\
      }\
      return true;\
    }\
  bool WMATH::interp_1D_##type(const pArray & Xorg, const pArray & Yorg,\
                               const pArray & Xout,       pArray & Yout) {\
    if(Xorg.GetLength() != Yorg.GetLength() ||\
       Xout.GetLength() != Yout.GetLength()) {\
      cout<<"interp_1D_"<<#type<<" ERROR: length error"<<endl;\
      cout<<"(org)pArray::"<<Xorg.GetName()<<": Length = "<<Xorg.GetLength()<<endl;\
      cout<<"(org)pArray::"<<Yorg.GetName()<<": Length = "<<Yorg.GetLength()<<endl;\
      cout<<"(out)pArray::"<<Xout.GetName()<<": Length = "<<Xout.GetLength()<<endl;\
      cout<<"(out)pArray::"<<Yout.GetName()<<": Length = "<<Yout.GetLength()<<endl;\
      return false;\
    }\
    if(Xorg.GetUnit() != Xout.GetUnit() ||\
       Yorg.GetUnit() != Yout.GetUnit()) {\
      cout<<"interp_1D_"<<#type<<" ERROR: unit error"<<endl;\
      cout<<"(org)pArray::"<<Xorg.GetName()<<": Unit = "<<Xorg.GetLength()<<endl;\
      cout<<"(out)pArray::"<<Xout.GetName()<<": Unit = "<<Xout.GetLength()<<endl;\
      cout<<"(org)pArray::"<<Yorg.GetName()<<": Unit = "<<Yorg.GetLength()<<endl;\
      cout<<"(out)pArray::"<<Yout.GetName()<<": Unit = "<<Yout.GetLength()<<endl;\
      return false;\
    }\
    int Norg = Xorg.GetLength();\
    int Nout = Xout.GetLength();\
    bool success_or_not = interp_1D_##type(Xorg.a, Yorg.a, Norg,\
                                           Xout.a, Yout.a, Nout);\
    return success_or_not;\
  }\
  template <typename T>\
    bool WMATH::interp_1D_##type(const vector<T> & Xorg, const vector<T> & Yorg,\
                                 const vector<T> & Xout,       vector<T> & Yout) {\
      if(Xorg.size() != Yorg.size() ||\
         Xout.size() != Yout.size()) {\
        cout<<"interp_1D_"<<#type<<" ERROR: length error"<<endl;\
        cout<<"(org)vector<float>::X: Length = "<<Xorg.size()<<endl;\
        cout<<"(org)vector<float>::Y: Length = "<<Yorg.size()<<endl;\
        cout<<"(out)vector<float>::X: Length = "<<Xout.size()<<endl;\
        cout<<"(out)vector<float>::Y: Length = "<<Yout.size()<<endl;\
        return false;\
      }\
      int Norg = Xorg.size();\
      int Nout = Xout.size();\
      bool success_or_not = interp_1D_##type(Xorg.data(), Yorg.data(), Norg,\
                                             Xout.data(), Yout.data(), Nout);\
      return success_or_not;\
    }

INTERP_1D_FUN(line)
INTERP_1D_FUN(log)
INTERP_1D_FUN(loglog)

#undef line_do
#undef log_do
#undef loglog_do
#undef INTERP_1D_FUN


bool WMATH::interp_1D_gsl(const double * Xorg, const double * Yorg, int Norg,
                          const double * Xout,       double * Yout, int Nout,
                          const gsl_interp_type * type) {
  gsl_interp_accel * acc = gsl_interp_accel_alloc();
  gsl_spline * spline = gsl_spline_alloc(type, Norg);
  if(acc==0 || spline==0) {
    cout<<"interp_1D_gsl ERROR: gsl initialization fail"<<endl;
    return false;
  }
  int init_state = gsl_spline_init(spline, Xorg, Yorg, Norg);
  if(init_state != 0) {
    cout<<"interp_1D_gsl ERROR: gsl initialization fail"<<endl;
    return false;
  }
    
  int nErrors = 0;
  for(int i=0; i<Nout; i++) {
    int out_range = gsl_spline_eval_e(spline, Xout[i], acc, Yout+i);
    if(out_range) {
      *(Yout+i) = 0;
      nErrors++;
    }
  }
  if(nErrors > 5) {
    cout<<"interp_1D_gsl WARRNING: more than 5 points out of range"<<endl;
  }

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  return true;
}
bool WMATH::interp_1D_gsl(const pArray & Xorg, const pArray & Yorg,
                          const pArray & Xout,       pArray & Yout,
                          const gsl_interp_type * type) {
  if(Xorg.GetLength() != Yorg.GetLength() ||
     Xout.GetLength() != Yout.GetLength()) {
    cout<<"interp_1D_gsl ERROR: length error"<<endl;
    cout<<"(org)pArray::"<<Xorg.GetName()<<": Length = "<<Xorg.GetLength()<<endl;
    cout<<"(org)pArray::"<<Yorg.GetName()<<": Length = "<<Yorg.GetLength()<<endl;
    cout<<"(out)pArray::"<<Xout.GetName()<<": Length = "<<Xout.GetLength()<<endl;
    cout<<"(out)pArray::"<<Yout.GetName()<<": Length = "<<Yout.GetLength()<<endl;
    return false;
  }
  if(Xorg.GetUnit() != Xout.GetUnit() ||
     Yorg.GetUnit() != Yout.GetUnit()) {
    cout<<"interp_1D_gsl ERROR: unit error"<<endl;
    cout<<"(org)pArray::"<<Xorg.GetName()<<": Unit = "<<Xorg.GetLength()<<endl;
    cout<<"(out)pArray::"<<Xout.GetName()<<": Unit = "<<Xout.GetLength()<<endl;
    cout<<"(org)pArray::"<<Yorg.GetName()<<": Unit = "<<Yorg.GetLength()<<endl;
    cout<<"(out)pArray::"<<Yout.GetName()<<": Unit = "<<Yout.GetLength()<<endl;
    return false;
  }

  int Norg = Xorg.GetLength();
  int Nout = Xout.GetLength();
  bool success_or_not = interp_1D_gsl(Xorg.a, Yorg.a, Norg,
                                      Xout.a, Yout.a, Nout);
  return success_or_not;
}
bool WMATH::interp_1D_gsl(const vector<double> & Xorg, const vector<double> & Yorg,
                          const vector<double> & Xout,       vector<double> & Yout,
                          const gsl_interp_type * type) {
  if(Xorg.size() != Yorg.size() ||
     Xout.size() != Yout.size()) {
    cout<<"interp_1D_gsl ERROR: length error"<<endl;
    cout<<"(org)vector<float>::X: Length = "<<Xorg.size()<<endl;
    cout<<"(org)vector<float>::Y: Length = "<<Yorg.size()<<endl;
    cout<<"(out)vector<float>::X: Length = "<<Xout.size()<<endl;
    cout<<"(out)vector<float>::Y: Length = "<<Yout.size()<<endl;
    return false;
  }
  int Norg = Xorg.size();
  int Nout = Xout.size();
  bool success_or_not = interp_1D_gsl(Xorg.data(), Yorg.data(), Norg,
                                      Xout.data(), Yout.data(), Nout);
  return success_or_not;
}

// 显式实例化
#define SUPPORT_1D(type, T)\
  template bool WMATH::interp_1D_##type(const T *, const T *, int,\
                                        const T *,       T *, int);\
  template bool WMATH::interp_1D_##type(const vector<T> &, const vector<T> &,\
                                        const vector<T> &,       vector<T> &);
SUPPORT_1D(line, double)
SUPPORT_1D(line, float)
SUPPORT_1D(log, double)
SUPPORT_1D(log, float)
SUPPORT_1D(loglog, double)
SUPPORT_1D(loglog, float)
#undef SUPPORT_1D
