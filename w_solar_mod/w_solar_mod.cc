#include <iostream>
#include <algorithm>    // std::copy
#include <cmath>
#include "w_solar_mod.h"
#include "w_interp.h"
using std::cout;
using std::endl;
using std::copy;

const double WMATH::solar_mod::eMass = 0.5109990615e-3;
const double WMATH::solar_mod::pMass = 939.e-3;

WMATH::solar_mod::solar_mod() {
  init(0, 0, 0);
}
WMATH::solar_mod::solar_mod(int _Z,
                            int _A,
                            double _phi) {
  init(_Z, _A, _phi);
}
WMATH::solar_mod & WMATH::solar_mod::set(int _Z,
                                         int _A,
                                         double _phi) {
  init(_Z, _A, _phi);
  return *this;
}
void WMATH::solar_mod::init(int _Z,
                            int _A,
                            double _phi) {
  Z = _Z;
  _A == 0 ? A = 1, mass = eMass : A = _A, mass = pMass;
  phi = _phi;
  if (phi < 0) {
    cout<<"[init_ERROR]: phi < 0"<<endl;
  }
}

// modulate
template <typename T>
  bool WMATH::solar_mod::modulate(const T * Eorg, const T * Forg, int Norg,
                                  const T * Eout,       T * Fout, int Nout) {
    if (phi < 0)
      return false;
    else if (phi == 0)
      return interp_1D_loglog(Eorg, Forg, Norg, Eout, Fout, Nout);

    // for phi > 0
    T * Etemp = new T[Nout];
    for (int i=0; i<Nout; i++) {
      Etemp[i] = Eout[i] + fabs(Z)*phi/A;
    }
    if (!interp_1D_loglog(Eorg, Forg, Norg, Etemp, Fout, Nout)) {
      cout<<"[modulate_ERROR]: unsuccess interp_1D_loglog"<<endl;
      delete [] Etemp;
      return false;
    }
    for (int i=0; i<Nout; i++) {
      Fout[i] *= Eout[i]*(Eout[i]+2*mass) / (Etemp[i]*(Etemp[i]+2*mass));
    }
    // modulate success
    delete [] Etemp;
    return true;
  }
template <typename T>
  bool WMATH::solar_mod::modulate(const T * Eorg,       T *Forg, int Norg) {
    T * Fout = new T[Norg];
    bool success_or_not = modulate(Eorg, Forg, Norg, Eorg, Fout, Norg);
    copy(Fout, Fout+Norg, Forg);
    delete [] Fout;
    return success_or_not;
  }
// for vector input
template <typename T>
  bool WMATH::solar_mod::modulate(const vector<T> & Eorg, const vector<T> & Forg,
                                  const vector<T> & Eout,       vector<T> & Fout) {
    if (Eorg.size() != Forg.size()) {
      cout<<"[modulate_ERROR](vector input): illegal length"<<endl;
      return false;
    }
    Fout.resize(Eout.size());
    return modulate(Eorg.data(), Forg.data(), Eorg.size(), 
                    Eout.data(), Fout.data(), Eout.size());
  }
template <typename T>
  bool WMATH::solar_mod::modulate(const vector<T> & Eorg,       vector<T> & Forg) {
    vector<T> Fout;
    bool success_or_not = modulate(Eorg, Forg, Eorg, Fout);
    Forg = Fout;
    return success_or_not;
  }
// for pArray input
bool WMATH::solar_mod::modulate(const pArray & Eorg, const pArray & Forg,
                                const pArray & Eout,       pArray & Fout) {
  if (Eorg.GetLength() != Forg.GetLength()) {
    cout<<"[modulate_ERROR](pArray input): illegal length"<<endl;
    return false;
  }
  Fout.resize(Eout.GetLength());
  return modulate(Eorg.a, Forg.a, Eorg.GetLength(), 
                  Eout.a, Fout.a, Eout.GetLength());
}
bool WMATH::solar_mod::modulate(const pArray & Eorg,       pArray & Forg) {
  pArray Fout(Forg);
  bool success_or_not = modulate(Eorg, Forg, Eorg, Fout);
  Forg = Fout;
  return success_or_not;
}

// 显式实例化
#define SUPPORT(T)\
  template bool WMATH::solar_mod::modulate(const T *, const T *, int,\
                                           const T *,       T *, int);\
  template bool WMATH::solar_mod::modulate(const T *,       T *, int);\
  template bool WMATH::solar_mod::modulate(const vector<T> &, const vector<T> &,\
                                           const vector<T> &,       vector<T> &);\
  template bool WMATH::solar_mod::modulate(const vector<T> &,       vector<T> &);

SUPPORT(double)
SUPPORT(float)
