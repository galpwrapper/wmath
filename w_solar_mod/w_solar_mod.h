#ifndef _W_SOLAR_MOD
#define _W_SOLAR_MOD
#include <vector>
#include "pArray.h"
#include "w_interp.h"
using std::vector;

namespace WMATH {
  
class solar_mod {
  public:
    solar_mod();
    solar_mod(int, int, double); // Z, A, phi
    solar_mod & set(int, int, double);

    // for array input
    template <typename T>
      bool modulate(const T *, const T *, int,
                    const T *,       T *, int);
    template <typename T>
      bool modulate(const T *,       T *, int);

    // for vector input
    template <typename T>
      bool modulate(const vector<T> &, const vector<T> &,
                    const vector<T> &,       vector<T> &);
    template <typename T>
      bool modulate(const vector<T> &,       vector<T> &);

    // for pArray input
    bool modulate(const pArray &, const pArray &,
                  const pArray &,       pArray &);
    bool modulate(const pArray &,       pArray &);
  private:
    int Z, A;
    double phi; // Gev
    double mass; // Gev
    const static double eMass; //Electron mass in Gev
    const static double pMass; //Proton mass in Gev
    void init(int, int, double);
};

}
#endif
