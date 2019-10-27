#ifndef CouplingFunctionSM
#define CouplingFunctionSM
#include "ModelParameters.h"
#include "clooptools.h"

namespace SM{
ComplexType dKappahWW(double m12, double m22, double m32);

ComplexType dKappahZZ(double m12, double m22, double m32);

ComplexType dKappahbb(double m12, double m22, double m32);

ComplexType dKappahcc(double m12, double m22, double m32);

ComplexType dKappahtautau(double m12, double m22, double m32);

ComplexType dKappahtt(double m12, double m22, double m32);

ComplexType dghgaga(double m12, double m22, double m32);

ComplexType dghgg(double m12, double m22, double m32);

} //end namespace SM
#endif
