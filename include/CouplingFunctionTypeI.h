#ifndef CouplingFunctionTypeI
#define CouplingFunctionTypeI
#include "ModelParameters.h"
#include "clooptools.h"

namespace TypeI{
ComplexType dKappahWW(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32);

ComplexType dKappahZZ(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32);

ComplexType dKappahbb(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32);

ComplexType dKappahcc(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32);

ComplexType dKappahtautau(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32);

ComplexType dKappahtt(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32);

ComplexType dghgaga(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32);

ComplexType dghgg(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32);

} //end namespace TypeI
#endif
