#ifndef CouplingFunctionSTU
#define CouplingFunctionSTU
#include "ModelParameters.h"
#include "clooptools.h"

namespace STU{
ComplexType SinTHDM(double MHH2, double MA02, double MHp2, double M2, double beta, double alp);

ComplexType TinTHDM(double MHH2, double MA02, double MHp2, double M2, double beta, double alp);

ComplexType UinTHDM(double MHH2, double MA02, double MHp2, double M2, double beta, double alp);

} //end namespace STU
#endif
