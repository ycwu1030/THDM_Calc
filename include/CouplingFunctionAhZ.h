#ifndef CouplingFunctionAhZ
#define CouplingFunctionAhZ
#include "ModelParameters.h"
#include "clooptools.h"

namespace TypeI{
ComplexType CAhZ(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32);
}

namespace TypeII{
ComplexType CAhZ(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32);
}

namespace TypeLS{
ComplexType CAhZ(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32);
}

namespace TypeFL{
ComplexType CAhZ(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32);
}
#endif
