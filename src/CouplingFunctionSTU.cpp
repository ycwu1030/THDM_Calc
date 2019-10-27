#include "CouplingFunctionSTU.h"

namespace STU{
ComplexType SinTHDM(double MHH2, double MA02, double MHp2, double M2, double beta, double alp)
{
return (0.3183098861837907*(B0i(bb00, 0, MHp2, MHp2) - B0i(bb00, MZ2, MHp2, MHp2) - 
   MZ2*B0i(bb0, 0, Mh02, MZ2)*(pow(Cos(alp - beta), 2)) + 
   MZ2*B0i(bb0, 0, MHH2, MZ2)*(pow(Cos(alp - beta), 2)) + 
   MZ2*B0i(bb0, MZ2, Mh02, MZ2)*(pow(Cos(alp - beta), 2)) - 
   MZ2*B0i(bb0, MZ2, MHH2, MZ2)*(pow(Cos(alp - beta), 2)) - 
   B0i(bb00, 0, MA02, Mh02)*(pow(Cos(alp - beta), 2)) + 
   B0i(bb00, 0, Mh02, MZ2)*(pow(Cos(alp - beta), 2)) - 
   B0i(bb00, 0, MHH2, MZ2)*(pow(Cos(alp - beta), 2)) + 
   B0i(bb00, MZ2, MA02, Mh02)*(pow(Cos(alp - beta), 2)) - 
   B0i(bb00, MZ2, Mh02, MZ2)*(pow(Cos(alp - beta), 2)) + 
   B0i(bb00, MZ2, MHH2, MZ2)*(pow(Cos(alp - beta), 2)) - 
   B0i(bb00, 0, MA02, MHH2)*(pow(Sin(alp - beta), 2)) + 
   B0i(bb00, MZ2, MA02, MHH2)*(pow(Sin(alp - beta), 2))))/MZ2
;
}

ComplexType TinTHDM(double MHH2, double MA02, double MHp2, double M2, double beta, double alp)
{
return (0.039788735772973836*(-16*MW2*MZ2*B0i(bb00, 0, MHp2, MHp2) + 
   8.*B0i(bb00, 0, MHp2, MHp2)*(pow(MW2, 2)) + 
   2.*B0i(bb00, 0, MA02, MHp2)*(pow(MZ2, 2)) + 
   B0i(bb00, 0, Mh02, MHp2)*(pow(MZ2, 2)) - B0i(bb00, 0, Mh02, MW2)*
    (pow(MZ2, 2)) + B0i(bb00, 0, Mh02, MZ2)*(pow(MZ2, 2)) + 
   B0i(bb00, 0, MHH2, MHp2)*(pow(MZ2, 2)) + B0i(bb00, 0, MHH2, MW2)*
    (pow(MZ2, 2)) - B0i(bb00, 0, MHH2, MZ2)*(pow(MZ2, 2)) + 
   6.*B0i(bb00, 0, MHp2, MHp2)*(pow(MZ2, 2)) + 
   B0i(bb00, 0, Mh02, MHp2)*Cos(2*(alp - beta))*(pow(MZ2, 2)) - 
   B0i(bb00, 0, Mh02, MW2)*Cos(2*(alp - beta))*(pow(MZ2, 2)) + 
   B0i(bb00, 0, Mh02, MZ2)*Cos(2*(alp - beta))*(pow(MZ2, 2)) - 
   B0i(bb00, 0, MHH2, MHp2)*Cos(2*(alp - beta))*(pow(MZ2, 2)) + 
   B0i(bb00, 0, MHH2, MW2)*Cos(2*(alp - beta))*(pow(MZ2, 2)) - 
   B0i(bb00, 0, MHH2, MZ2)*Cos(2*(alp - beta))*(pow(MZ2, 2)) - 
   4.*A0i(aa0, MHp2)*(pow(MZ2, 2))*(pow(SW2, 2)) + 
   2.*MW2*B0i(bb0, 0, Mh02, MW2)*(pow(MZ2, 2))*
    (pow(Cos(alp - beta), 2)) - 2*MW2*B0i(bb0, 0, MHH2, MW2)*
    (pow(MZ2, 2))*(pow(Cos(alp - beta), 2)) - 
   2.*B0i(bb00, 0, MA02, Mh02)*(pow(MZ2, 2))*
    (pow(Cos(alp - beta), 2)) - 2.*B0i(bb0, 0, Mh02, MZ2)*
    (pow(MZ2, 3))*(pow(Cos(alp - beta), 2)) + 
   2.*B0i(bb0, 0, MHH2, MZ2)*(pow(MZ2, 3))*(pow(Cos(alp - beta), 2)) - 
   2.*B0i(bb00, 0, MA02, MHH2)*(pow(MZ2, 2))*
    (pow(Sin(alp - beta), 2))))/(MW2*SW2*(pow(MZ2, 2)))
;
}

ComplexType UinTHDM(double MHH2, double MA02, double MHp2, double M2, double beta, double alp)
{
return (0.15915494309189535*(-2*MZ2*B0i(bb00, 0, MA02, MHp2) - 
   2*(MW2 - 2*MZ2)*B0i(bb00, 0, MHp2, MHp2) + 
   2*MZ2*B0i(bb00, MW2, MA02, MHp2) + MZ2*B0i(bb00, MW2, Mh02, MHp2) - 
   MZ2*B0i(bb00, MW2, Mh02, MW2) + MZ2*B0i(bb00, MW2, MHH2, MHp2) + 
   MZ2*B0i(bb00, MW2, MHH2, MW2) - 4*MZ2*B0i(bb00, MW2, MHp2, MHp2) - 
   MW2*B0i(bb00, MZ2, MA02, Mh02) - MW2*B0i(bb00, MZ2, MA02, MHH2) + 
   MW2*B0i(bb00, MZ2, Mh02, MZ2) - MW2*B0i(bb00, MZ2, MHH2, MZ2) + 
   2*MW2*B0i(bb00, MZ2, MHp2, MHp2) + MZ2*B0i(bb00, MW2, Mh02, MHp2)*
    Cos(2*(alp - beta)) - MZ2*B0i(bb00, MW2, Mh02, MW2)*Cos(2*(alp - beta)) - 
   MZ2*B0i(bb00, MW2, MHH2, MHp2)*Cos(2*(alp - beta)) + 
   MZ2*B0i(bb00, MW2, MHH2, MW2)*Cos(2*(alp - beta)) - 
   MW2*B0i(bb00, MZ2, MA02, Mh02)*Cos(2*(alp - beta)) + 
   MW2*B0i(bb00, MZ2, MA02, MHH2)*Cos(2*(alp - beta)) + 
   MW2*B0i(bb00, MZ2, Mh02, MZ2)*Cos(2*(alp - beta)) - 
   MW2*B0i(bb00, MZ2, MHH2, MZ2)*Cos(2*(alp - beta)) - 
   2*MW2*MZ2*B0i(bb0, 0, Mh02, MW2)*(pow(Cos(alp - beta), 2)) + 
   2*MW2*MZ2*B0i(bb0, 0, Mh02, MZ2)*(pow(Cos(alp - beta), 2)) + 
   2*MW2*MZ2*B0i(bb0, 0, MHH2, MW2)*(pow(Cos(alp - beta), 2)) - 
   2*MW2*MZ2*B0i(bb0, 0, MHH2, MZ2)*(pow(Cos(alp - beta), 2)) + 
   2*MW2*MZ2*B0i(bb0, MW2, Mh02, MW2)*(pow(Cos(alp - beta), 2)) - 
   2*MW2*MZ2*B0i(bb0, MW2, MHH2, MW2)*(pow(Cos(alp - beta), 2)) - 
   2*MW2*MZ2*B0i(bb0, MZ2, Mh02, MZ2)*(pow(Cos(alp - beta), 2)) + 
   2*MW2*MZ2*B0i(bb0, MZ2, MHH2, MZ2)*(pow(Cos(alp - beta), 2)) + 
   2*MW2*B0i(bb00, 0, MA02, Mh02)*(pow(Cos(alp - beta), 2)) - 
   2*MZ2*B0i(bb00, 0, Mh02, MHp2)*(pow(Cos(alp - beta), 2)) + 
   2*MZ2*B0i(bb00, 0, Mh02, MW2)*(pow(Cos(alp - beta), 2)) - 
   2*MW2*B0i(bb00, 0, Mh02, MZ2)*(pow(Cos(alp - beta), 2)) - 
   2*MZ2*B0i(bb00, 0, MHH2, MW2)*(pow(Cos(alp - beta), 2)) + 
   2*MW2*B0i(bb00, 0, MHH2, MZ2)*(pow(Cos(alp - beta), 2)) + 
   2*MW2*B0i(bb00, 0, MA02, MHH2)*(pow(Sin(alp - beta), 2)) - 
   2*MZ2*B0i(bb00, 0, MHH2, MHp2)*(pow(Sin(alp - beta), 2))))/(MW2*MZ2)
;
}

} //end namespace STU
