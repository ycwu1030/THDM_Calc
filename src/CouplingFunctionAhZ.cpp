#include "CouplingFunctionAhZ.h"

namespace TypeI{
ComplexType CAhZ(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32)
{
return 0. - (0.5*Cos(alp - beta))/(CW*SW) + 
 (0.009947183943243459*Alfa*MZ2*B0i(bb0, m12, Mh02, MZ2)*Cos(alp - beta))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*B0i(bb0, m12, MHp2, MW2)*
   Cos(alp - beta))/(CW*SW) + (0.009947183943243459*Alfa*MZ2*
   B0i(bb0, m22, MA02, MZ2)*Cos(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*B0i(bb0, m22, MHp2, MW2)*Cos(alp - beta))/
  (CW*SW) - (0.07957747154594767*Alfa*CW*B0i(bb0, m32, MW2, MW2)*
   Cos(alp - beta))/(SW*SW2) - (0.009947183943243459*Alfa*MZ2*
   B0i(bb1, m12, Mh02, MZ2)*Cos(alp - beta))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*B0i(bb1, m12, MHp2, MW2)*Cos(alp - beta))/
  (CW*SW) - (0.009947183943243459*Alfa*MZ2*B0i(bb1, m22, MA02, MZ2)*
   Cos(alp - beta))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*B0i(bb1, m22, MHp2, MW2)*Cos(alp - beta))/
  (CW*SW) + (0.07957747154594767*Alfa*CW*C0i(cc00, m22, m32, m12, MHp2, 
    MW2, MW2)*Cos(alp - beta))/(SW*SW2) + 
 (0.07957747154594767*Alfa*(2*MW2 - MZ2)*C0i(cc00, m32, m12, m22, MHp2, 
    MHp2, MW2)*Cos(alp - beta))/(CW*MZ2*SW*SW2) - 
 (0.019894367886486918*Alfa*(MW2*(2*(MA02 - MHp2)*(-Mh02 + MHp2) - 
      (MA02 + Mh02 - 2*MHp2 - m12 + m22 + 2*m32)*MW2) + 
    ((MA02 - MHp2)*(Mh02 - MHp2) + (MA02 + Mh02 - 2*MHp2)*MW2)*MZ2)*
   C0i(cc1, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/
  (MW*MW2*MZ*SW*SW2) + (0.019894367886486918*Alfa*CW*(m12 + 3*m22 - m32)*
   C0i(cc11, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/(SW*SW2) + 
 (0.039788735772973836*Alfa*CW*(2*(m12 + m22) - m32)*
   C0i(cc12, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/(SW*SW2) - 
 (0.039788735772973836*Alfa*(m12 - m22)*(2*MW2 - MZ2)*
   C0i(cc12, m32, m12, m22, MHp2, MHp2, MW2)*Cos(alp - beta))/
  (CW*MZ2*SW*SW2) - (0.019894367886486918*Alfa*
   (MW2*(2*(MA02 - MHp2)*(-Mh02 + MHp2) - (MA02 + Mh02 - 2*MHp2 + m12 - 
        m22 + 2*m32)*MW2) + ((MA02 - MHp2)*(Mh02 - MHp2) + 
      (MA02 + Mh02 - 2*MHp2)*MW2)*MZ2)*C0i(cc2, m22, m32, m12, MHp2, MW2, 
    MW2)*Cos(alp - beta))/(MW*MW2*MZ*SW*SW2) - 
 (0.019894367886486918*Alfa*((MA02 - MHp2)*(-Mh02 + MHp2) + 
    MW2*(m12 - m22 - m32 + MW2))*(2*MW2 - MZ2)*
   C0i(cc2, m32, m12, m22, MHp2, MHp2, MW2)*Cos(alp - beta))/
  (CW*MW2*MZ2*SW*SW2) + (0.019894367886486918*Alfa*CW*(3*m12 + m22 - m32)*
   C0i(cc22, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/(SW*SW2) + 
 (0.019894367886486918*Alfa*(m12 + 3*m22 - m32)*(2*MW2 - MZ2)*
   C0i(cc22, m32, m12, m22, MHp2, MHp2, MW2)*Cos(alp - beta))/
  (CW*MZ2*SW*SW2) + (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MB2, MB2, MB2)*
   Cos(alp)*Cot(beta)*Csc(beta)*(pow(MB2, 2)))/
  (-8*MW*MW2*MZ2*Pi*SW + 8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MC2, MC2, MC2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MC2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MD2, MD2, MD2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MD2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, ME2, ME2, ME2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(ME2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, ML2, ML2, ML2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(ML2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MM2, MM2, MM2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MM2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MS2, MS2, MS2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MS2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MT2, MT2, MT2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MT2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MU2, MU2, MU2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MU2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MB2*MZ*MZ2*B0i(bb0, m32, MB2, MB2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MC2*MZ*MZ2*B0i(bb0, m32, MC2, MC2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MD2*MZ*MZ2*B0i(bb0, m32, MD2, MD2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ME2*MZ*MZ2*B0i(bb0, m32, ME2, ME2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ML2*MZ*MZ2*B0i(bb0, m32, ML2, ML2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MM2*MZ*MZ2*B0i(bb0, m32, MM2, MM2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MS2*MZ*MZ2*B0i(bb0, m32, MS2, MS2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MT2*MZ*MZ2*B0i(bb0, m32, MT2, MT2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MU2*MZ*MZ2*B0i(bb0, m32, MU2, MU2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MB2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MB2, MB2, MB2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MC2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MC2, MC2, MC2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MD2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MD2, MD2, MD2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ME2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, ME2, ME2, ME2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ML2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, ML2, ML2, ML2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MM2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MM2, MM2, MM2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MS2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MS2, MS2, MS2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MT2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MT2, MT2, MT2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MU2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MU2, MU2, MU2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MB2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MB2, MB2, MB2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MC2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MC2, MC2, MC2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MD2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MD2, MD2, MD2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ME2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, ME2, ME2, ME2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ML2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, ML2, ML2, ML2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MM2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MM2, MM2, MM2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MS2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MS2, MS2, MS2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MT2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MT2, MT2, MT2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*m12*MU2*MZ*MZ2*C0i(cc2, m22, m32, m12, MU2, MU2, MU2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (0.019894367886486918*Alfa*MZ2*C0i(cc0, m22, m32, m12, MHp2, MW2, MW2)*
   Cos(alp - beta)*(Mh02 - MHp2 - ((4*MHp2 + m32)*MW2 + 
      (MA02 - MHp2)*(MW2*(2*MHp2 + MW2) - (MHp2 + MW2)*MZ2 + 
        Mh02*(-2*MW2 + MZ2))*(pow(MW2, -1)))*(pow(-MW2 + MZ2, -1))))/
  (MW*MZ*SW) - (0.039788735772973836*Alfa*MZ2*C0i(cc00, m32, m12, m22, 
    MA02, Mh02, MZ2)*(pow(Cos(alp - beta), 3)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(m12 - m22)*MZ2*C0i(cc12, m32, m12, m22, 
    MA02, Mh02, MZ2)*(pow(Cos(alp - beta), 3)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*(m12 + 3*m22 - m32)*MZ2*
   C0i(cc22, m32, m12, m22, MA02, Mh02, MZ2)*
   (pow(Cos(alp - beta), 3)))/(CW*MW2*SW*SW2) + 
 (0.03125*Alfa2*v2*C0i(cc2, m32, m12, m22, MA02, Mh02, MZ2)*
   (MZ2*(m12 - m22 - m32 + MZ2) + (pow(MA02 - Mh02, 2)))*
   (pow(Cos(alp - beta), 3))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0234375*Alfa2*v2*C0i(cc1, m12, m22, m32, MA02, Mh02, Mh02)*
   Cos(alp - beta)*(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh02)*Cos(alp + beta))*((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta))*(pow(Csc(2*beta), 2))*
   (pow(Csc(w), 2)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*v2*C0i(cc0, m22, m32, m12, MA02, MA02, Mh02)*
   Cos(alp - beta)*(pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
      (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta), 2))*
   (pow(Csc(2*beta), 2))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*v2*C0i(cc1, m22, m32, m12, MA02, MA02, Mh02)*
   Cos(alp - beta)*(pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
      (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta), 2))*
   (pow(Csc(2*beta), 2))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*v2*C0i(cc2, m22, m32, m12, MA02, MA02, Mh02)*
   Cos(alp - beta)*(pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
      (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta), 2))*
   (pow(Csc(2*beta), 2))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0049735919716217296*Alfa*Mh02*(MA02 - Mh02 + MZ2)*B0i(bb0, 0, Mh02, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MA02*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*MHH2*(MA02 - MHH2 + MZ2)*B0i(bb0, 0, MHH2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MA02*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*(Mh02*(-Mh02 + MZ2) + MA02*(Mh02 + MZ2))*
   B0i(bb0, MA02, Mh02, MZ2)*Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/
  (CW*MA02*MW2*SW*SW2) + (0.0049735919716217296*Alfa*
   (MHH2*(-MHH2 + MZ2) + MA02*(MHH2 + MZ2))*B0i(bb0, MA02, MHH2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MA02*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MZ2*B0i(bb1, MA02, Mh02, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MZ2*B0i(bb1, MA02, MHH2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*MZ2*(-MA02 + Mh02 + MZ2)*
   C0i(cc0, m12, m22, m32, Mh02, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*MZ2*C0i(cc00, m12, m22, m32, Mh02, MZ2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.039788735772973836*Alfa*MZ2*C0i(cc00, m12, m22, m32, MHH2, MZ2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*MZ2*C0i(cc00, m32, m12, m22, MA02, MHH2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*(3*m12 + m22 - m32)*MZ2*
   C0i(cc11, m12, m22, m32, Mh02, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*(3*m12 + m22 - m32)*MZ2*
   C0i(cc11, m12, m22, m32, MHH2, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*(m12 - m22)*MZ2*C0i(cc12, m12, m22, m32, 
    Mh02, MZ2, MZ2)*Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*(m12 - m22)*MZ2*
   C0i(cc12, m12, m22, m32, MHH2, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(m12 - m22)*MZ2*C0i(cc12, m32, m12, m22, 
    MA02, MHH2, MZ2)*Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2) - (0.009947183943243459*Alfa*(m12 + 3*m22 - m32)*MZ2*
   C0i(cc22, m32, m12, m22, MA02, MHH2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*C0i(cc1, m12, m22, m32, Mh02, MZ2, MZ2)*
   Cos(alp - beta)*((MA02 - Mh02)*Mh02 + (MA02 - Mh02 + m12 - m22 + m32)*
     MZ2 + (pow(MZ2, 2)))*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2) - (0.0625*Alfa2*MZ2*(-MA02 + MHH2 + MZ2)*v2*
   C0i(cc0, m12, m22, m32, MHH2, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Csc(w), 2))*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*(-(Mh02*MHH2) + MA02*(Mh02 + MZ2) + 
    MZ2*(-MHH2 + m12 - m22 + m32 + MZ2))*v2*C0i(cc1, m12, m22, m32, 
    MHH2, MZ2, MZ2)*Cos(alp - beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*v2*C0i(cc2, m32, m12, m22, MA02, MHH2, MZ2)*
   Cos(alp - beta)*(Mh02*MHH2 - MA02*(Mh02 + MHH2) + 
    MZ2*(m12 - m22 - m32 + MZ2) + (pow(MA02, 2)))*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-46*MW2*MZ2 + 20*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(A0i(aa0, MB2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-70*MW2*MZ2 + 56*(pow(MW2, 2)) + 23*(pow(MZ2, 2)))*
   Re(A0i(aa0, MC2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-46*MW2*MZ2 + 20*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(A0i(aa0, MD2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-22*MW2*MZ2 + 12*(pow(MW2, 2)) + 9*(pow(MZ2, 2)))*
   Re(A0i(aa0, ME2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*Re(A0i(aa0, MHp2)))/
  (CW*MW2*MZ2*SW*SW2) - (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-22*MW2*MZ2 + 12*(pow(MW2, 2)) + 9*(pow(MZ2, 2)))*
   Re(A0i(aa0, ML2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-22*MW2*MZ2 + 12*(pow(MW2, 2)) + 9*(pow(MZ2, 2)))*
   Re(A0i(aa0, MM2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-46*MW2*MZ2 + 20*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(A0i(aa0, MS2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-70*MW2*MZ2 + 56*(pow(MW2, 2)) + 23*(pow(MZ2, 2)))*
   Re(A0i(aa0, MT2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-70*MW2*MZ2 + 56*(pow(MW2, 2)) + 23*(pow(MZ2, 2)))*
   Re(A0i(aa0, MU2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(4*MW2 - 3*MZ2)*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(A0i(aa0, MW2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*Re(A0i(aa0, MZ2)))/
  (CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) + 
 (0.17904931097838225*Alfa*Cos(alp - beta)*Re(B0i(bb0, 0, MW2, MW2)))/
  (CW*SW) - (0.0049735919716217296*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb0, MA02, Mh02, MZ2)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MA02, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.009947183943243459*Alfa*Cos(alp - beta)*
   Re(B0i(bb0, MA02, MHp2, MW2)))/(CW*SW*SW2) + 
 (0.0024867959858108648*Alfa*Cos(alp - beta)*(Mh02*MHH2 - 2*Mh02*MZ2 + 
    MHH2*MZ2 - MA02*(Mh02 + MHH2 + MZ2) + Cos(2*(alp - beta))*
     (MHH2*(-Mh02 + MZ2) + MA02*(Mh02 + MHH2 + MZ2) - (pow(MA02, 2))) + 
    (pow(MA02, 2)))*Re(B0i(bb0, Mh02, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.0049735919716217296*Alfa*Cos(alp - beta)*
   ((Mh02 - MHp2)*(MHH2 - MHp2) + (-2*Mh02 + MHH2 - MHp2)*MW2 + 
    ((Mh02 - MHp2)*(-MHH2 + MHp2) + (MHH2 + MHp2)*MW2)*Cos(2*(alp - beta)))*
   Re(B0i(bb0, Mh02, MHp2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*(Mh02*MHH2 + 2*MW2*(-MHH2 + 6*MW2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, MW2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.0024867959858108648*Alfa*(Mh02*MHH2 + 2*MZ2*(-MHH2 + 6*MZ2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*Cos(alp - beta)*(MHH2*(-Mh02 + MZ2) + 
    MA02*(Mh02 + MHH2 + MZ2) - (pow(MA02, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MHH2, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*((Mh02 - MHp2)*(MHH2 - MHp2) - (MHH2 + MHp2)*MW2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MHp2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*(Mh02*MHH2 + 2*MW2*(-MHH2 + 6*MW2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MW2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0024867959858108648*Alfa*(Mh02*MHH2 + 2*MZ2*(-MHH2 + 6*MZ2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.07957747154594767*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, 0, MW2)))/(CW*MZ2*SW*SW2) + 
 (0.05968310365946075*Alfa*(MT2 - MW2)*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, MB2, MT2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.05968310365946075*Alfa*MC2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, MC2, MS2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.05968310365946075*Alfa*(MU2 - MW2)*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, MD2, MU2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MW2, Mh02, MW2)))/
  (CW*MZ2*SW*(pow(SW2, 2))) - (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   (pow(Cos(alp - beta), 3))*Re(B0i(bb0, MW2, MHH2, MW2)))/
  (CW*MZ2*SW*(pow(SW2, 2))) + (0.019894367886486918*Alfa*Cos(alp - beta)*
   (6*MZ2*(pow(MW2, 2)) + 4*(pow(MW2, 3)) - 6*MW2*(pow(MZ2, 2)) + 
    (pow(MZ2, 3)))*Re(B0i(bb0, MW2, MW2, MZ2)))/
  (CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MB2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MB2, MB2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MC2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MC2, MC2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MD2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MD2, MD2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*ME2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, ME2, ME2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MZ2, Mh02, MZ2)))/
  (CW*MW2*SW*(pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   (pow(Cos(alp - beta), 3))*Re(B0i(bb0, MZ2, MHH2, MZ2)))/
  (CW*MW2*SW*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*ML2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, ML2, ML2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*MM2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MM2, MM2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MS2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MS2, MS2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MT2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MT2, MT2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MU2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MU2, MU2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-8*MW2*MZ2 + 13*(pow(MW2, 2)) + 2*(pow(MZ2, 2)))*
   Re(B0i(bb0, MZ2, MW2, MW2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.026525823848649224*Alfa*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MB2, MB2)))/(CW*MW2*MZ2*SW) - 
 (0.05305164769729845*Alfa*(8*MW2 - 5*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MC2, MC2)))/(CW*MW2*MZ2*SW) - 
 (0.026525823848649224*Alfa*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MD2, MD2)))/(CW*MW2*MZ2*SW) - 
 (0.07957747154594767*Alfa*(4*MW2 - 3*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, ME2, ME2)))/(CW*MW2*MZ2*SW) + 
 (0.07957747154594767*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MHp2, MHp2)))/(CW*MW2*MZ2*SW) - 
 (0.07957747154594767*Alfa*(4*MW2 - 3*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, ML2, ML2)))/(CW*MW2*MZ2*SW) - 
 (0.07957747154594767*Alfa*(4*MW2 - 3*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MM2, MM2)))/(CW*MW2*MZ2*SW) - 
 (0.026525823848649224*Alfa*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MS2, MS2)))/(CW*MW2*MZ2*SW) - 
 (0.05305164769729845*Alfa*(8*MW2 - 5*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MT2, MT2)))/(CW*MW2*MZ2*SW) - 
 (0.05305164769729845*Alfa*(8*MW2 - 5*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MU2, MU2)))/(CW*MW2*MZ2*SW) + 
 (0.07957747154594767*Alfa*(6*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MW2, MW2)))/(CW*MW2*MZ2*SW) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, ME2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, ML2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, MM2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.15915494309189535*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, MW2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MA02, MHp2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.1193662073189215*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MB2, MT2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.1193662073189215*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MC2, MS2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.1193662073189215*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MD2, MU2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MW2, Mh02, MHp2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MW2, Mh02, MW2)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MW2, MHH2, MHp2)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MW2, MHH2, MW2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(8*MW2 + MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MW2, MZ2)))/(CW*MW2*SW*(pow(MZ2, 2))*
   (pow(SW2, 2))) + (0.05968310365946075*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*Re(B0i(bb00, MZ2, 0, 0)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MZ2, MA02, Mh02)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MZ2, MA02, MHH2)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.006631455962162306*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MB2, MB2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.006631455962162306*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 
    17*(pow(MZ2, 2)))*Re(B0i(bb00, MZ2, MC2, MC2)))/
  (CW*MW2*SW*(pow(MZ2, 3))*(pow(SW2, 2))) + 
 (0.006631455962162306*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MD2, MD2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, ME2, ME2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) - (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb00, MZ2, Mh02, MZ2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MZ2, MHH2, MZ2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*(pow(2*MW2 - MZ2, 3))*
   Re(B0i(bb00, MZ2, MHp2, MHp2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, ML2, ML2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MM2, MM2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.006631455962162306*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MS2, MS2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.006631455962162306*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 
    17*(pow(MZ2, 2)))*Re(B0i(bb00, MZ2, MT2, MT2)))/
  (CW*MW2*SW*(pow(MZ2, 3))*(pow(SW2, 2))) + 
 (0.006631455962162306*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MU2, MU2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) - (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-4*MW2*MZ2 + 12*(pow(MW2, 2)) + (pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MW2, MW2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.026525823848649224*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MB2, MB2)))/(CW*SW) + 
 (0.1061032953945969*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MC2, MC2)))/
  (CW*SW) + (0.026525823848649224*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MD2, MD2)))/(CW*SW) + 
 (0.07957747154594767*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, ME2, ME2)))/
  (CW*SW) + (0.07957747154594767*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, ML2, ML2)))/(CW*SW) + 
 (0.07957747154594767*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MM2, MM2)))/
  (CW*SW) + (0.026525823848649224*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MS2, MS2)))/(CW*SW) + 
 (0.1061032953945969*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MT2, MT2)))/
  (CW*SW) + (0.1061032953945969*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MU2, MU2)))/(CW*SW) + 
 (0.039788735772973836*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MW2, MW2)))/
  (CW*SW) - (0.05968310365946075*Alfa*MB2*Cos(alp - beta)*
   (pow(Cot(beta), 2))*Re(B0i(bb1, MA02, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MC2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MD2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*ME2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, ME2, ME2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb1, MA02, Mh02, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MZ2*Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb1, MA02, MHH2, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*Re(B0i(bb1, MA02, MHp2, MW2)))/
  (CW*SW*SW2) - (0.019894367886486918*Alfa*ML2*Cos(alp - beta)*
   (pow(Cot(beta), 2))*Re(B0i(bb1, MA02, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MM2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MS2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MT2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MU2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MU2, MU2)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*MZ2*Cos(alp - beta)*
   (-2*Mh02 + MHH2 + MHH2*Cos(2*(alp - beta)))*Re(B0i(bb1, Mh02, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-2*Mh02 + MHH2 + MHH2*Cos(2*(alp - beta)))*Re(B0i(bb1, Mh02, MHp2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) - (0.019894367886486918*Alfa*MHH2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, Mh02, MW2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) - (0.009947183943243459*Alfa*MHH2*MZ2*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb1, Mh02, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MHH2*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, MHH2, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MHH2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, MHH2, MHp2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) + (0.019894367886486918*Alfa*MHH2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, MHH2, MW2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) + (0.009947183943243459*Alfa*MHH2*MZ2*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb1, MHH2, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, ME2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, ML2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, MM2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, MW2)))/(CW*MZ2*SW*SW2) - 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MB2, MT2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MC2, MS2)))/(CW*MZ2*SW*(pow(SW2, 2))) - 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MD2, MU2)))/(CW*MZ2*SW*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*MW2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MW2, MZ2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*Cos(alp - beta)*Re(B0i(bb1, MZ2, 0, 0)))/
  (CW*SW*(pow(SW2, 2))) - (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MC2, MC2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, ME2, ME2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, ML2, ML2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MM2, MM2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MT2, MT2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MU2, MU2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*Cos(alp - beta)*(pow(MW2, 2))*
   Re(B0i(bb1, MZ2, MW2, MW2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) + 
 (0.0012433979929054324*Alfa*Cos(alp - beta)*
   (pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*Mh02)*
       Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, MA02, MA02, Mh02)))/(CW*MW2*SW*SW2) + 
 (0.0012433979929054324*Alfa*Cos(alp - beta)*(pow(Csc(2*beta), 2))*
   (pow((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
       Sin(alp + beta), 2))*Re(B0i(dbb0, MA02, MA02, MHH2)))/
  (CW*MW2*SW*SW2) - (0.0049735919716217296*Alfa*
   ((MA02 + Mh02)*MZ2 - (pow(MA02 - Mh02, 2)))*
   (pow(Cos(alp - beta), 3))*Re(B0i(dbb0, MA02, Mh02, MZ2)))/
  (CW*MW2*SW*SW2) - (0.0049735919716217296*Alfa*Cos(alp - beta)*
   ((MA02 + MHH2)*MZ2 - (pow(MA02 - MHH2, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, MA02, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) + (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-((MA02 + MHp2)*MW2) + (pow(MA02 - MHp2, 2)))*
   Re(B0i(dbb0, MA02, MHp2, MW2)))/(CW*MW2*SW*SW2) + 
 (0.0006216989964527162*Alfa*Cos(alp - beta)*
   (pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*Mh02)*
       Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh02, MA02, MA02)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*((MA02 + Mh02)*MZ2 - (pow(MA02 - Mh02, 2)))*
   (pow(Cos(alp - beta), 3))*Re(B0i(dbb0, Mh02, MA02, MZ2)))/
  (CW*MW2*SW*SW2) + (0.005595290968074445*Alfa*Cos(alp - beta)*
   (pow(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
      (2*M2 - 3*Mh02)*Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh02, Mh02, Mh02)))/(CW*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*(pow(Cos(alp - beta), 3))*
   (pow(M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*Sin(2*alp), 2))*
   Re(B0i(dbb0, Mh02, Mh02, MHH2)))/(CW*MW2*SW*SW2) + 
 (0.0024867959858108648*Alfa*Cos(alp - beta)*(pow(Csc(2*beta), 2))*
   (pow(Sin(alp - beta), 2))*
   (pow((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta), 2))*
   Re(B0i(dbb0, Mh02, MHH2, MHH2)))/(CW*MW2*SW*SW2) + 
 (0.0012433979929054324*Alfa*Cos(alp - beta)*
   (pow((Mh02 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh02 + 2*MHp2)*
       Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh02, MHp2, MHp2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*(-((Mh02 + MHp2)*MW2) + 
    (pow(Mh02 - MHp2, 2)))*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb0, Mh02, MHp2, MW2)))/(CW*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*Cos(alp - beta)*
   (-2*Mh02*MW2 + (pow(Mh02, 2)) + 12*(pow(MW2, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, Mh02, MW2, MW2)))/
  (CW*MW2*SW*SW2) + (0.0024867959858108648*Alfa*Cos(alp - beta)*
   (-2*Mh02*MZ2 + (pow(Mh02, 2)) + 12*(pow(MZ2, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, Mh02, MZ2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.029841551829730376*Alfa*MB2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MC2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MD2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*ME2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, ME2, ME2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*(pow(MZ2, 2))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, MZ2, Mh02, MZ2)))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*(pow(MZ2, 2))*
   (pow(Cos(alp - beta), 3))*Re(B0i(dbb0, MZ2, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.009947183943243459*Alfa*ML2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MM2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MS2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MT2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MU2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MU2, MU2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*(9*MW2 - 2*MZ2)*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MW2, MW2)))/(CW*SW*SW2) - 
 (0.05305164769729845*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MB2, MB2)))/
  (CW*SW) - (0.2122065907891938*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, MC2, MC2)))/(CW*SW) - 
 (0.05305164769729845*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MD2, MD2)))/
  (CW*SW) - (0.15915494309189535*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, ME2, ME2)))/(CW*SW) + 
 (0.07957747154594767*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MHp2, MHp2)))/
  (CW*SW) - (0.15915494309189535*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, ML2, ML2)))/(CW*SW) - 
 (0.15915494309189535*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MM2, MM2)))/
  (CW*SW) - (0.05305164769729845*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, MS2, MS2)))/(CW*SW) - 
 (0.2122065907891938*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MT2, MT2)))/
  (CW*SW) - (0.2122065907891938*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, MU2, MU2)))/(CW*SW) + 
 (0.238732414637843*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MW2, MW2)))/
  (CW*SW) + (0.05968310365946075*Alfa*MZ2*Cos(alp - beta)*
   Re(B0i(dbb00, MZ2, 0, 0)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb00, MZ2, MA02, Mh02)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(dbb00, MZ2, MA02, MHH2)))/(CW*MW2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MB2, MB2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MC2, MC2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MD2, MD2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, ME2, ME2)))/(CW*MW2*MZ2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(dbb00, MZ2, Mh02, MZ2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb00, MZ2, MHH2, MZ2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*(pow(-2*MW2 + MZ2, 2))*
   Re(B0i(dbb00, MZ2, MHp2, MHp2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, ML2, ML2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MM2, MM2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MS2, MS2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MT2, MT2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MU2, MU2)))/(CW*MW2*MZ2*SW*SW2) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 12*(pow(MW2, 2)) + (pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MW2, MW2)))/(CW*MW2*MZ2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MB2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MC2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MD2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MA02*ME2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, ME2, ME2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MA02*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb1, MA02, Mh02, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MA02*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb1, MA02, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*MA02*Cos(alp - beta)*
   Re(B0i(dbb1, MA02, MHp2, MW2)))/(CW*SW*SW2) - 
 (0.019894367886486918*Alfa*MA02*ML2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MA02*MM2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MS2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MT2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MU2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MU2, MU2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*Mh02*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb1, Mh02, MA02, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*Mh02*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb1, Mh02, MHp2, MW2)))/(CW*SW*SW2) + 
 (0.019894367886486918*Alfa*Mh02*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb1, Mh02, MW2, MW2)))/(CW*SW*SW2) + 
 (0.009947183943243459*Alfa*Mh02*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb1, Mh02, MZ2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.029841551829730376*Alfa*Cos(alp - beta)*
   (pow(MZ2, 2))*Re(B0i(dbb1, MZ2, 0, 0)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, ME2, ME2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MU2, MU2)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*MW2*Cos(alp - beta)*Re(B0i(dbb1, MZ2, MW2, MW2)))/
  (CW*SW*SW2) + (0.05968310365946075*Alfa*MB2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MB2, MB2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MC2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MC2, MC2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MD2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MD2, MD2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*ME2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, ME2, ME2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*ML2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, ML2, ML2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*MM2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MM2, MM2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MS2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MS2, MS2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MT2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MT2, MT2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MU2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MU2, MU2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*(pow(MB2, 2))*
   Re(B0i(bb0, Mh02, MB2, MB2))*Sin(alp)*(Cos(alp) - Cot(beta)*Sin(alp)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*
   (pow(MC2, 2))*Re(B0i(bb0, Mh02, MC2, MC2))*Sin(alp)*
   (Cos(alp) - Cot(beta)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*(pow(MD2, 2))*
   Re(B0i(bb0, Mh02, MD2, MD2))*Sin(alp)*(Cos(alp) - Cot(beta)*Sin(alp)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*
   (pow(MS2, 2))*Re(B0i(bb0, Mh02, MS2, MS2))*Sin(alp)*
   (Cos(alp) - Cot(beta)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*(pow(MT2, 2))*
   Re(B0i(bb0, Mh02, MT2, MT2))*Sin(alp)*(Cos(alp) - Cot(beta)*Sin(alp)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*
   (pow(MU2, 2))*Re(B0i(bb0, Mh02, MU2, MU2))*Sin(alp)*
   (Cos(alp) - Cot(beta)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (1.*Re(B0i(bb0, Mh02, ME2, ME2))*
   (-(Alfa*Cos(alp)*Cot(beta)*Csc(beta)*(pow(ME2, 2))*
      (pow(Sin(alp), 2))) + Alfa*Csc(beta)*(pow(ME2, 2))*
     (pow(Cos(alp), 2))*Sin(alp)))/(8*CW*Mh02*MW2*Pi*SW*SW2 - 
   8*CW*MHH2*MW2*Pi*SW*SW2) + (1.*Re(B0i(bb0, Mh02, ML2, ML2))*
   (-(Alfa*Cos(alp)*Cot(beta)*Csc(beta)*(pow(ML2, 2))*
      (pow(Sin(alp), 2))) + Alfa*Csc(beta)*(pow(ML2, 2))*
     (pow(Cos(alp), 2))*Sin(alp)))/(8*CW*Mh02*MW2*Pi*SW*SW2 - 
   8*CW*MHH2*MW2*Pi*SW*SW2) + (1.*Re(B0i(bb0, Mh02, MM2, MM2))*
   (-(Alfa*Cos(alp)*Cot(beta)*Csc(beta)*(pow(MM2, 2))*
      (pow(Sin(alp), 2))) + Alfa*Csc(beta)*(pow(MM2, 2))*
     (pow(Cos(alp), 2))*Sin(alp)))/(8*CW*Mh02*MW2*Pi*SW*SW2 - 
   8*CW*MHH2*MW2*Pi*SW*SW2) + (0.03125*Alfa2*(-MA02 + Mh02 + MZ2)*v2*
   C0i(cc0, m22, m32, m12, Mh02, MHH2, MZ2)*(pow(Cos(alp - beta), 3))*
   (pow(Csc(w), 2))*(M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*
     Sin(2*alp)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.03125*Alfa2*(MA02 - Mh02 + MZ2)*v2*C0i(cc1, m22, m32, m12, Mh02, MHH2, 
    MZ2)*(pow(Cos(alp - beta), 3))*(pow(Csc(w), 2))*
   (M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*Sin(2*alp)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.03125*Alfa2*(MA02 - Mh02 + MZ2)*v2*C0i(cc2, m22, m32, m12, Mh02, MHH2, 
    MZ2)*(pow(Cos(alp - beta), 3))*(pow(Csc(w), 2))*
   (M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*Sin(2*alp)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.05968310365946075*Alfa*MB2*B0i(bb1, MA02, MB2, MB2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MC2*B0i(bb1, MA02, MC2, MC2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MD2*B0i(bb1, MA02, MD2, MD2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*ME2*B0i(bb1, MA02, ME2, ME2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*ML2*B0i(bb1, MA02, ML2, ML2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*MM2*B0i(bb1, MA02, MM2, MM2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MS2*B0i(bb1, MA02, MS2, MS2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MT2*B0i(bb1, MA02, MT2, MT2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MU2*B0i(bb1, MA02, MU2, MU2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*(pow(MB2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MB2, MB2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.05968310365946075*Alfa*(pow(MC2, 2))*
   (pow(Csc(beta), 2))*Re(B0i(bb0, MHH2, MC2, MC2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*(pow(MD2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MD2, MD2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(pow(ME2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, ME2, ME2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(pow(ML2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, ML2, ML2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(pow(MM2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MM2, MM2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.05968310365946075*Alfa*(pow(MS2, 2))*
   (pow(Csc(beta), 2))*Re(B0i(bb0, MHH2, MS2, MS2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*(pow(MT2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MT2, MT2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.05968310365946075*Alfa*(pow(MU2, 2))*
   (pow(Csc(beta), 2))*Re(B0i(bb0, MHH2, MU2, MU2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MB2*MHH2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MB2, MB2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MC2*MHH2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MC2, MC2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MD2*MHH2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MD2, MD2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*ME2*MHH2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, ME2, ME2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MHH2*ML2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, ML2, ML2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MHH2*MM2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MM2, MM2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MHH2*MS2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MS2, MS2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MHH2*MT2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MT2, MT2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MHH2*MU2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MU2, MU2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.0012433979929054324*Alfa*(MA02 - Mh02)*
   B0i(bb0, 0, MA02, Mh02)*((2*MA02 - Mh02)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh02)*Cos(alp + beta))*Csc(2*beta)*
   Sin(2*(alp - beta)))/(CW*MA02*MW2*SW*SW2) - 
 (0.0012433979929054324*Alfa*(MA02 - Mh02)*B0i(bb0, MA02, MA02, Mh02)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(2*beta)*Sin(2*(alp - beta)))/(CW*MA02*MW2*SW*SW2) + 
 (0.003730193978716297*Alfa*(-MA02 + Mh02 + MZ2)*C0i(cc0, m22, m32, m12, 
    Mh02, Mh02, MZ2)*(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh02)*Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) - (0.003730193978716297*Alfa*(MA02 - Mh02 + MZ2)*
   C0i(cc1, m22, m32, m12, Mh02, Mh02, MZ2)*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh02)*Cos(3*alp - beta) + (2*M2 - 3*Mh02)*Cos(alp + beta))*
   Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/(CW*MW2*SW*SW2) - 
 (0.003730193978716297*Alfa*(MA02 - Mh02 + MZ2)*C0i(cc2, m22, m32, m12, 
    Mh02, Mh02, MZ2)*(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh02)*Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) + (0.0012433979929054324*Alfa*(MA02 - Mh02 + MZ2)*
   C0i(cc0, m12, m32, m22, MA02, Mh02, MZ2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) - (0.0012433979929054324*Alfa*(-MA02 + Mh02 + MZ2)*
   C0i(cc1, m12, m32, m22, MA02, Mh02, MZ2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) - (0.0012433979929054324*Alfa*(-MA02 + Mh02 + MZ2)*
   C0i(cc2, m12, m32, m22, MA02, Mh02, MZ2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) - (0.03125*Alfa2*(-MA02 + MHH2 + MZ2)*v2*
   C0i(cc0, m22, m32, m12, MHH2, MHH2, MZ2)*Cos(alp - beta)*Csc(2*beta)*
   (pow(Csc(w), 2))*(pow(Sin(alp - beta), 2))*
   ((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*(MA02 - MHH2 + MZ2)*v2*C0i(cc1, m22, m32, m12, MHH2, MHH2, 
    MZ2)*Cos(alp - beta)*Csc(2*beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2))*((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*(MA02 - MHH2 + MZ2)*v2*C0i(cc2, m22, m32, m12, MHH2, MHH2, 
    MZ2)*Cos(alp - beta)*Csc(2*beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2))*((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.009947183943243459*Alfa*MZ2*C0i(cc0, m22, m12, m32, Mh02, MHH2, MZ2)*
   Cos(alp - beta)*Csc(beta)*(pow(Sin(alp - beta), 2))*Sec(beta)*
   ((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*(MA02 - MHH2 + MZ2)*C0i(cc1, m22, m12, m32, 
    Mh02, MHH2, MZ2)*Cos(alp - beta)*Csc(beta)*(pow(Sin(alp - beta), 2))*
   Sec(beta)*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*MW2*SW*SW2) + (0.0078125*Alfa2*v2*C0i(cc1, m12, m22, m32, MA02, 
    Mh02, MHH2)*((2*MA02 - Mh02)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh02)*Cos(alp + beta))*(pow(Csc(2*beta), 2))*
   (pow(Csc(w), 2))*Sin(2*(alp - beta))*
   ((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0018650969893581485*Alfa*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh02)*Cos(3*alp - beta) + (2*M2 - 3*Mh02)*Cos(alp + beta))*
   (pow(Csc(2*beta), 2))*Re(B0i(bb0, Mh02, Mh02, Mh02))*
   Sin(2*(alp - beta))*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0018650969893581485*Alfa*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh02)*Cos(3*alp - beta) + (2*M2 - 3*Mh02)*Cos(alp + beta))*
   (pow(Csc(2*beta), 2))*Re(B0i(bb0, MHH2, Mh02, Mh02))*
   Sin(2*(alp - beta))*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - (0.0024867959858108648*Alfa*Cos(alp - beta)*
   (pow(Csc(2*beta), 2))*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, Mh02, MHH2))*(9*M2*(Mh02 + MHH2) - 
    (2*Mh02 + MHH2)*(Mh02 + 2*MHH2) + (3*M2 - Mh02 - 2*MHH2)*
     (3*M2 - 2*Mh02 - MHH2)*Cos(4*alp) - 8*(pow(M2, 2)) + 
    M2*(-(M2*Cos(4*beta)) + 2*(Mh02 - MHH2)*Sin(2*alp)*Sin(2*beta))))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.0024867959858108648*Alfa*Cos(alp - beta)*
   (pow(Csc(2*beta), 2))*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, Mh02, MHH2))*(9*M2*(Mh02 + MHH2) - 
    (2*Mh02 + MHH2)*(Mh02 + 2*MHH2) + (3*M2 - Mh02 - 2*MHH2)*
     (3*M2 - 2*Mh02 - MHH2)*Cos(4*alp) - 8*(pow(M2, 2)) + 
    M2*(-(M2*Cos(4*beta)) + 2*(Mh02 - MHH2)*Sin(2*alp)*Sin(2*beta))))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.003730193978716297*Alfa*Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, MHH2, MHH2))*(M2 + (3*M2 - Mh02 - 2*MHH2)*Csc(2*beta)*
     Sin(2*alp))*(M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.003730193978716297*Alfa*Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MHH2, MHH2))*(M2 + (3*M2 - Mh02 - 2*MHH2)*Csc(2*beta)*
     Sin(2*alp))*(M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.015625*Alfa2*v2*C0i(cc1, m12, m22, m32, MA02, MHH2, MHH2)*Csc(2*beta)*
   (pow(Csc(w), 2))*(pow(Sin(alp - beta), 2))*
   (M2 + (3*M2 - Mh02 - 2*MHH2)*Csc(2*beta)*Sin(2*alp))*
   ((2*MA02 - MHH2)*Sin(alp - 3*beta) + (4*M2 - 2*MA02 - 3*MHH2)*
     Sin(alp + beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.015625*Alfa2*v2*C0i(cc2, m32, m22, m12, MA02, Mh02, MHH2)*
   (pow(Cos(alp - beta), 2))*(pow(Csc(2*beta), 2))*
   (pow(Csc(w), 2))*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta))*
   ((2*MA02 - MHH2)*Sin(alp - 3*beta) + (4*M2 - 2*MA02 - 3*MHH2)*
     Sin(alp + beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0024867959858108648*Alfa*(MA02 - MHH2)*B0i(bb0, 0, MA02, MHH2)*
   Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*MA02*MW2*SW*SW2) + 
 (0.0024867959858108648*Alfa*(MA02 - MHH2)*B0i(bb0, MA02, MA02, MHH2)*
   Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*MA02*MW2*SW*SW2) + 
 (0.0078125*Alfa2*(MA02 - Mh02 + MZ2)*v2*C0i(cc0, m12, m32, m22, MA02, 
    MHH2, MZ2)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Sec(beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*(-MA02 + Mh02 + MZ2)*v2*C0i(cc1, m12, m32, m22, MA02, 
    MHH2, MZ2)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Sec(beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*(-MA02 + Mh02 + MZ2)*v2*C0i(cc2, m12, m32, m22, MA02, 
    MHH2, MZ2)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Sec(beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0078125*Alfa2*v2*C0i(cc0, m22, m32, m12, MA02, MA02, MHH2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*(pow(Csc(2*beta), 2))*(pow(Csc(w), 2))*
   Sin(alp - beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0078125*Alfa2*v2*C0i(cc1, m22, m32, m12, MA02, MA02, MHH2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*(pow(Csc(2*beta), 2))*(pow(Csc(w), 2))*
   Sin(alp - beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0078125*Alfa2*v2*C0i(cc2, m22, m32, m12, MA02, MA02, MHH2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*(pow(Csc(2*beta), 2))*(pow(Csc(w), 2))*
   Sin(alp - beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.00015542474911317905*Alfa*((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh02, MA02, MA02))*Sin(alp - beta)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.00015542474911317905*Alfa*((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MA02, MA02))*Sin(alp - beta)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0003108494982263581*Alfa*((Mh02 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh02 + 2*MHp2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh02, MHp2, MHp2))*Sin(alp - beta)*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.0003108494982263581*Alfa*((Mh02 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh02 + 2*MHp2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MHp2, MHp2))*Sin(alp - beta)*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MB2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MB2, MB2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MC2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MC2, MC2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MD2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MD2, MD2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*Csc(beta)*(pow(ME2, 2))*
   (pow(Cos(alp), 3))*Re(B0i(dbb0, Mh02, ME2, ME2))*
   (Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*Csc(beta)*(pow(ML2, 2))*
   (pow(Cos(alp), 3))*Re(B0i(dbb0, Mh02, ML2, ML2))*
   (Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*Csc(beta)*(pow(MM2, 2))*
   (pow(Cos(alp), 3))*Re(B0i(dbb0, Mh02, MM2, MM2))*
   (Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MS2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MS2, MS2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MT2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MT2, MT2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MU2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MU2, MU2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MB2*Mh02*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MB2, MB2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MC2*Mh02*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MC2, MC2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MD2*Mh02*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MD2, MD2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*ME2*Mh02*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, ME2, ME2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*Mh02*ML2*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, ML2, ML2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*Mh02*MM2*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MM2, MM2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*Mh02*MS2*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MS2, MS2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*Mh02*MT2*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MT2, MT2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*Mh02*MU2*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MU2, MU2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2)
;
}

} //end namespace TypeI

namespace TypeII{
ComplexType CAhZ(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32)
{
return 0. - (0.5*Cos(alp - beta))/(CW*SW) + 
 (0.009947183943243459*Alfa*MZ2*B0i(bb0, m12, Mh02, MZ2)*Cos(alp - beta))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*B0i(bb0, m12, MHp2, MW2)*
   Cos(alp - beta))/(CW*SW) + (0.009947183943243459*Alfa*MZ2*
   B0i(bb0, m22, MA02, MZ2)*Cos(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*B0i(bb0, m22, MHp2, MW2)*Cos(alp - beta))/
  (CW*SW) - (0.07957747154594767*Alfa*CW*B0i(bb0, m32, MW2, MW2)*
   Cos(alp - beta))/(SW*SW2) - (0.009947183943243459*Alfa*MZ2*
   B0i(bb1, m12, Mh02, MZ2)*Cos(alp - beta))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*B0i(bb1, m12, MHp2, MW2)*Cos(alp - beta))/
  (CW*SW) - (0.009947183943243459*Alfa*MZ2*B0i(bb1, m22, MA02, MZ2)*
   Cos(alp - beta))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*B0i(bb1, m22, MHp2, MW2)*Cos(alp - beta))/
  (CW*SW) + (0.07957747154594767*Alfa*CW*C0i(cc00, m22, m32, m12, MHp2, 
    MW2, MW2)*Cos(alp - beta))/(SW*SW2) + 
 (0.07957747154594767*Alfa*(2*MW2 - MZ2)*C0i(cc00, m32, m12, m22, MHp2, 
    MHp2, MW2)*Cos(alp - beta))/(CW*MZ2*SW*SW2) - 
 (0.019894367886486918*Alfa*(MW2*(2*(MA02 - MHp2)*(-Mh02 + MHp2) - 
      (MA02 + Mh02 - 2*MHp2 - m12 + m22 + 2*m32)*MW2) + 
    ((MA02 - MHp2)*(Mh02 - MHp2) + (MA02 + Mh02 - 2*MHp2)*MW2)*MZ2)*
   C0i(cc1, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/
  (MW*MW2*MZ*SW*SW2) + (0.019894367886486918*Alfa*CW*(m12 + 3*m22 - m32)*
   C0i(cc11, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/(SW*SW2) + 
 (0.039788735772973836*Alfa*CW*(2*(m12 + m22) - m32)*
   C0i(cc12, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/(SW*SW2) - 
 (0.039788735772973836*Alfa*(m12 - m22)*(2*MW2 - MZ2)*
   C0i(cc12, m32, m12, m22, MHp2, MHp2, MW2)*Cos(alp - beta))/
  (CW*MZ2*SW*SW2) - (0.019894367886486918*Alfa*
   (MW2*(2*(MA02 - MHp2)*(-Mh02 + MHp2) - (MA02 + Mh02 - 2*MHp2 + m12 - 
        m22 + 2*m32)*MW2) + ((MA02 - MHp2)*(Mh02 - MHp2) + 
      (MA02 + Mh02 - 2*MHp2)*MW2)*MZ2)*C0i(cc2, m22, m32, m12, MHp2, MW2, 
    MW2)*Cos(alp - beta))/(MW*MW2*MZ*SW*SW2) - 
 (0.019894367886486918*Alfa*((MA02 - MHp2)*(-Mh02 + MHp2) + 
    MW2*(m12 - m22 - m32 + MW2))*(2*MW2 - MZ2)*
   C0i(cc2, m32, m12, m22, MHp2, MHp2, MW2)*Cos(alp - beta))/
  (CW*MW2*MZ2*SW*SW2) + (0.019894367886486918*Alfa*CW*(3*m12 + m22 - m32)*
   C0i(cc22, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/(SW*SW2) + 
 (0.019894367886486918*Alfa*(m12 + 3*m22 - m32)*(2*MW2 - MZ2)*
   C0i(cc22, m32, m12, m22, MHp2, MHp2, MW2)*Cos(alp - beta))/
  (CW*MZ2*SW*SW2) + (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MC2, MC2, MC2)*
   Cos(alp)*Cot(beta)*Csc(beta)*(pow(MC2, 2)))/
  (-8*MW*MW2*MZ2*Pi*SW + 8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MT2, MT2, MT2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MT2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MU2, MU2, MU2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MU2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MC2*MZ*MZ2*B0i(bb0, m32, MC2, MC2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MT2*MZ*MZ2*B0i(bb0, m32, MT2, MT2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MU2*MZ*MZ2*B0i(bb0, m32, MU2, MU2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MC2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MC2, MC2, MC2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MT2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MT2, MT2, MT2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MU2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MU2, MU2, MU2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MC2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MC2, MC2, MC2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MT2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MT2, MT2, MT2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*m12*MU2*MZ*MZ2*C0i(cc2, m22, m32, m12, MU2, MU2, MU2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (0.019894367886486918*Alfa*MZ2*C0i(cc0, m22, m32, m12, MHp2, MW2, MW2)*
   Cos(alp - beta)*(Mh02 - MHp2 - ((4*MHp2 + m32)*MW2 + 
      (MA02 - MHp2)*(MW2*(2*MHp2 + MW2) - (MHp2 + MW2)*MZ2 + 
        Mh02*(-2*MW2 + MZ2))*(pow(MW2, -1)))*(pow(-MW2 + MZ2, -1))))/
  (MW*MZ*SW) - (0.039788735772973836*Alfa*MZ2*C0i(cc00, m32, m12, m22, 
    MA02, Mh02, MZ2)*(pow(Cos(alp - beta), 3)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(m12 - m22)*MZ2*C0i(cc12, m32, m12, m22, 
    MA02, Mh02, MZ2)*(pow(Cos(alp - beta), 3)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*(m12 + 3*m22 - m32)*MZ2*
   C0i(cc22, m32, m12, m22, MA02, Mh02, MZ2)*
   (pow(Cos(alp - beta), 3)))/(CW*MW2*SW*SW2) + 
 (0.03125*Alfa2*v2*C0i(cc2, m32, m12, m22, MA02, Mh02, MZ2)*
   (MZ2*(m12 - m22 - m32 + MZ2) + (pow(MA02 - Mh02, 2)))*
   (pow(Cos(alp - beta), 3))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0234375*Alfa2*v2*C0i(cc1, m12, m22, m32, MA02, Mh02, Mh02)*
   Cos(alp - beta)*(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh02)*Cos(alp + beta))*((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta))*(pow(Csc(2*beta), 2))*
   (pow(Csc(w), 2)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*v2*C0i(cc0, m22, m32, m12, MA02, MA02, Mh02)*
   Cos(alp - beta)*(pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
      (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta), 2))*
   (pow(Csc(2*beta), 2))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*v2*C0i(cc1, m22, m32, m12, MA02, MA02, Mh02)*
   Cos(alp - beta)*(pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
      (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta), 2))*
   (pow(Csc(2*beta), 2))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*v2*C0i(cc2, m22, m32, m12, MA02, MA02, Mh02)*
   Cos(alp - beta)*(pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
      (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta), 2))*
   (pow(Csc(2*beta), 2))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0049735919716217296*Alfa*Mh02*(MA02 - Mh02 + MZ2)*B0i(bb0, 0, Mh02, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MA02*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*MHH2*(MA02 - MHH2 + MZ2)*B0i(bb0, 0, MHH2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MA02*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*(Mh02*(-Mh02 + MZ2) + MA02*(Mh02 + MZ2))*
   B0i(bb0, MA02, Mh02, MZ2)*Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/
  (CW*MA02*MW2*SW*SW2) + (0.0049735919716217296*Alfa*
   (MHH2*(-MHH2 + MZ2) + MA02*(MHH2 + MZ2))*B0i(bb0, MA02, MHH2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MA02*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MZ2*B0i(bb1, MA02, Mh02, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MZ2*B0i(bb1, MA02, MHH2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*MZ2*(-MA02 + Mh02 + MZ2)*
   C0i(cc0, m12, m22, m32, Mh02, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*MZ2*C0i(cc00, m12, m22, m32, Mh02, MZ2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.039788735772973836*Alfa*MZ2*C0i(cc00, m12, m22, m32, MHH2, MZ2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*MZ2*C0i(cc00, m32, m12, m22, MA02, MHH2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*(3*m12 + m22 - m32)*MZ2*
   C0i(cc11, m12, m22, m32, Mh02, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*(3*m12 + m22 - m32)*MZ2*
   C0i(cc11, m12, m22, m32, MHH2, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*(m12 - m22)*MZ2*C0i(cc12, m12, m22, m32, 
    Mh02, MZ2, MZ2)*Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*(m12 - m22)*MZ2*
   C0i(cc12, m12, m22, m32, MHH2, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(m12 - m22)*MZ2*C0i(cc12, m32, m12, m22, 
    MA02, MHH2, MZ2)*Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2) - (0.009947183943243459*Alfa*(m12 + 3*m22 - m32)*MZ2*
   C0i(cc22, m32, m12, m22, MA02, MHH2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*C0i(cc1, m12, m22, m32, Mh02, MZ2, MZ2)*
   Cos(alp - beta)*((MA02 - Mh02)*Mh02 + (MA02 - Mh02 + m12 - m22 + m32)*
     MZ2 + (pow(MZ2, 2)))*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2) - (0.0625*Alfa2*MZ2*(-MA02 + MHH2 + MZ2)*v2*
   C0i(cc0, m12, m22, m32, MHH2, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Csc(w), 2))*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*(-(Mh02*MHH2) + MA02*(Mh02 + MZ2) + 
    MZ2*(-MHH2 + m12 - m22 + m32 + MZ2))*v2*C0i(cc1, m12, m22, m32, 
    MHH2, MZ2, MZ2)*Cos(alp - beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*v2*C0i(cc2, m32, m12, m22, MA02, MHH2, MZ2)*
   Cos(alp - beta)*(Mh02*MHH2 - MA02*(Mh02 + MHH2) + 
    MZ2*(m12 - m22 - m32 + MZ2) + (pow(MA02, 2)))*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-46*MW2*MZ2 + 20*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(A0i(aa0, MB2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-70*MW2*MZ2 + 56*(pow(MW2, 2)) + 23*(pow(MZ2, 2)))*
   Re(A0i(aa0, MC2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-46*MW2*MZ2 + 20*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(A0i(aa0, MD2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-22*MW2*MZ2 + 12*(pow(MW2, 2)) + 9*(pow(MZ2, 2)))*
   Re(A0i(aa0, ME2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*Re(A0i(aa0, MHp2)))/
  (CW*MW2*MZ2*SW*SW2) - (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-22*MW2*MZ2 + 12*(pow(MW2, 2)) + 9*(pow(MZ2, 2)))*
   Re(A0i(aa0, ML2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-22*MW2*MZ2 + 12*(pow(MW2, 2)) + 9*(pow(MZ2, 2)))*
   Re(A0i(aa0, MM2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-46*MW2*MZ2 + 20*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(A0i(aa0, MS2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-70*MW2*MZ2 + 56*(pow(MW2, 2)) + 23*(pow(MZ2, 2)))*
   Re(A0i(aa0, MT2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-70*MW2*MZ2 + 56*(pow(MW2, 2)) + 23*(pow(MZ2, 2)))*
   Re(A0i(aa0, MU2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(4*MW2 - 3*MZ2)*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(A0i(aa0, MW2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*Re(A0i(aa0, MZ2)))/
  (CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) + 
 (0.17904931097838225*Alfa*Cos(alp - beta)*Re(B0i(bb0, 0, MW2, MW2)))/
  (CW*SW) - (0.0049735919716217296*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb0, MA02, Mh02, MZ2)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MA02, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.009947183943243459*Alfa*Cos(alp - beta)*
   Re(B0i(bb0, MA02, MHp2, MW2)))/(CW*SW*SW2) + 
 (0.0024867959858108648*Alfa*Cos(alp - beta)*(Mh02*MHH2 - 2*Mh02*MZ2 + 
    MHH2*MZ2 - MA02*(Mh02 + MHH2 + MZ2) + Cos(2*(alp - beta))*
     (MHH2*(-Mh02 + MZ2) + MA02*(Mh02 + MHH2 + MZ2) - (pow(MA02, 2))) + 
    (pow(MA02, 2)))*Re(B0i(bb0, Mh02, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.0049735919716217296*Alfa*Cos(alp - beta)*
   ((Mh02 - MHp2)*(MHH2 - MHp2) + (-2*Mh02 + MHH2 - MHp2)*MW2 + 
    ((Mh02 - MHp2)*(-MHH2 + MHp2) + (MHH2 + MHp2)*MW2)*Cos(2*(alp - beta)))*
   Re(B0i(bb0, Mh02, MHp2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*(Mh02*MHH2 + 2*MW2*(-MHH2 + 6*MW2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, MW2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.0024867959858108648*Alfa*(Mh02*MHH2 + 2*MZ2*(-MHH2 + 6*MZ2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*Cos(alp - beta)*(MHH2*(-Mh02 + MZ2) + 
    MA02*(Mh02 + MHH2 + MZ2) - (pow(MA02, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MHH2, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*((Mh02 - MHp2)*(MHH2 - MHp2) - (MHH2 + MHp2)*MW2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MHp2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*(Mh02*MHH2 + 2*MW2*(-MHH2 + 6*MW2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MW2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0024867959858108648*Alfa*(Mh02*MHH2 + 2*MZ2*(-MHH2 + 6*MZ2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.07957747154594767*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, 0, MW2)))/(CW*MZ2*SW*SW2) + 
 (0.05968310365946075*Alfa*(MT2 - MW2)*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, MB2, MT2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.05968310365946075*Alfa*MC2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, MC2, MS2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.05968310365946075*Alfa*(MU2 - MW2)*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, MD2, MU2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MW2, Mh02, MW2)))/
  (CW*MZ2*SW*(pow(SW2, 2))) - (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   (pow(Cos(alp - beta), 3))*Re(B0i(bb0, MW2, MHH2, MW2)))/
  (CW*MZ2*SW*(pow(SW2, 2))) + (0.019894367886486918*Alfa*Cos(alp - beta)*
   (6*MZ2*(pow(MW2, 2)) + 4*(pow(MW2, 3)) - 6*MW2*(pow(MZ2, 2)) + 
    (pow(MZ2, 3)))*Re(B0i(bb0, MW2, MW2, MZ2)))/
  (CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MB2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MB2, MB2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MC2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MC2, MC2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MD2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MD2, MD2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*ME2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, ME2, ME2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MZ2, Mh02, MZ2)))/
  (CW*MW2*SW*(pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   (pow(Cos(alp - beta), 3))*Re(B0i(bb0, MZ2, MHH2, MZ2)))/
  (CW*MW2*SW*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*ML2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, ML2, ML2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*MM2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MM2, MM2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MS2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MS2, MS2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MT2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MT2, MT2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MU2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MU2, MU2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-8*MW2*MZ2 + 13*(pow(MW2, 2)) + 2*(pow(MZ2, 2)))*
   Re(B0i(bb0, MZ2, MW2, MW2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.026525823848649224*Alfa*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MB2, MB2)))/(CW*MW2*MZ2*SW) - 
 (0.05305164769729845*Alfa*(8*MW2 - 5*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MC2, MC2)))/(CW*MW2*MZ2*SW) - 
 (0.026525823848649224*Alfa*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MD2, MD2)))/(CW*MW2*MZ2*SW) - 
 (0.07957747154594767*Alfa*(4*MW2 - 3*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, ME2, ME2)))/(CW*MW2*MZ2*SW) + 
 (0.07957747154594767*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MHp2, MHp2)))/(CW*MW2*MZ2*SW) - 
 (0.07957747154594767*Alfa*(4*MW2 - 3*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, ML2, ML2)))/(CW*MW2*MZ2*SW) - 
 (0.07957747154594767*Alfa*(4*MW2 - 3*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MM2, MM2)))/(CW*MW2*MZ2*SW) - 
 (0.026525823848649224*Alfa*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MS2, MS2)))/(CW*MW2*MZ2*SW) - 
 (0.05305164769729845*Alfa*(8*MW2 - 5*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MT2, MT2)))/(CW*MW2*MZ2*SW) - 
 (0.05305164769729845*Alfa*(8*MW2 - 5*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MU2, MU2)))/(CW*MW2*MZ2*SW) + 
 (0.07957747154594767*Alfa*(6*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MW2, MW2)))/(CW*MW2*MZ2*SW) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, ME2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, ML2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, MM2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.15915494309189535*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, MW2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MA02, MHp2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.1193662073189215*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MB2, MT2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.1193662073189215*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MC2, MS2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.1193662073189215*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MD2, MU2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MW2, Mh02, MHp2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MW2, Mh02, MW2)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MW2, MHH2, MHp2)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MW2, MHH2, MW2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(8*MW2 + MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MW2, MZ2)))/(CW*MW2*SW*(pow(MZ2, 2))*
   (pow(SW2, 2))) + (0.05968310365946075*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*Re(B0i(bb00, MZ2, 0, 0)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MZ2, MA02, Mh02)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MZ2, MA02, MHH2)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.006631455962162306*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MB2, MB2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.006631455962162306*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 
    17*(pow(MZ2, 2)))*Re(B0i(bb00, MZ2, MC2, MC2)))/
  (CW*MW2*SW*(pow(MZ2, 3))*(pow(SW2, 2))) + 
 (0.006631455962162306*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MD2, MD2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, ME2, ME2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) - (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb00, MZ2, Mh02, MZ2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MZ2, MHH2, MZ2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*(pow(2*MW2 - MZ2, 3))*
   Re(B0i(bb00, MZ2, MHp2, MHp2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, ML2, ML2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MM2, MM2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.006631455962162306*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MS2, MS2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.006631455962162306*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 
    17*(pow(MZ2, 2)))*Re(B0i(bb00, MZ2, MT2, MT2)))/
  (CW*MW2*SW*(pow(MZ2, 3))*(pow(SW2, 2))) + 
 (0.006631455962162306*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MU2, MU2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) - (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-4*MW2*MZ2 + 12*(pow(MW2, 2)) + (pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MW2, MW2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.026525823848649224*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MB2, MB2)))/(CW*SW) + 
 (0.1061032953945969*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MC2, MC2)))/
  (CW*SW) + (0.026525823848649224*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MD2, MD2)))/(CW*SW) + 
 (0.07957747154594767*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, ME2, ME2)))/
  (CW*SW) + (0.07957747154594767*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, ML2, ML2)))/(CW*SW) + 
 (0.07957747154594767*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MM2, MM2)))/
  (CW*SW) + (0.026525823848649224*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MS2, MS2)))/(CW*SW) + 
 (0.1061032953945969*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MT2, MT2)))/
  (CW*SW) + (0.1061032953945969*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MU2, MU2)))/(CW*SW) + 
 (0.039788735772973836*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MW2, MW2)))/
  (CW*SW) - (0.05968310365946075*Alfa*MB2*Cos(alp - beta)*
   (pow(Tan(beta), 2))*Re(B0i(bb1, MA02, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MC2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MD2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(bb1, MA02, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*ME2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(bb1, MA02, ME2, ME2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb1, MA02, Mh02, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MZ2*Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb1, MA02, MHH2, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*Re(B0i(bb1, MA02, MHp2, MW2)))/
  (CW*SW*SW2) - (0.019894367886486918*Alfa*ML2*Cos(alp - beta)*
   (pow(Tan(beta), 2))*Re(B0i(bb1, MA02, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MM2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(bb1, MA02, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MS2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(bb1, MA02, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MT2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MU2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MU2, MU2)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*MZ2*Cos(alp - beta)*
   (-2*Mh02 + MHH2 + MHH2*Cos(2*(alp - beta)))*Re(B0i(bb1, Mh02, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-2*Mh02 + MHH2 + MHH2*Cos(2*(alp - beta)))*Re(B0i(bb1, Mh02, MHp2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) - (0.019894367886486918*Alfa*MHH2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, Mh02, MW2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) - (0.009947183943243459*Alfa*MHH2*MZ2*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb1, Mh02, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MHH2*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, MHH2, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MHH2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, MHH2, MHp2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) + (0.019894367886486918*Alfa*MHH2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, MHH2, MW2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) + (0.009947183943243459*Alfa*MHH2*MZ2*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb1, MHH2, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, ME2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, ML2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, MM2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, MW2)))/(CW*MZ2*SW*SW2) - 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MB2, MT2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MC2, MS2)))/(CW*MZ2*SW*(pow(SW2, 2))) - 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MD2, MU2)))/(CW*MZ2*SW*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*MW2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MW2, MZ2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*Cos(alp - beta)*Re(B0i(bb1, MZ2, 0, 0)))/
  (CW*SW*(pow(SW2, 2))) - (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MC2, MC2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, ME2, ME2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, ML2, ML2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MM2, MM2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MT2, MT2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MU2, MU2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*Cos(alp - beta)*(pow(MW2, 2))*
   Re(B0i(bb1, MZ2, MW2, MW2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) + 
 (0.0012433979929054324*Alfa*Cos(alp - beta)*
   (pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*Mh02)*
       Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, MA02, MA02, Mh02)))/(CW*MW2*SW*SW2) + 
 (0.0012433979929054324*Alfa*Cos(alp - beta)*(pow(Csc(2*beta), 2))*
   (pow((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
       Sin(alp + beta), 2))*Re(B0i(dbb0, MA02, MA02, MHH2)))/
  (CW*MW2*SW*SW2) - (0.0049735919716217296*Alfa*
   ((MA02 + Mh02)*MZ2 - (pow(MA02 - Mh02, 2)))*
   (pow(Cos(alp - beta), 3))*Re(B0i(dbb0, MA02, Mh02, MZ2)))/
  (CW*MW2*SW*SW2) - (0.0049735919716217296*Alfa*Cos(alp - beta)*
   ((MA02 + MHH2)*MZ2 - (pow(MA02 - MHH2, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, MA02, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) + (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-((MA02 + MHp2)*MW2) + (pow(MA02 - MHp2, 2)))*
   Re(B0i(dbb0, MA02, MHp2, MW2)))/(CW*MW2*SW*SW2) + 
 (0.0006216989964527162*Alfa*Cos(alp - beta)*
   (pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*Mh02)*
       Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh02, MA02, MA02)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*((MA02 + Mh02)*MZ2 - (pow(MA02 - Mh02, 2)))*
   (pow(Cos(alp - beta), 3))*Re(B0i(dbb0, Mh02, MA02, MZ2)))/
  (CW*MW2*SW*SW2) + (0.005595290968074445*Alfa*Cos(alp - beta)*
   (pow(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
      (2*M2 - 3*Mh02)*Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh02, Mh02, Mh02)))/(CW*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*(pow(Cos(alp - beta), 3))*
   (pow(M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*Sin(2*alp), 2))*
   Re(B0i(dbb0, Mh02, Mh02, MHH2)))/(CW*MW2*SW*SW2) + 
 (0.0024867959858108648*Alfa*Cos(alp - beta)*(pow(Csc(2*beta), 2))*
   (pow(Sin(alp - beta), 2))*
   (pow((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta), 2))*
   Re(B0i(dbb0, Mh02, MHH2, MHH2)))/(CW*MW2*SW*SW2) + 
 (0.0012433979929054324*Alfa*Cos(alp - beta)*
   (pow((Mh02 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh02 + 2*MHp2)*
       Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh02, MHp2, MHp2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*(-((Mh02 + MHp2)*MW2) + 
    (pow(Mh02 - MHp2, 2)))*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb0, Mh02, MHp2, MW2)))/(CW*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*Cos(alp - beta)*
   (-2*Mh02*MW2 + (pow(Mh02, 2)) + 12*(pow(MW2, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, Mh02, MW2, MW2)))/
  (CW*MW2*SW*SW2) + (0.0024867959858108648*Alfa*Cos(alp - beta)*
   (-2*Mh02*MZ2 + (pow(Mh02, 2)) + 12*(pow(MZ2, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, Mh02, MZ2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.029841551829730376*Alfa*MB2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MC2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MD2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*ME2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, ME2, ME2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*(pow(MZ2, 2))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, MZ2, Mh02, MZ2)))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*(pow(MZ2, 2))*
   (pow(Cos(alp - beta), 3))*Re(B0i(dbb0, MZ2, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.009947183943243459*Alfa*ML2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MM2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MS2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MT2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MU2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MU2, MU2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*(9*MW2 - 2*MZ2)*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MW2, MW2)))/(CW*SW*SW2) - 
 (0.05305164769729845*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MB2, MB2)))/
  (CW*SW) - (0.2122065907891938*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, MC2, MC2)))/(CW*SW) - 
 (0.05305164769729845*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MD2, MD2)))/
  (CW*SW) - (0.15915494309189535*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, ME2, ME2)))/(CW*SW) + 
 (0.07957747154594767*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MHp2, MHp2)))/
  (CW*SW) - (0.15915494309189535*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, ML2, ML2)))/(CW*SW) - 
 (0.15915494309189535*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MM2, MM2)))/
  (CW*SW) - (0.05305164769729845*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, MS2, MS2)))/(CW*SW) - 
 (0.2122065907891938*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MT2, MT2)))/
  (CW*SW) - (0.2122065907891938*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, MU2, MU2)))/(CW*SW) + 
 (0.238732414637843*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MW2, MW2)))/
  (CW*SW) + (0.05968310365946075*Alfa*MZ2*Cos(alp - beta)*
   Re(B0i(dbb00, MZ2, 0, 0)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb00, MZ2, MA02, Mh02)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(dbb00, MZ2, MA02, MHH2)))/(CW*MW2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MB2, MB2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MC2, MC2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MD2, MD2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, ME2, ME2)))/(CW*MW2*MZ2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(dbb00, MZ2, Mh02, MZ2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb00, MZ2, MHH2, MZ2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*(pow(-2*MW2 + MZ2, 2))*
   Re(B0i(dbb00, MZ2, MHp2, MHp2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, ML2, ML2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MM2, MM2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MS2, MS2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MT2, MT2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MU2, MU2)))/(CW*MW2*MZ2*SW*SW2) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 12*(pow(MW2, 2)) + (pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MW2, MW2)))/(CW*MW2*MZ2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MB2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(dbb1, MA02, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MC2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MD2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(dbb1, MA02, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MA02*ME2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(dbb1, MA02, ME2, ME2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MA02*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb1, MA02, Mh02, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MA02*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb1, MA02, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*MA02*Cos(alp - beta)*
   Re(B0i(dbb1, MA02, MHp2, MW2)))/(CW*SW*SW2) - 
 (0.019894367886486918*Alfa*MA02*ML2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(dbb1, MA02, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MA02*MM2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(dbb1, MA02, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MS2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(dbb1, MA02, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MT2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MU2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MU2, MU2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*Mh02*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb1, Mh02, MA02, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*Mh02*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb1, Mh02, MHp2, MW2)))/(CW*SW*SW2) + 
 (0.019894367886486918*Alfa*Mh02*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb1, Mh02, MW2, MW2)))/(CW*SW*SW2) + 
 (0.009947183943243459*Alfa*Mh02*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb1, Mh02, MZ2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.029841551829730376*Alfa*Cos(alp - beta)*
   (pow(MZ2, 2))*Re(B0i(dbb1, MZ2, 0, 0)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, ME2, ME2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MU2, MU2)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*MW2*Cos(alp - beta)*Re(B0i(dbb1, MZ2, MW2, MW2)))/
  (CW*SW*SW2) + (0.05968310365946075*Alfa*MC2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MC2, MC2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MT2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MT2, MT2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MU2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MU2, MU2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*(pow(MC2, 2))*
   Re(B0i(bb0, Mh02, MC2, MC2))*Sin(alp)*(Cos(alp) - Cot(beta)*Sin(alp)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*
   (pow(MT2, 2))*Re(B0i(bb0, Mh02, MT2, MT2))*Sin(alp)*
   (Cos(alp) - Cot(beta)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*(pow(MU2, 2))*
   Re(B0i(bb0, Mh02, MU2, MU2))*Sin(alp)*(Cos(alp) - Cot(beta)*Sin(alp)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.03125*Alfa2*(-MA02 + Mh02 + MZ2)*v2*C0i(cc0, m22, m32, m12, Mh02, 
    MHH2, MZ2)*(pow(Cos(alp - beta), 3))*(pow(Csc(w), 2))*
   (M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*Sin(2*alp)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.03125*Alfa2*(MA02 - Mh02 + MZ2)*v2*C0i(cc1, m22, m32, m12, Mh02, MHH2, 
    MZ2)*(pow(Cos(alp - beta), 3))*(pow(Csc(w), 2))*
   (M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*Sin(2*alp)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.03125*Alfa2*(MA02 - Mh02 + MZ2)*v2*C0i(cc2, m22, m32, m12, Mh02, MHH2, 
    MZ2)*(pow(Cos(alp - beta), 3))*(pow(Csc(w), 2))*
   (M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*Sin(2*alp)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.05968310365946075*Alfa*MC2*B0i(bb1, MA02, MC2, MC2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MT2*B0i(bb1, MA02, MT2, MT2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MU2*B0i(bb1, MA02, MU2, MU2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*(pow(MB2, 2))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, Mh02, MB2, MB2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.05968310365946075*Alfa*(pow(MD2, 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh02, MD2, MD2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*(pow(MS2, 2))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, Mh02, MS2, MS2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - (0.05968310365946075*Alfa*(pow(MB2, 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MB2, MB2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*(pow(MC2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MC2, MC2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - (0.05968310365946075*Alfa*(pow(MD2, 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MD2, MD2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*(pow(ME2, 2))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, ME2, ME2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*(pow(ML2, 2))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, ML2, ML2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*(pow(MM2, 2))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, MM2, MM2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - (0.05968310365946075*Alfa*(pow(MS2, 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MS2, MS2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*(pow(MT2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MT2, MT2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.05968310365946075*Alfa*(pow(MU2, 2))*
   (pow(Csc(beta), 2))*Re(B0i(bb0, MHH2, MU2, MU2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MB2*MHH2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, MHH2, MB2, MB2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MC2*MHH2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MC2, MC2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MD2*MHH2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, MHH2, MD2, MD2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*ME2*MHH2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, MHH2, ME2, ME2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MHH2*ML2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, MHH2, ML2, ML2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MHH2*MM2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, MHH2, MM2, MM2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MHH2*MS2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, MHH2, MS2, MS2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MHH2*MT2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MT2, MT2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MHH2*MU2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MU2, MU2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.0012433979929054324*Alfa*(MA02 - Mh02)*
   B0i(bb0, 0, MA02, Mh02)*((2*MA02 - Mh02)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh02)*Cos(alp + beta))*Csc(2*beta)*
   Sin(2*(alp - beta)))/(CW*MA02*MW2*SW*SW2) - 
 (0.0012433979929054324*Alfa*(MA02 - Mh02)*B0i(bb0, MA02, MA02, Mh02)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(2*beta)*Sin(2*(alp - beta)))/(CW*MA02*MW2*SW*SW2) + 
 (0.003730193978716297*Alfa*(-MA02 + Mh02 + MZ2)*C0i(cc0, m22, m32, m12, 
    Mh02, Mh02, MZ2)*(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh02)*Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) - (0.003730193978716297*Alfa*(MA02 - Mh02 + MZ2)*
   C0i(cc1, m22, m32, m12, Mh02, Mh02, MZ2)*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh02)*Cos(3*alp - beta) + (2*M2 - 3*Mh02)*Cos(alp + beta))*
   Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/(CW*MW2*SW*SW2) - 
 (0.003730193978716297*Alfa*(MA02 - Mh02 + MZ2)*C0i(cc2, m22, m32, m12, 
    Mh02, Mh02, MZ2)*(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh02)*Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) + (0.0012433979929054324*Alfa*(MA02 - Mh02 + MZ2)*
   C0i(cc0, m12, m32, m22, MA02, Mh02, MZ2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) - (0.0012433979929054324*Alfa*(-MA02 + Mh02 + MZ2)*
   C0i(cc1, m12, m32, m22, MA02, Mh02, MZ2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) - (0.0012433979929054324*Alfa*(-MA02 + Mh02 + MZ2)*
   C0i(cc2, m12, m32, m22, MA02, Mh02, MZ2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) + (0.05968310365946075*Alfa*MB2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, Mh02, MB2, MB2))*Sin(alp)*(MHH2*Cos(alp - beta)*Sin(alp) - 
    Mh02*Sin(beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MD2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, Mh02, MD2, MD2))*Sin(alp)*(MHH2*Cos(alp - beta)*Sin(alp) - 
    Mh02*Sin(beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*ME2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, Mh02, ME2, ME2))*Sin(alp)*(MHH2*Cos(alp - beta)*Sin(alp) - 
    Mh02*Sin(beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*ML2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, Mh02, ML2, ML2))*Sin(alp)*(MHH2*Cos(alp - beta)*Sin(alp) - 
    Mh02*Sin(beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*MM2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, Mh02, MM2, MM2))*Sin(alp)*(MHH2*Cos(alp - beta)*Sin(alp) - 
    Mh02*Sin(beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MS2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, Mh02, MS2, MS2))*Sin(alp)*(MHH2*Cos(alp - beta)*Sin(alp) - 
    Mh02*Sin(beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.03125*Alfa2*(-MA02 + MHH2 + MZ2)*v2*C0i(cc0, m22, m32, m12, MHH2, 
    MHH2, MZ2)*Cos(alp - beta)*Csc(2*beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2))*((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*(MA02 - MHH2 + MZ2)*v2*C0i(cc1, m22, m32, m12, MHH2, MHH2, 
    MZ2)*Cos(alp - beta)*Csc(2*beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2))*((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*(MA02 - MHH2 + MZ2)*v2*C0i(cc2, m22, m32, m12, MHH2, MHH2, 
    MZ2)*Cos(alp - beta)*Csc(2*beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2))*((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.009947183943243459*Alfa*MZ2*C0i(cc0, m22, m12, m32, Mh02, MHH2, MZ2)*
   Cos(alp - beta)*Csc(beta)*(pow(Sin(alp - beta), 2))*Sec(beta)*
   ((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*(MA02 - MHH2 + MZ2)*C0i(cc1, m22, m12, m32, 
    Mh02, MHH2, MZ2)*Cos(alp - beta)*Csc(beta)*(pow(Sin(alp - beta), 2))*
   Sec(beta)*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*MW2*SW*SW2) + (0.0078125*Alfa2*v2*C0i(cc1, m12, m22, m32, MA02, 
    Mh02, MHH2)*((2*MA02 - Mh02)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh02)*Cos(alp + beta))*(pow(Csc(2*beta), 2))*
   (pow(Csc(w), 2))*Sin(2*(alp - beta))*
   ((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0018650969893581485*Alfa*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh02)*Cos(3*alp - beta) + (2*M2 - 3*Mh02)*Cos(alp + beta))*
   (pow(Csc(2*beta), 2))*Re(B0i(bb0, Mh02, Mh02, Mh02))*
   Sin(2*(alp - beta))*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0018650969893581485*Alfa*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh02)*Cos(3*alp - beta) + (2*M2 - 3*Mh02)*Cos(alp + beta))*
   (pow(Csc(2*beta), 2))*Re(B0i(bb0, MHH2, Mh02, Mh02))*
   Sin(2*(alp - beta))*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - (0.0024867959858108648*Alfa*Cos(alp - beta)*
   (pow(Csc(2*beta), 2))*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, Mh02, MHH2))*(9*M2*(Mh02 + MHH2) - 
    (2*Mh02 + MHH2)*(Mh02 + 2*MHH2) + (3*M2 - Mh02 - 2*MHH2)*
     (3*M2 - 2*Mh02 - MHH2)*Cos(4*alp) - 8*(pow(M2, 2)) + 
    M2*(-(M2*Cos(4*beta)) + 2*(Mh02 - MHH2)*Sin(2*alp)*Sin(2*beta))))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.0024867959858108648*Alfa*Cos(alp - beta)*
   (pow(Csc(2*beta), 2))*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, Mh02, MHH2))*(9*M2*(Mh02 + MHH2) - 
    (2*Mh02 + MHH2)*(Mh02 + 2*MHH2) + (3*M2 - Mh02 - 2*MHH2)*
     (3*M2 - 2*Mh02 - MHH2)*Cos(4*alp) - 8*(pow(M2, 2)) + 
    M2*(-(M2*Cos(4*beta)) + 2*(Mh02 - MHH2)*Sin(2*alp)*Sin(2*beta))))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.003730193978716297*Alfa*Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, MHH2, MHH2))*(M2 + (3*M2 - Mh02 - 2*MHH2)*Csc(2*beta)*
     Sin(2*alp))*(M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.003730193978716297*Alfa*Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MHH2, MHH2))*(M2 + (3*M2 - Mh02 - 2*MHH2)*Csc(2*beta)*
     Sin(2*alp))*(M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.015625*Alfa2*v2*C0i(cc1, m12, m22, m32, MA02, MHH2, MHH2)*Csc(2*beta)*
   (pow(Csc(w), 2))*(pow(Sin(alp - beta), 2))*
   (M2 + (3*M2 - Mh02 - 2*MHH2)*Csc(2*beta)*Sin(2*alp))*
   ((2*MA02 - MHH2)*Sin(alp - 3*beta) + (4*M2 - 2*MA02 - 3*MHH2)*
     Sin(alp + beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.015625*Alfa2*v2*C0i(cc2, m32, m22, m12, MA02, Mh02, MHH2)*
   (pow(Cos(alp - beta), 2))*(pow(Csc(2*beta), 2))*
   (pow(Csc(w), 2))*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta))*
   ((2*MA02 - MHH2)*Sin(alp - 3*beta) + (4*M2 - 2*MA02 - 3*MHH2)*
     Sin(alp + beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0024867959858108648*Alfa*(MA02 - MHH2)*B0i(bb0, 0, MA02, MHH2)*
   Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*MA02*MW2*SW*SW2) + 
 (0.0024867959858108648*Alfa*(MA02 - MHH2)*B0i(bb0, MA02, MA02, MHH2)*
   Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*MA02*MW2*SW*SW2) + 
 (0.0078125*Alfa2*(MA02 - Mh02 + MZ2)*v2*C0i(cc0, m12, m32, m22, MA02, 
    MHH2, MZ2)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Sec(beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*(-MA02 + Mh02 + MZ2)*v2*C0i(cc1, m12, m32, m22, MA02, 
    MHH2, MZ2)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Sec(beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*(-MA02 + Mh02 + MZ2)*v2*C0i(cc2, m12, m32, m22, MA02, 
    MHH2, MZ2)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Sec(beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0078125*Alfa2*v2*C0i(cc0, m22, m32, m12, MA02, MA02, MHH2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*(pow(Csc(2*beta), 2))*(pow(Csc(w), 2))*
   Sin(alp - beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0078125*Alfa2*v2*C0i(cc1, m22, m32, m12, MA02, MA02, MHH2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*(pow(Csc(2*beta), 2))*(pow(Csc(w), 2))*
   Sin(alp - beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0078125*Alfa2*v2*C0i(cc2, m22, m32, m12, MA02, MA02, MHH2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*(pow(Csc(2*beta), 2))*(pow(Csc(w), 2))*
   Sin(alp - beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.00015542474911317905*Alfa*((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh02, MA02, MA02))*Sin(alp - beta)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.00015542474911317905*Alfa*((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MA02, MA02))*Sin(alp - beta)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0003108494982263581*Alfa*((Mh02 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh02 + 2*MHp2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh02, MHp2, MHp2))*Sin(alp - beta)*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.0003108494982263581*Alfa*((Mh02 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh02 + 2*MHp2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MHp2, MHp2))*Sin(alp - beta)*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MC2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MC2, MC2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MT2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MT2, MT2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MU2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MU2, MU2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MC2*Mh02*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MC2, MC2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*Mh02*MT2*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MT2, MT2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*Mh02*MU2*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MU2, MU2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MB2, MB2, MB2)*(pow(MB2, 2))*
   Sec(beta)*Sin(alp)*Tan(beta))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MD2, MD2, MD2)*(pow(MD2, 2))*
   Sec(beta)*Sin(alp)*Tan(beta))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, ME2, ME2, ME2)*(pow(ME2, 2))*
   Sec(beta)*Sin(alp)*Tan(beta))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, ML2, ML2, ML2)*(pow(ML2, 2))*
   Sec(beta)*Sin(alp)*Tan(beta))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MM2, MM2, MM2)*(pow(MM2, 2))*
   Sec(beta)*Sin(alp)*Tan(beta))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MS2, MS2, MS2)*(pow(MS2, 2))*
   Sec(beta)*Sin(alp)*Tan(beta))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MB2*MZ*MZ2*B0i(bb0, m32, MB2, MB2)*Sec(beta)*Sin(alp)*Tan(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MD2*MZ*MZ2*B0i(bb0, m32, MD2, MD2)*Sec(beta)*Sin(alp)*Tan(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ME2*MZ*MZ2*B0i(bb0, m32, ME2, ME2)*Sec(beta)*Sin(alp)*Tan(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ML2*MZ*MZ2*B0i(bb0, m32, ML2, ML2)*Sec(beta)*Sin(alp)*Tan(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MM2*MZ*MZ2*B0i(bb0, m32, MM2, MM2)*Sec(beta)*Sin(alp)*Tan(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MS2*MZ*MZ2*B0i(bb0, m32, MS2, MS2)*Sec(beta)*Sin(alp)*Tan(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MB2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MB2, MB2, MB2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MD2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MD2, MD2, MD2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ME2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, ME2, ME2, ME2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ML2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, ML2, ML2, ML2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MM2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MM2, MM2, MM2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MS2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MS2, MS2, MS2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MB2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MB2, MB2, MB2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MD2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MD2, MD2, MD2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ME2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, ME2, ME2, ME2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ML2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, ML2, ML2, ML2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MM2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MM2, MM2, MM2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MS2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MS2, MS2, MS2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) - 
 (0.05968310365946075*Alfa*MB2*B0i(bb1, MA02, MB2, MB2)*Sin(alp - beta)*
   Tan(beta))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MD2*B0i(bb1, MA02, MD2, MD2)*Sin(alp - beta)*
   Tan(beta))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*ME2*B0i(bb1, MA02, ME2, ME2)*Sin(alp - beta)*
   Tan(beta))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*ML2*B0i(bb1, MA02, ML2, ML2)*Sin(alp - beta)*
   Tan(beta))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MM2*B0i(bb1, MA02, MM2, MM2)*Sin(alp - beta)*
   Tan(beta))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MS2*B0i(bb1, MA02, MS2, MS2)*Sin(alp - beta)*
   Tan(beta))/(CW*MW2*SW*SW2) - (0.1193662073189215*Alfa*(pow(MB2, 2))*
   (pow(Sin(alp), 3))*Re(B0i(dbb0, Mh02, MB2, MB2))*Sec(beta)*
   (Cot(alp) + Tan(beta)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*(pow(MD2, 2))*(pow(Sin(alp), 3))*
   Re(B0i(dbb0, Mh02, MD2, MD2))*Sec(beta)*(Cot(alp) + Tan(beta)))/
  (CW*MW2*SW*SW2) - (0.039788735772973836*Alfa*(pow(ME2, 2))*
   (pow(Sin(alp), 3))*Re(B0i(dbb0, Mh02, ME2, ME2))*Sec(beta)*
   (Cot(alp) + Tan(beta)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*(pow(ML2, 2))*(pow(Sin(alp), 3))*
   Re(B0i(dbb0, Mh02, ML2, ML2))*Sec(beta)*(Cot(alp) + Tan(beta)))/
  (CW*MW2*SW*SW2) - (0.039788735772973836*Alfa*(pow(MM2, 2))*
   (pow(Sin(alp), 3))*Re(B0i(dbb0, Mh02, MM2, MM2))*Sec(beta)*
   (Cot(alp) + Tan(beta)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*(pow(MS2, 2))*(pow(Sin(alp), 3))*
   Re(B0i(dbb0, Mh02, MS2, MS2))*Sec(beta)*(Cot(alp) + Tan(beta)))/
  (CW*MW2*SW*SW2) - (0.05968310365946075*Alfa*MB2*Mh02*(pow(Sin(alp), 3))*
   Re(B0i(dbb1, Mh02, MB2, MB2))*Sec(beta)*(Cot(alp) + Tan(beta)))/
  (CW*MW2*SW*SW2) - (0.05968310365946075*Alfa*MD2*Mh02*(pow(Sin(alp), 3))*
   Re(B0i(dbb1, Mh02, MD2, MD2))*Sec(beta)*(Cot(alp) + Tan(beta)))/
  (CW*MW2*SW*SW2) - (0.019894367886486918*Alfa*ME2*Mh02*
   (pow(Sin(alp), 3))*Re(B0i(dbb1, Mh02, ME2, ME2))*Sec(beta)*
   (Cot(alp) + Tan(beta)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*Mh02*ML2*(pow(Sin(alp), 3))*
   Re(B0i(dbb1, Mh02, ML2, ML2))*Sec(beta)*(Cot(alp) + Tan(beta)))/
  (CW*MW2*SW*SW2) - (0.019894367886486918*Alfa*Mh02*MM2*
   (pow(Sin(alp), 3))*Re(B0i(dbb1, Mh02, MM2, MM2))*Sec(beta)*
   (Cot(alp) + Tan(beta)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*Mh02*MS2*(pow(Sin(alp), 3))*
   Re(B0i(dbb1, Mh02, MS2, MS2))*Sec(beta)*(Cot(alp) + Tan(beta)))/
  (CW*MW2*SW*SW2) + (1.*Re(B0i(bb0, Mh02, ME2, ME2))*
   (Alfa*Cos(alp)*(pow(ME2, 2))*(pow(Sin(alp), 2))*Sec(beta) - 
    Alfa*(pow(ME2, 2))*(pow(Cos(alp), 2))*Sec(beta)*Sin(alp)*
     Tan(beta)))/(8*CW*Mh02*MW2*Pi*SW*SW2 - 8*CW*MHH2*MW2*Pi*SW*SW2) + 
 (1.*Re(B0i(bb0, Mh02, ML2, ML2))*(Alfa*Cos(alp)*(pow(ML2, 2))*
     (pow(Sin(alp), 2))*Sec(beta) - Alfa*(pow(ML2, 2))*
     (pow(Cos(alp), 2))*Sec(beta)*Sin(alp)*Tan(beta)))/
  (8*CW*Mh02*MW2*Pi*SW*SW2 - 8*CW*MHH2*MW2*Pi*SW*SW2) + 
 (1.*Re(B0i(bb0, Mh02, MM2, MM2))*(Alfa*Cos(alp)*(pow(MM2, 2))*
     (pow(Sin(alp), 2))*Sec(beta) - Alfa*(pow(MM2, 2))*
     (pow(Cos(alp), 2))*Sec(beta)*Sin(alp)*Tan(beta)))/
  (8*CW*Mh02*MW2*Pi*SW*SW2 - 8*CW*MHH2*MW2*Pi*SW*SW2)
;
}

} //end namespace TypeII

namespace TypeLS{
ComplexType CAhZ(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32)
{
return 0. - (0.5*Cos(alp - beta))/(CW*SW) + 
 (0.009947183943243459*Alfa*MZ2*B0i(bb0, m12, Mh02, MZ2)*Cos(alp - beta))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*B0i(bb0, m12, MHp2, MW2)*
   Cos(alp - beta))/(CW*SW) + (0.009947183943243459*Alfa*MZ2*
   B0i(bb0, m22, MA02, MZ2)*Cos(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*B0i(bb0, m22, MHp2, MW2)*Cos(alp - beta))/
  (CW*SW) - (0.07957747154594767*Alfa*CW*B0i(bb0, m32, MW2, MW2)*
   Cos(alp - beta))/(SW*SW2) - (0.009947183943243459*Alfa*MZ2*
   B0i(bb1, m12, Mh02, MZ2)*Cos(alp - beta))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*B0i(bb1, m12, MHp2, MW2)*Cos(alp - beta))/
  (CW*SW) - (0.009947183943243459*Alfa*MZ2*B0i(bb1, m22, MA02, MZ2)*
   Cos(alp - beta))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*B0i(bb1, m22, MHp2, MW2)*Cos(alp - beta))/
  (CW*SW) + (0.07957747154594767*Alfa*CW*C0i(cc00, m22, m32, m12, MHp2, 
    MW2, MW2)*Cos(alp - beta))/(SW*SW2) + 
 (0.07957747154594767*Alfa*(2*MW2 - MZ2)*C0i(cc00, m32, m12, m22, MHp2, 
    MHp2, MW2)*Cos(alp - beta))/(CW*MZ2*SW*SW2) - 
 (0.019894367886486918*Alfa*(MW2*(2*(MA02 - MHp2)*(-Mh02 + MHp2) - 
      (MA02 + Mh02 - 2*MHp2 - m12 + m22 + 2*m32)*MW2) + 
    ((MA02 - MHp2)*(Mh02 - MHp2) + (MA02 + Mh02 - 2*MHp2)*MW2)*MZ2)*
   C0i(cc1, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/
  (MW*MW2*MZ*SW*SW2) + (0.019894367886486918*Alfa*CW*(m12 + 3*m22 - m32)*
   C0i(cc11, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/(SW*SW2) + 
 (0.039788735772973836*Alfa*CW*(2*(m12 + m22) - m32)*
   C0i(cc12, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/(SW*SW2) - 
 (0.039788735772973836*Alfa*(m12 - m22)*(2*MW2 - MZ2)*
   C0i(cc12, m32, m12, m22, MHp2, MHp2, MW2)*Cos(alp - beta))/
  (CW*MZ2*SW*SW2) - (0.019894367886486918*Alfa*
   (MW2*(2*(MA02 - MHp2)*(-Mh02 + MHp2) - (MA02 + Mh02 - 2*MHp2 + m12 - 
        m22 + 2*m32)*MW2) + ((MA02 - MHp2)*(Mh02 - MHp2) + 
      (MA02 + Mh02 - 2*MHp2)*MW2)*MZ2)*C0i(cc2, m22, m32, m12, MHp2, MW2, 
    MW2)*Cos(alp - beta))/(MW*MW2*MZ*SW*SW2) - 
 (0.019894367886486918*Alfa*((MA02 - MHp2)*(-Mh02 + MHp2) + 
    MW2*(m12 - m22 - m32 + MW2))*(2*MW2 - MZ2)*
   C0i(cc2, m32, m12, m22, MHp2, MHp2, MW2)*Cos(alp - beta))/
  (CW*MW2*MZ2*SW*SW2) + (0.019894367886486918*Alfa*CW*(3*m12 + m22 - m32)*
   C0i(cc22, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/(SW*SW2) + 
 (0.019894367886486918*Alfa*(m12 + 3*m22 - m32)*(2*MW2 - MZ2)*
   C0i(cc22, m32, m12, m22, MHp2, MHp2, MW2)*Cos(alp - beta))/
  (CW*MZ2*SW*SW2) + (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MB2, MB2, MB2)*
   Cos(alp)*Cot(beta)*Csc(beta)*(pow(MB2, 2)))/
  (-8*MW*MW2*MZ2*Pi*SW + 8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MC2, MC2, MC2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MC2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MD2, MD2, MD2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MD2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MS2, MS2, MS2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MS2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MT2, MT2, MT2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MT2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MU2, MU2, MU2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MU2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MB2*MZ*MZ2*B0i(bb0, m32, MB2, MB2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MC2*MZ*MZ2*B0i(bb0, m32, MC2, MC2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MD2*MZ*MZ2*B0i(bb0, m32, MD2, MD2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MS2*MZ*MZ2*B0i(bb0, m32, MS2, MS2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MT2*MZ*MZ2*B0i(bb0, m32, MT2, MT2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MU2*MZ*MZ2*B0i(bb0, m32, MU2, MU2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MB2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MB2, MB2, MB2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MC2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MC2, MC2, MC2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MD2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MD2, MD2, MD2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MS2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MS2, MS2, MS2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MT2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MT2, MT2, MT2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MU2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MU2, MU2, MU2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MB2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MB2, MB2, MB2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MC2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MC2, MC2, MC2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MD2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MD2, MD2, MD2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MS2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MS2, MS2, MS2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MT2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MT2, MT2, MT2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*m12*MU2*MZ*MZ2*C0i(cc2, m22, m32, m12, MU2, MU2, MU2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (0.019894367886486918*Alfa*MZ2*C0i(cc0, m22, m32, m12, MHp2, MW2, MW2)*
   Cos(alp - beta)*(Mh02 - MHp2 - ((4*MHp2 + m32)*MW2 + 
      (MA02 - MHp2)*(MW2*(2*MHp2 + MW2) - (MHp2 + MW2)*MZ2 + 
        Mh02*(-2*MW2 + MZ2))*(pow(MW2, -1)))*(pow(-MW2 + MZ2, -1))))/
  (MW*MZ*SW) - (0.039788735772973836*Alfa*MZ2*C0i(cc00, m32, m12, m22, 
    MA02, Mh02, MZ2)*(pow(Cos(alp - beta), 3)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(m12 - m22)*MZ2*C0i(cc12, m32, m12, m22, 
    MA02, Mh02, MZ2)*(pow(Cos(alp - beta), 3)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*(m12 + 3*m22 - m32)*MZ2*
   C0i(cc22, m32, m12, m22, MA02, Mh02, MZ2)*
   (pow(Cos(alp - beta), 3)))/(CW*MW2*SW*SW2) + 
 (0.03125*Alfa2*v2*C0i(cc2, m32, m12, m22, MA02, Mh02, MZ2)*
   (MZ2*(m12 - m22 - m32 + MZ2) + (pow(MA02 - Mh02, 2)))*
   (pow(Cos(alp - beta), 3))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0234375*Alfa2*v2*C0i(cc1, m12, m22, m32, MA02, Mh02, Mh02)*
   Cos(alp - beta)*(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh02)*Cos(alp + beta))*((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta))*(pow(Csc(2*beta), 2))*
   (pow(Csc(w), 2)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*v2*C0i(cc0, m22, m32, m12, MA02, MA02, Mh02)*
   Cos(alp - beta)*(pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
      (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta), 2))*
   (pow(Csc(2*beta), 2))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*v2*C0i(cc1, m22, m32, m12, MA02, MA02, Mh02)*
   Cos(alp - beta)*(pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
      (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta), 2))*
   (pow(Csc(2*beta), 2))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*v2*C0i(cc2, m22, m32, m12, MA02, MA02, Mh02)*
   Cos(alp - beta)*(pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
      (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta), 2))*
   (pow(Csc(2*beta), 2))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0049735919716217296*Alfa*Mh02*(MA02 - Mh02 + MZ2)*B0i(bb0, 0, Mh02, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MA02*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*MHH2*(MA02 - MHH2 + MZ2)*B0i(bb0, 0, MHH2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MA02*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*(Mh02*(-Mh02 + MZ2) + MA02*(Mh02 + MZ2))*
   B0i(bb0, MA02, Mh02, MZ2)*Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/
  (CW*MA02*MW2*SW*SW2) + (0.0049735919716217296*Alfa*
   (MHH2*(-MHH2 + MZ2) + MA02*(MHH2 + MZ2))*B0i(bb0, MA02, MHH2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MA02*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MZ2*B0i(bb1, MA02, Mh02, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MZ2*B0i(bb1, MA02, MHH2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*MZ2*(-MA02 + Mh02 + MZ2)*
   C0i(cc0, m12, m22, m32, Mh02, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*MZ2*C0i(cc00, m12, m22, m32, Mh02, MZ2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.039788735772973836*Alfa*MZ2*C0i(cc00, m12, m22, m32, MHH2, MZ2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*MZ2*C0i(cc00, m32, m12, m22, MA02, MHH2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*(3*m12 + m22 - m32)*MZ2*
   C0i(cc11, m12, m22, m32, Mh02, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*(3*m12 + m22 - m32)*MZ2*
   C0i(cc11, m12, m22, m32, MHH2, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*(m12 - m22)*MZ2*C0i(cc12, m12, m22, m32, 
    Mh02, MZ2, MZ2)*Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*(m12 - m22)*MZ2*
   C0i(cc12, m12, m22, m32, MHH2, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(m12 - m22)*MZ2*C0i(cc12, m32, m12, m22, 
    MA02, MHH2, MZ2)*Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2) - (0.009947183943243459*Alfa*(m12 + 3*m22 - m32)*MZ2*
   C0i(cc22, m32, m12, m22, MA02, MHH2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*C0i(cc1, m12, m22, m32, Mh02, MZ2, MZ2)*
   Cos(alp - beta)*((MA02 - Mh02)*Mh02 + (MA02 - Mh02 + m12 - m22 + m32)*
     MZ2 + (pow(MZ2, 2)))*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2) - (0.0625*Alfa2*MZ2*(-MA02 + MHH2 + MZ2)*v2*
   C0i(cc0, m12, m22, m32, MHH2, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Csc(w), 2))*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*(-(Mh02*MHH2) + MA02*(Mh02 + MZ2) + 
    MZ2*(-MHH2 + m12 - m22 + m32 + MZ2))*v2*C0i(cc1, m12, m22, m32, 
    MHH2, MZ2, MZ2)*Cos(alp - beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*v2*C0i(cc2, m32, m12, m22, MA02, MHH2, MZ2)*
   Cos(alp - beta)*(Mh02*MHH2 - MA02*(Mh02 + MHH2) + 
    MZ2*(m12 - m22 - m32 + MZ2) + (pow(MA02, 2)))*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-46*MW2*MZ2 + 20*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(A0i(aa0, MB2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-70*MW2*MZ2 + 56*(pow(MW2, 2)) + 23*(pow(MZ2, 2)))*
   Re(A0i(aa0, MC2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-46*MW2*MZ2 + 20*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(A0i(aa0, MD2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-22*MW2*MZ2 + 12*(pow(MW2, 2)) + 9*(pow(MZ2, 2)))*
   Re(A0i(aa0, ME2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*Re(A0i(aa0, MHp2)))/
  (CW*MW2*MZ2*SW*SW2) - (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-22*MW2*MZ2 + 12*(pow(MW2, 2)) + 9*(pow(MZ2, 2)))*
   Re(A0i(aa0, ML2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-22*MW2*MZ2 + 12*(pow(MW2, 2)) + 9*(pow(MZ2, 2)))*
   Re(A0i(aa0, MM2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-46*MW2*MZ2 + 20*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(A0i(aa0, MS2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-70*MW2*MZ2 + 56*(pow(MW2, 2)) + 23*(pow(MZ2, 2)))*
   Re(A0i(aa0, MT2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-70*MW2*MZ2 + 56*(pow(MW2, 2)) + 23*(pow(MZ2, 2)))*
   Re(A0i(aa0, MU2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(4*MW2 - 3*MZ2)*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(A0i(aa0, MW2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*Re(A0i(aa0, MZ2)))/
  (CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) + 
 (0.17904931097838225*Alfa*Cos(alp - beta)*Re(B0i(bb0, 0, MW2, MW2)))/
  (CW*SW) - (0.0049735919716217296*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb0, MA02, Mh02, MZ2)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MA02, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.009947183943243459*Alfa*Cos(alp - beta)*
   Re(B0i(bb0, MA02, MHp2, MW2)))/(CW*SW*SW2) + 
 (0.0024867959858108648*Alfa*Cos(alp - beta)*(Mh02*MHH2 - 2*Mh02*MZ2 + 
    MHH2*MZ2 - MA02*(Mh02 + MHH2 + MZ2) + Cos(2*(alp - beta))*
     (MHH2*(-Mh02 + MZ2) + MA02*(Mh02 + MHH2 + MZ2) - (pow(MA02, 2))) + 
    (pow(MA02, 2)))*Re(B0i(bb0, Mh02, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.0049735919716217296*Alfa*Cos(alp - beta)*
   ((Mh02 - MHp2)*(MHH2 - MHp2) + (-2*Mh02 + MHH2 - MHp2)*MW2 + 
    ((Mh02 - MHp2)*(-MHH2 + MHp2) + (MHH2 + MHp2)*MW2)*Cos(2*(alp - beta)))*
   Re(B0i(bb0, Mh02, MHp2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*(Mh02*MHH2 + 2*MW2*(-MHH2 + 6*MW2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, MW2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.0024867959858108648*Alfa*(Mh02*MHH2 + 2*MZ2*(-MHH2 + 6*MZ2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*Cos(alp - beta)*(MHH2*(-Mh02 + MZ2) + 
    MA02*(Mh02 + MHH2 + MZ2) - (pow(MA02, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MHH2, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*((Mh02 - MHp2)*(MHH2 - MHp2) - (MHH2 + MHp2)*MW2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MHp2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*(Mh02*MHH2 + 2*MW2*(-MHH2 + 6*MW2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MW2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0024867959858108648*Alfa*(Mh02*MHH2 + 2*MZ2*(-MHH2 + 6*MZ2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.07957747154594767*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, 0, MW2)))/(CW*MZ2*SW*SW2) + 
 (0.05968310365946075*Alfa*(MT2 - MW2)*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, MB2, MT2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.05968310365946075*Alfa*MC2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, MC2, MS2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.05968310365946075*Alfa*(MU2 - MW2)*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, MD2, MU2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MW2, Mh02, MW2)))/
  (CW*MZ2*SW*(pow(SW2, 2))) - (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   (pow(Cos(alp - beta), 3))*Re(B0i(bb0, MW2, MHH2, MW2)))/
  (CW*MZ2*SW*(pow(SW2, 2))) + (0.019894367886486918*Alfa*Cos(alp - beta)*
   (6*MZ2*(pow(MW2, 2)) + 4*(pow(MW2, 3)) - 6*MW2*(pow(MZ2, 2)) + 
    (pow(MZ2, 3)))*Re(B0i(bb0, MW2, MW2, MZ2)))/
  (CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MB2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MB2, MB2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MC2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MC2, MC2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MD2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MD2, MD2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*ME2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, ME2, ME2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MZ2, Mh02, MZ2)))/
  (CW*MW2*SW*(pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   (pow(Cos(alp - beta), 3))*Re(B0i(bb0, MZ2, MHH2, MZ2)))/
  (CW*MW2*SW*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*ML2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, ML2, ML2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*MM2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MM2, MM2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MS2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MS2, MS2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MT2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MT2, MT2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MU2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MU2, MU2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-8*MW2*MZ2 + 13*(pow(MW2, 2)) + 2*(pow(MZ2, 2)))*
   Re(B0i(bb0, MZ2, MW2, MW2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.026525823848649224*Alfa*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MB2, MB2)))/(CW*MW2*MZ2*SW) - 
 (0.05305164769729845*Alfa*(8*MW2 - 5*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MC2, MC2)))/(CW*MW2*MZ2*SW) - 
 (0.026525823848649224*Alfa*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MD2, MD2)))/(CW*MW2*MZ2*SW) - 
 (0.07957747154594767*Alfa*(4*MW2 - 3*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, ME2, ME2)))/(CW*MW2*MZ2*SW) + 
 (0.07957747154594767*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MHp2, MHp2)))/(CW*MW2*MZ2*SW) - 
 (0.07957747154594767*Alfa*(4*MW2 - 3*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, ML2, ML2)))/(CW*MW2*MZ2*SW) - 
 (0.07957747154594767*Alfa*(4*MW2 - 3*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MM2, MM2)))/(CW*MW2*MZ2*SW) - 
 (0.026525823848649224*Alfa*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MS2, MS2)))/(CW*MW2*MZ2*SW) - 
 (0.05305164769729845*Alfa*(8*MW2 - 5*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MT2, MT2)))/(CW*MW2*MZ2*SW) - 
 (0.05305164769729845*Alfa*(8*MW2 - 5*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MU2, MU2)))/(CW*MW2*MZ2*SW) + 
 (0.07957747154594767*Alfa*(6*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MW2, MW2)))/(CW*MW2*MZ2*SW) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, ME2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, ML2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, MM2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.15915494309189535*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, MW2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MA02, MHp2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.1193662073189215*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MB2, MT2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.1193662073189215*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MC2, MS2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.1193662073189215*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MD2, MU2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MW2, Mh02, MHp2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MW2, Mh02, MW2)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MW2, MHH2, MHp2)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MW2, MHH2, MW2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(8*MW2 + MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MW2, MZ2)))/(CW*MW2*SW*(pow(MZ2, 2))*
   (pow(SW2, 2))) + (0.05968310365946075*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*Re(B0i(bb00, MZ2, 0, 0)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MZ2, MA02, Mh02)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MZ2, MA02, MHH2)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.006631455962162306*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MB2, MB2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.006631455962162306*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 
    17*(pow(MZ2, 2)))*Re(B0i(bb00, MZ2, MC2, MC2)))/
  (CW*MW2*SW*(pow(MZ2, 3))*(pow(SW2, 2))) + 
 (0.006631455962162306*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MD2, MD2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, ME2, ME2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) - (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb00, MZ2, Mh02, MZ2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MZ2, MHH2, MZ2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*(pow(2*MW2 - MZ2, 3))*
   Re(B0i(bb00, MZ2, MHp2, MHp2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, ML2, ML2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MM2, MM2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.006631455962162306*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MS2, MS2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.006631455962162306*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 
    17*(pow(MZ2, 2)))*Re(B0i(bb00, MZ2, MT2, MT2)))/
  (CW*MW2*SW*(pow(MZ2, 3))*(pow(SW2, 2))) + 
 (0.006631455962162306*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MU2, MU2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) - (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-4*MW2*MZ2 + 12*(pow(MW2, 2)) + (pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MW2, MW2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.026525823848649224*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MB2, MB2)))/(CW*SW) + 
 (0.1061032953945969*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MC2, MC2)))/
  (CW*SW) + (0.026525823848649224*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MD2, MD2)))/(CW*SW) + 
 (0.07957747154594767*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, ME2, ME2)))/
  (CW*SW) + (0.07957747154594767*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, ML2, ML2)))/(CW*SW) + 
 (0.07957747154594767*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MM2, MM2)))/
  (CW*SW) + (0.026525823848649224*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MS2, MS2)))/(CW*SW) + 
 (0.1061032953945969*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MT2, MT2)))/
  (CW*SW) + (0.1061032953945969*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MU2, MU2)))/(CW*SW) + 
 (0.039788735772973836*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MW2, MW2)))/
  (CW*SW) - (0.05968310365946075*Alfa*MB2*Cos(alp - beta)*
   (pow(Cot(beta), 2))*Re(B0i(bb1, MA02, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MC2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MD2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*ME2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(bb1, MA02, ME2, ME2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb1, MA02, Mh02, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MZ2*Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb1, MA02, MHH2, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*Re(B0i(bb1, MA02, MHp2, MW2)))/
  (CW*SW*SW2) - (0.019894367886486918*Alfa*ML2*Cos(alp - beta)*
   (pow(Tan(beta), 2))*Re(B0i(bb1, MA02, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MM2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(bb1, MA02, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MS2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MT2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MU2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MU2, MU2)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*MZ2*Cos(alp - beta)*
   (-2*Mh02 + MHH2 + MHH2*Cos(2*(alp - beta)))*Re(B0i(bb1, Mh02, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-2*Mh02 + MHH2 + MHH2*Cos(2*(alp - beta)))*Re(B0i(bb1, Mh02, MHp2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) - (0.019894367886486918*Alfa*MHH2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, Mh02, MW2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) - (0.009947183943243459*Alfa*MHH2*MZ2*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb1, Mh02, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MHH2*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, MHH2, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MHH2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, MHH2, MHp2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) + (0.019894367886486918*Alfa*MHH2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, MHH2, MW2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) + (0.009947183943243459*Alfa*MHH2*MZ2*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb1, MHH2, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, ME2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, ML2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, MM2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, MW2)))/(CW*MZ2*SW*SW2) - 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MB2, MT2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MC2, MS2)))/(CW*MZ2*SW*(pow(SW2, 2))) - 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MD2, MU2)))/(CW*MZ2*SW*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*MW2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MW2, MZ2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*Cos(alp - beta)*Re(B0i(bb1, MZ2, 0, 0)))/
  (CW*SW*(pow(SW2, 2))) - (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MC2, MC2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, ME2, ME2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, ML2, ML2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MM2, MM2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MT2, MT2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MU2, MU2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*Cos(alp - beta)*(pow(MW2, 2))*
   Re(B0i(bb1, MZ2, MW2, MW2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) + 
 (0.0012433979929054324*Alfa*Cos(alp - beta)*
   (pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*Mh02)*
       Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, MA02, MA02, Mh02)))/(CW*MW2*SW*SW2) + 
 (0.0012433979929054324*Alfa*Cos(alp - beta)*(pow(Csc(2*beta), 2))*
   (pow((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
       Sin(alp + beta), 2))*Re(B0i(dbb0, MA02, MA02, MHH2)))/
  (CW*MW2*SW*SW2) - (0.0049735919716217296*Alfa*
   ((MA02 + Mh02)*MZ2 - (pow(MA02 - Mh02, 2)))*
   (pow(Cos(alp - beta), 3))*Re(B0i(dbb0, MA02, Mh02, MZ2)))/
  (CW*MW2*SW*SW2) - (0.0049735919716217296*Alfa*Cos(alp - beta)*
   ((MA02 + MHH2)*MZ2 - (pow(MA02 - MHH2, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, MA02, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) + (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-((MA02 + MHp2)*MW2) + (pow(MA02 - MHp2, 2)))*
   Re(B0i(dbb0, MA02, MHp2, MW2)))/(CW*MW2*SW*SW2) + 
 (0.0006216989964527162*Alfa*Cos(alp - beta)*
   (pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*Mh02)*
       Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh02, MA02, MA02)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*((MA02 + Mh02)*MZ2 - (pow(MA02 - Mh02, 2)))*
   (pow(Cos(alp - beta), 3))*Re(B0i(dbb0, Mh02, MA02, MZ2)))/
  (CW*MW2*SW*SW2) + (0.005595290968074445*Alfa*Cos(alp - beta)*
   (pow(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
      (2*M2 - 3*Mh02)*Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh02, Mh02, Mh02)))/(CW*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*(pow(Cos(alp - beta), 3))*
   (pow(M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*Sin(2*alp), 2))*
   Re(B0i(dbb0, Mh02, Mh02, MHH2)))/(CW*MW2*SW*SW2) + 
 (0.0024867959858108648*Alfa*Cos(alp - beta)*(pow(Csc(2*beta), 2))*
   (pow(Sin(alp - beta), 2))*
   (pow((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta), 2))*
   Re(B0i(dbb0, Mh02, MHH2, MHH2)))/(CW*MW2*SW*SW2) + 
 (0.0012433979929054324*Alfa*Cos(alp - beta)*
   (pow((Mh02 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh02 + 2*MHp2)*
       Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh02, MHp2, MHp2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*(-((Mh02 + MHp2)*MW2) + 
    (pow(Mh02 - MHp2, 2)))*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb0, Mh02, MHp2, MW2)))/(CW*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*Cos(alp - beta)*
   (-2*Mh02*MW2 + (pow(Mh02, 2)) + 12*(pow(MW2, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, Mh02, MW2, MW2)))/
  (CW*MW2*SW*SW2) + (0.0024867959858108648*Alfa*Cos(alp - beta)*
   (-2*Mh02*MZ2 + (pow(Mh02, 2)) + 12*(pow(MZ2, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, Mh02, MZ2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.029841551829730376*Alfa*MB2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MC2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MD2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*ME2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, ME2, ME2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*(pow(MZ2, 2))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, MZ2, Mh02, MZ2)))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*(pow(MZ2, 2))*
   (pow(Cos(alp - beta), 3))*Re(B0i(dbb0, MZ2, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.009947183943243459*Alfa*ML2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MM2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MS2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MT2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MU2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MU2, MU2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*(9*MW2 - 2*MZ2)*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MW2, MW2)))/(CW*SW*SW2) - 
 (0.05305164769729845*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MB2, MB2)))/
  (CW*SW) - (0.2122065907891938*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, MC2, MC2)))/(CW*SW) - 
 (0.05305164769729845*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MD2, MD2)))/
  (CW*SW) - (0.15915494309189535*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, ME2, ME2)))/(CW*SW) + 
 (0.07957747154594767*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MHp2, MHp2)))/
  (CW*SW) - (0.15915494309189535*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, ML2, ML2)))/(CW*SW) - 
 (0.15915494309189535*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MM2, MM2)))/
  (CW*SW) - (0.05305164769729845*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, MS2, MS2)))/(CW*SW) - 
 (0.2122065907891938*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MT2, MT2)))/
  (CW*SW) - (0.2122065907891938*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, MU2, MU2)))/(CW*SW) + 
 (0.238732414637843*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MW2, MW2)))/
  (CW*SW) + (0.05968310365946075*Alfa*MZ2*Cos(alp - beta)*
   Re(B0i(dbb00, MZ2, 0, 0)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb00, MZ2, MA02, Mh02)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(dbb00, MZ2, MA02, MHH2)))/(CW*MW2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MB2, MB2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MC2, MC2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MD2, MD2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, ME2, ME2)))/(CW*MW2*MZ2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(dbb00, MZ2, Mh02, MZ2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb00, MZ2, MHH2, MZ2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*(pow(-2*MW2 + MZ2, 2))*
   Re(B0i(dbb00, MZ2, MHp2, MHp2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, ML2, ML2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MM2, MM2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MS2, MS2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MT2, MT2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MU2, MU2)))/(CW*MW2*MZ2*SW*SW2) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 12*(pow(MW2, 2)) + (pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MW2, MW2)))/(CW*MW2*MZ2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MB2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MC2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MD2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MA02*ME2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(dbb1, MA02, ME2, ME2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MA02*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb1, MA02, Mh02, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MA02*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb1, MA02, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*MA02*Cos(alp - beta)*
   Re(B0i(dbb1, MA02, MHp2, MW2)))/(CW*SW*SW2) - 
 (0.019894367886486918*Alfa*MA02*ML2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(dbb1, MA02, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MA02*MM2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(dbb1, MA02, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MS2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MT2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MU2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MU2, MU2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*Mh02*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb1, Mh02, MA02, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*Mh02*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb1, Mh02, MHp2, MW2)))/(CW*SW*SW2) + 
 (0.019894367886486918*Alfa*Mh02*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb1, Mh02, MW2, MW2)))/(CW*SW*SW2) + 
 (0.009947183943243459*Alfa*Mh02*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb1, Mh02, MZ2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.029841551829730376*Alfa*Cos(alp - beta)*
   (pow(MZ2, 2))*Re(B0i(dbb1, MZ2, 0, 0)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, ME2, ME2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MU2, MU2)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*MW2*Cos(alp - beta)*Re(B0i(dbb1, MZ2, MW2, MW2)))/
  (CW*SW*SW2) + (0.05968310365946075*Alfa*MB2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MB2, MB2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MC2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MC2, MC2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MD2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MD2, MD2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MS2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MS2, MS2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MT2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MT2, MT2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MU2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MU2, MU2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*(pow(MB2, 2))*
   Re(B0i(bb0, Mh02, MB2, MB2))*Sin(alp)*(Cos(alp) - Cot(beta)*Sin(alp)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*
   (pow(MC2, 2))*Re(B0i(bb0, Mh02, MC2, MC2))*Sin(alp)*
   (Cos(alp) - Cot(beta)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*(pow(MD2, 2))*
   Re(B0i(bb0, Mh02, MD2, MD2))*Sin(alp)*(Cos(alp) - Cot(beta)*Sin(alp)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*
   (pow(MS2, 2))*Re(B0i(bb0, Mh02, MS2, MS2))*Sin(alp)*
   (Cos(alp) - Cot(beta)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*(pow(MT2, 2))*
   Re(B0i(bb0, Mh02, MT2, MT2))*Sin(alp)*(Cos(alp) - Cot(beta)*Sin(alp)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*
   (pow(MU2, 2))*Re(B0i(bb0, Mh02, MU2, MU2))*Sin(alp)*
   (Cos(alp) - Cot(beta)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.03125*Alfa2*(-MA02 + Mh02 + MZ2)*v2*C0i(cc0, m22, m32, m12, Mh02, 
    MHH2, MZ2)*(pow(Cos(alp - beta), 3))*(pow(Csc(w), 2))*
   (M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*Sin(2*alp)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.03125*Alfa2*(MA02 - Mh02 + MZ2)*v2*C0i(cc1, m22, m32, m12, Mh02, MHH2, 
    MZ2)*(pow(Cos(alp - beta), 3))*(pow(Csc(w), 2))*
   (M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*Sin(2*alp)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.03125*Alfa2*(MA02 - Mh02 + MZ2)*v2*C0i(cc2, m22, m32, m12, Mh02, MHH2, 
    MZ2)*(pow(Cos(alp - beta), 3))*(pow(Csc(w), 2))*
   (M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*Sin(2*alp)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.05968310365946075*Alfa*MB2*B0i(bb1, MA02, MB2, MB2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MC2*B0i(bb1, MA02, MC2, MC2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MD2*B0i(bb1, MA02, MD2, MD2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MS2*B0i(bb1, MA02, MS2, MS2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MT2*B0i(bb1, MA02, MT2, MT2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MU2*B0i(bb1, MA02, MU2, MU2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*(pow(MB2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MB2, MB2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.05968310365946075*Alfa*(pow(MC2, 2))*
   (pow(Csc(beta), 2))*Re(B0i(bb0, MHH2, MC2, MC2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*(pow(MD2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MD2, MD2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*(pow(ME2, 2))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, ME2, ME2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*(pow(ML2, 2))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, ML2, ML2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*(pow(MM2, 2))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, MM2, MM2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.05968310365946075*Alfa*(pow(MS2, 2))*
   (pow(Csc(beta), 2))*Re(B0i(bb0, MHH2, MS2, MS2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*(pow(MT2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MT2, MT2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.05968310365946075*Alfa*(pow(MU2, 2))*
   (pow(Csc(beta), 2))*Re(B0i(bb0, MHH2, MU2, MU2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MB2*MHH2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MB2, MB2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MC2*MHH2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MC2, MC2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MD2*MHH2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MD2, MD2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*ME2*MHH2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, MHH2, ME2, ME2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MHH2*ML2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, MHH2, ML2, ML2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MHH2*MM2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, MHH2, MM2, MM2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MHH2*MS2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MS2, MS2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MHH2*MT2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MT2, MT2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MHH2*MU2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MU2, MU2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.0012433979929054324*Alfa*(MA02 - Mh02)*
   B0i(bb0, 0, MA02, Mh02)*((2*MA02 - Mh02)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh02)*Cos(alp + beta))*Csc(2*beta)*
   Sin(2*(alp - beta)))/(CW*MA02*MW2*SW*SW2) - 
 (0.0012433979929054324*Alfa*(MA02 - Mh02)*B0i(bb0, MA02, MA02, Mh02)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(2*beta)*Sin(2*(alp - beta)))/(CW*MA02*MW2*SW*SW2) + 
 (0.003730193978716297*Alfa*(-MA02 + Mh02 + MZ2)*C0i(cc0, m22, m32, m12, 
    Mh02, Mh02, MZ2)*(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh02)*Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) - (0.003730193978716297*Alfa*(MA02 - Mh02 + MZ2)*
   C0i(cc1, m22, m32, m12, Mh02, Mh02, MZ2)*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh02)*Cos(3*alp - beta) + (2*M2 - 3*Mh02)*Cos(alp + beta))*
   Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/(CW*MW2*SW*SW2) - 
 (0.003730193978716297*Alfa*(MA02 - Mh02 + MZ2)*C0i(cc2, m22, m32, m12, 
    Mh02, Mh02, MZ2)*(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh02)*Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) + (0.0012433979929054324*Alfa*(MA02 - Mh02 + MZ2)*
   C0i(cc0, m12, m32, m22, MA02, Mh02, MZ2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) - (0.0012433979929054324*Alfa*(-MA02 + Mh02 + MZ2)*
   C0i(cc1, m12, m32, m22, MA02, Mh02, MZ2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) - (0.0012433979929054324*Alfa*(-MA02 + Mh02 + MZ2)*
   C0i(cc2, m12, m32, m22, MA02, Mh02, MZ2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*ME2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, Mh02, ME2, ME2))*Sin(alp)*(MHH2*Cos(alp - beta)*Sin(alp) - 
    Mh02*Sin(beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*ML2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, Mh02, ML2, ML2))*Sin(alp)*(MHH2*Cos(alp - beta)*Sin(alp) - 
    Mh02*Sin(beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*MM2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, Mh02, MM2, MM2))*Sin(alp)*(MHH2*Cos(alp - beta)*Sin(alp) - 
    Mh02*Sin(beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.03125*Alfa2*(-MA02 + MHH2 + MZ2)*v2*C0i(cc0, m22, m32, m12, MHH2, 
    MHH2, MZ2)*Cos(alp - beta)*Csc(2*beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2))*((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*(MA02 - MHH2 + MZ2)*v2*C0i(cc1, m22, m32, m12, MHH2, MHH2, 
    MZ2)*Cos(alp - beta)*Csc(2*beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2))*((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*(MA02 - MHH2 + MZ2)*v2*C0i(cc2, m22, m32, m12, MHH2, MHH2, 
    MZ2)*Cos(alp - beta)*Csc(2*beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2))*((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.009947183943243459*Alfa*MZ2*C0i(cc0, m22, m12, m32, Mh02, MHH2, MZ2)*
   Cos(alp - beta)*Csc(beta)*(pow(Sin(alp - beta), 2))*Sec(beta)*
   ((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*(MA02 - MHH2 + MZ2)*C0i(cc1, m22, m12, m32, 
    Mh02, MHH2, MZ2)*Cos(alp - beta)*Csc(beta)*(pow(Sin(alp - beta), 2))*
   Sec(beta)*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*MW2*SW*SW2) + (0.0078125*Alfa2*v2*C0i(cc1, m12, m22, m32, MA02, 
    Mh02, MHH2)*((2*MA02 - Mh02)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh02)*Cos(alp + beta))*(pow(Csc(2*beta), 2))*
   (pow(Csc(w), 2))*Sin(2*(alp - beta))*
   ((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0018650969893581485*Alfa*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh02)*Cos(3*alp - beta) + (2*M2 - 3*Mh02)*Cos(alp + beta))*
   (pow(Csc(2*beta), 2))*Re(B0i(bb0, Mh02, Mh02, Mh02))*
   Sin(2*(alp - beta))*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0018650969893581485*Alfa*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh02)*Cos(3*alp - beta) + (2*M2 - 3*Mh02)*Cos(alp + beta))*
   (pow(Csc(2*beta), 2))*Re(B0i(bb0, MHH2, Mh02, Mh02))*
   Sin(2*(alp - beta))*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - (0.0024867959858108648*Alfa*Cos(alp - beta)*
   (pow(Csc(2*beta), 2))*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, Mh02, MHH2))*(9*M2*(Mh02 + MHH2) - 
    (2*Mh02 + MHH2)*(Mh02 + 2*MHH2) + (3*M2 - Mh02 - 2*MHH2)*
     (3*M2 - 2*Mh02 - MHH2)*Cos(4*alp) - 8*(pow(M2, 2)) + 
    M2*(-(M2*Cos(4*beta)) + 2*(Mh02 - MHH2)*Sin(2*alp)*Sin(2*beta))))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.0024867959858108648*Alfa*Cos(alp - beta)*
   (pow(Csc(2*beta), 2))*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, Mh02, MHH2))*(9*M2*(Mh02 + MHH2) - 
    (2*Mh02 + MHH2)*(Mh02 + 2*MHH2) + (3*M2 - Mh02 - 2*MHH2)*
     (3*M2 - 2*Mh02 - MHH2)*Cos(4*alp) - 8*(pow(M2, 2)) + 
    M2*(-(M2*Cos(4*beta)) + 2*(Mh02 - MHH2)*Sin(2*alp)*Sin(2*beta))))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.003730193978716297*Alfa*Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, MHH2, MHH2))*(M2 + (3*M2 - Mh02 - 2*MHH2)*Csc(2*beta)*
     Sin(2*alp))*(M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.003730193978716297*Alfa*Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MHH2, MHH2))*(M2 + (3*M2 - Mh02 - 2*MHH2)*Csc(2*beta)*
     Sin(2*alp))*(M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.015625*Alfa2*v2*C0i(cc1, m12, m22, m32, MA02, MHH2, MHH2)*Csc(2*beta)*
   (pow(Csc(w), 2))*(pow(Sin(alp - beta), 2))*
   (M2 + (3*M2 - Mh02 - 2*MHH2)*Csc(2*beta)*Sin(2*alp))*
   ((2*MA02 - MHH2)*Sin(alp - 3*beta) + (4*M2 - 2*MA02 - 3*MHH2)*
     Sin(alp + beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.015625*Alfa2*v2*C0i(cc2, m32, m22, m12, MA02, Mh02, MHH2)*
   (pow(Cos(alp - beta), 2))*(pow(Csc(2*beta), 2))*
   (pow(Csc(w), 2))*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta))*
   ((2*MA02 - MHH2)*Sin(alp - 3*beta) + (4*M2 - 2*MA02 - 3*MHH2)*
     Sin(alp + beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0024867959858108648*Alfa*(MA02 - MHH2)*B0i(bb0, 0, MA02, MHH2)*
   Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*MA02*MW2*SW*SW2) + 
 (0.0024867959858108648*Alfa*(MA02 - MHH2)*B0i(bb0, MA02, MA02, MHH2)*
   Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*MA02*MW2*SW*SW2) + 
 (0.0078125*Alfa2*(MA02 - Mh02 + MZ2)*v2*C0i(cc0, m12, m32, m22, MA02, 
    MHH2, MZ2)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Sec(beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*(-MA02 + Mh02 + MZ2)*v2*C0i(cc1, m12, m32, m22, MA02, 
    MHH2, MZ2)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Sec(beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*(-MA02 + Mh02 + MZ2)*v2*C0i(cc2, m12, m32, m22, MA02, 
    MHH2, MZ2)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Sec(beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0078125*Alfa2*v2*C0i(cc0, m22, m32, m12, MA02, MA02, MHH2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*(pow(Csc(2*beta), 2))*(pow(Csc(w), 2))*
   Sin(alp - beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0078125*Alfa2*v2*C0i(cc1, m22, m32, m12, MA02, MA02, MHH2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*(pow(Csc(2*beta), 2))*(pow(Csc(w), 2))*
   Sin(alp - beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0078125*Alfa2*v2*C0i(cc2, m22, m32, m12, MA02, MA02, MHH2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*(pow(Csc(2*beta), 2))*(pow(Csc(w), 2))*
   Sin(alp - beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.00015542474911317905*Alfa*((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh02, MA02, MA02))*Sin(alp - beta)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.00015542474911317905*Alfa*((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MA02, MA02))*Sin(alp - beta)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0003108494982263581*Alfa*((Mh02 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh02 + 2*MHp2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh02, MHp2, MHp2))*Sin(alp - beta)*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.0003108494982263581*Alfa*((Mh02 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh02 + 2*MHp2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MHp2, MHp2))*Sin(alp - beta)*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MB2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MB2, MB2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MC2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MC2, MC2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MD2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MD2, MD2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MS2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MS2, MS2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MT2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MT2, MT2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MU2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MU2, MU2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MB2*Mh02*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MB2, MB2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MC2*Mh02*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MC2, MC2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MD2*Mh02*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MD2, MD2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*Mh02*MS2*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MS2, MS2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*Mh02*MT2*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MT2, MT2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*Mh02*MU2*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MU2, MU2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) + 
 (1.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, ME2, ME2, ME2)*(pow(ME2, 2))*
   Sec(beta)*Sin(alp)*Tan(beta))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, ML2, ML2, ML2)*(pow(ML2, 2))*
   Sec(beta)*Sin(alp)*Tan(beta))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MM2, MM2, MM2)*(pow(MM2, 2))*
   Sec(beta)*Sin(alp)*Tan(beta))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ME2*MZ*MZ2*B0i(bb0, m32, ME2, ME2)*Sec(beta)*Sin(alp)*Tan(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ML2*MZ*MZ2*B0i(bb0, m32, ML2, ML2)*Sec(beta)*Sin(alp)*Tan(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MM2*MZ*MZ2*B0i(bb0, m32, MM2, MM2)*Sec(beta)*Sin(alp)*Tan(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ME2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, ME2, ME2, ME2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ML2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, ML2, ML2, ML2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MM2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MM2, MM2, MM2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ME2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, ME2, ME2, ME2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ML2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, ML2, ML2, ML2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MM2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MM2, MM2, MM2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) - 
 (0.019894367886486918*Alfa*ME2*B0i(bb1, MA02, ME2, ME2)*Sin(alp - beta)*
   Tan(beta))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*ML2*B0i(bb1, MA02, ML2, ML2)*Sin(alp - beta)*
   Tan(beta))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MM2*B0i(bb1, MA02, MM2, MM2)*Sin(alp - beta)*
   Tan(beta))/(CW*MW2*SW*SW2) - (0.039788735772973836*Alfa*(pow(ME2, 2))*
   (pow(Sin(alp), 3))*Re(B0i(dbb0, Mh02, ME2, ME2))*Sec(beta)*
   (Cot(alp) + Tan(beta)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*(pow(ML2, 2))*(pow(Sin(alp), 3))*
   Re(B0i(dbb0, Mh02, ML2, ML2))*Sec(beta)*(Cot(alp) + Tan(beta)))/
  (CW*MW2*SW*SW2) - (0.039788735772973836*Alfa*(pow(MM2, 2))*
   (pow(Sin(alp), 3))*Re(B0i(dbb0, Mh02, MM2, MM2))*Sec(beta)*
   (Cot(alp) + Tan(beta)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*ME2*Mh02*(pow(Sin(alp), 3))*
   Re(B0i(dbb1, Mh02, ME2, ME2))*Sec(beta)*(Cot(alp) + Tan(beta)))/
  (CW*MW2*SW*SW2) - (0.019894367886486918*Alfa*Mh02*ML2*
   (pow(Sin(alp), 3))*Re(B0i(dbb1, Mh02, ML2, ML2))*Sec(beta)*
   (Cot(alp) + Tan(beta)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*Mh02*MM2*(pow(Sin(alp), 3))*
   Re(B0i(dbb1, Mh02, MM2, MM2))*Sec(beta)*(Cot(alp) + Tan(beta)))/
  (CW*MW2*SW*SW2) + (1.*Re(B0i(bb0, Mh02, ME2, ME2))*
   (Alfa*Cos(alp)*(pow(ME2, 2))*(pow(Sin(alp), 2))*Sec(beta) - 
    Alfa*(pow(ME2, 2))*(pow(Cos(alp), 2))*Sec(beta)*Sin(alp)*
     Tan(beta)))/(8*CW*Mh02*MW2*Pi*SW*SW2 - 8*CW*MHH2*MW2*Pi*SW*SW2) + 
 (1.*Re(B0i(bb0, Mh02, ML2, ML2))*(Alfa*Cos(alp)*(pow(ML2, 2))*
     (pow(Sin(alp), 2))*Sec(beta) - Alfa*(pow(ML2, 2))*
     (pow(Cos(alp), 2))*Sec(beta)*Sin(alp)*Tan(beta)))/
  (8*CW*Mh02*MW2*Pi*SW*SW2 - 8*CW*MHH2*MW2*Pi*SW*SW2) + 
 (1.*Re(B0i(bb0, Mh02, MM2, MM2))*(Alfa*Cos(alp)*(pow(MM2, 2))*
     (pow(Sin(alp), 2))*Sec(beta) - Alfa*(pow(MM2, 2))*
     (pow(Cos(alp), 2))*Sec(beta)*Sin(alp)*Tan(beta)))/
  (8*CW*Mh02*MW2*Pi*SW*SW2 - 8*CW*MHH2*MW2*Pi*SW*SW2)
;
}

} //end namespace TypeLS

namespace TypeFL{
ComplexType CAhZ(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32)
{
return 0. - (0.5*Cos(alp - beta))/(CW*SW) + 
 (0.009947183943243459*Alfa*MZ2*B0i(bb0, m12, Mh02, MZ2)*Cos(alp - beta))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*B0i(bb0, m12, MHp2, MW2)*
   Cos(alp - beta))/(CW*SW) + (0.009947183943243459*Alfa*MZ2*
   B0i(bb0, m22, MA02, MZ2)*Cos(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*B0i(bb0, m22, MHp2, MW2)*Cos(alp - beta))/
  (CW*SW) - (0.07957747154594767*Alfa*CW*B0i(bb0, m32, MW2, MW2)*
   Cos(alp - beta))/(SW*SW2) - (0.009947183943243459*Alfa*MZ2*
   B0i(bb1, m12, Mh02, MZ2)*Cos(alp - beta))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*B0i(bb1, m12, MHp2, MW2)*Cos(alp - beta))/
  (CW*SW) - (0.009947183943243459*Alfa*MZ2*B0i(bb1, m22, MA02, MZ2)*
   Cos(alp - beta))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*B0i(bb1, m22, MHp2, MW2)*Cos(alp - beta))/
  (CW*SW) + (0.07957747154594767*Alfa*CW*C0i(cc00, m22, m32, m12, MHp2, 
    MW2, MW2)*Cos(alp - beta))/(SW*SW2) + 
 (0.07957747154594767*Alfa*(2*MW2 - MZ2)*C0i(cc00, m32, m12, m22, MHp2, 
    MHp2, MW2)*Cos(alp - beta))/(CW*MZ2*SW*SW2) - 
 (0.019894367886486918*Alfa*(MW2*(2*(MA02 - MHp2)*(-Mh02 + MHp2) - 
      (MA02 + Mh02 - 2*MHp2 - m12 + m22 + 2*m32)*MW2) + 
    ((MA02 - MHp2)*(Mh02 - MHp2) + (MA02 + Mh02 - 2*MHp2)*MW2)*MZ2)*
   C0i(cc1, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/
  (MW*MW2*MZ*SW*SW2) + (0.019894367886486918*Alfa*CW*(m12 + 3*m22 - m32)*
   C0i(cc11, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/(SW*SW2) + 
 (0.039788735772973836*Alfa*CW*(2*(m12 + m22) - m32)*
   C0i(cc12, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/(SW*SW2) - 
 (0.039788735772973836*Alfa*(m12 - m22)*(2*MW2 - MZ2)*
   C0i(cc12, m32, m12, m22, MHp2, MHp2, MW2)*Cos(alp - beta))/
  (CW*MZ2*SW*SW2) - (0.019894367886486918*Alfa*
   (MW2*(2*(MA02 - MHp2)*(-Mh02 + MHp2) - (MA02 + Mh02 - 2*MHp2 + m12 - 
        m22 + 2*m32)*MW2) + ((MA02 - MHp2)*(Mh02 - MHp2) + 
      (MA02 + Mh02 - 2*MHp2)*MW2)*MZ2)*C0i(cc2, m22, m32, m12, MHp2, MW2, 
    MW2)*Cos(alp - beta))/(MW*MW2*MZ*SW*SW2) - 
 (0.019894367886486918*Alfa*((MA02 - MHp2)*(-Mh02 + MHp2) + 
    MW2*(m12 - m22 - m32 + MW2))*(2*MW2 - MZ2)*
   C0i(cc2, m32, m12, m22, MHp2, MHp2, MW2)*Cos(alp - beta))/
  (CW*MW2*MZ2*SW*SW2) + (0.019894367886486918*Alfa*CW*(3*m12 + m22 - m32)*
   C0i(cc22, m22, m32, m12, MHp2, MW2, MW2)*Cos(alp - beta))/(SW*SW2) + 
 (0.019894367886486918*Alfa*(m12 + 3*m22 - m32)*(2*MW2 - MZ2)*
   C0i(cc22, m32, m12, m22, MHp2, MHp2, MW2)*Cos(alp - beta))/
  (CW*MZ2*SW*SW2) + (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MC2, MC2, MC2)*
   Cos(alp)*Cot(beta)*Csc(beta)*(pow(MC2, 2)))/
  (-8*MW*MW2*MZ2*Pi*SW + 8*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, ME2, ME2, ME2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(ME2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, ML2, ML2, ML2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(ML2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MM2, MM2, MM2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MM2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MT2, MT2, MT2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MT2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MU2, MU2, MU2)*Cos(alp)*Cot(beta)*
   Csc(beta)*(pow(MU2, 2)))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MC2*MZ*MZ2*B0i(bb0, m32, MC2, MC2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ME2*MZ*MZ2*B0i(bb0, m32, ME2, ME2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ML2*MZ*MZ2*B0i(bb0, m32, ML2, ML2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MM2*MZ*MZ2*B0i(bb0, m32, MM2, MM2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MT2*MZ*MZ2*B0i(bb0, m32, MT2, MT2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MU2*MZ*MZ2*B0i(bb0, m32, MU2, MU2)*Cos(alp)*Cot(beta)*Csc(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MC2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MC2, MC2, MC2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ME2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, ME2, ME2, ME2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ML2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, ML2, ML2, ML2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MM2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MM2, MM2, MM2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MT2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MT2, MT2, MT2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MU2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MU2, MU2, MU2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MC2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MC2, MC2, MC2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ME2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, ME2, ME2, ME2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*ML2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, ML2, ML2, ML2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (1.*Alfa*MM2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MM2, MM2, MM2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MT2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MT2, MT2, MT2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*m12*MU2*MZ*MZ2*C0i(cc2, m22, m32, m12, MU2, MU2, MU2)*Cos(alp)*
   Cot(beta)*Csc(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (0.019894367886486918*Alfa*MZ2*C0i(cc0, m22, m32, m12, MHp2, MW2, MW2)*
   Cos(alp - beta)*(Mh02 - MHp2 - ((4*MHp2 + m32)*MW2 + 
      (MA02 - MHp2)*(MW2*(2*MHp2 + MW2) - (MHp2 + MW2)*MZ2 + 
        Mh02*(-2*MW2 + MZ2))*(pow(MW2, -1)))*(pow(-MW2 + MZ2, -1))))/
  (MW*MZ*SW) - (0.039788735772973836*Alfa*MZ2*C0i(cc00, m32, m12, m22, 
    MA02, Mh02, MZ2)*(pow(Cos(alp - beta), 3)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(m12 - m22)*MZ2*C0i(cc12, m32, m12, m22, 
    MA02, Mh02, MZ2)*(pow(Cos(alp - beta), 3)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*(m12 + 3*m22 - m32)*MZ2*
   C0i(cc22, m32, m12, m22, MA02, Mh02, MZ2)*
   (pow(Cos(alp - beta), 3)))/(CW*MW2*SW*SW2) + 
 (0.03125*Alfa2*v2*C0i(cc2, m32, m12, m22, MA02, Mh02, MZ2)*
   (MZ2*(m12 - m22 - m32 + MZ2) + (pow(MA02 - Mh02, 2)))*
   (pow(Cos(alp - beta), 3))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0234375*Alfa2*v2*C0i(cc1, m12, m22, m32, MA02, Mh02, Mh02)*
   Cos(alp - beta)*(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh02)*Cos(alp + beta))*((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta))*(pow(Csc(2*beta), 2))*
   (pow(Csc(w), 2)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*v2*C0i(cc0, m22, m32, m12, MA02, MA02, Mh02)*
   Cos(alp - beta)*(pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
      (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta), 2))*
   (pow(Csc(2*beta), 2))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*v2*C0i(cc1, m22, m32, m12, MA02, MA02, Mh02)*
   Cos(alp - beta)*(pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
      (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta), 2))*
   (pow(Csc(2*beta), 2))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*v2*C0i(cc2, m22, m32, m12, MA02, MA02, Mh02)*
   Cos(alp - beta)*(pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
      (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta), 2))*
   (pow(Csc(2*beta), 2))*(pow(Csc(w), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0049735919716217296*Alfa*Mh02*(MA02 - Mh02 + MZ2)*B0i(bb0, 0, Mh02, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MA02*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*MHH2*(MA02 - MHH2 + MZ2)*B0i(bb0, 0, MHH2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MA02*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*(Mh02*(-Mh02 + MZ2) + MA02*(Mh02 + MZ2))*
   B0i(bb0, MA02, Mh02, MZ2)*Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/
  (CW*MA02*MW2*SW*SW2) + (0.0049735919716217296*Alfa*
   (MHH2*(-MHH2 + MZ2) + MA02*(MHH2 + MZ2))*B0i(bb0, MA02, MHH2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MA02*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MZ2*B0i(bb1, MA02, Mh02, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MZ2*B0i(bb1, MA02, MHH2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*MZ2*(-MA02 + Mh02 + MZ2)*
   C0i(cc0, m12, m22, m32, Mh02, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*MZ2*C0i(cc00, m12, m22, m32, Mh02, MZ2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.039788735772973836*Alfa*MZ2*C0i(cc00, m12, m22, m32, MHH2, MZ2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*MZ2*C0i(cc00, m32, m12, m22, MA02, MHH2, MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*(3*m12 + m22 - m32)*MZ2*
   C0i(cc11, m12, m22, m32, Mh02, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*(3*m12 + m22 - m32)*MZ2*
   C0i(cc11, m12, m22, m32, MHH2, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*(m12 - m22)*MZ2*C0i(cc12, m12, m22, m32, 
    Mh02, MZ2, MZ2)*Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*(m12 - m22)*MZ2*
   C0i(cc12, m12, m22, m32, MHH2, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(m12 - m22)*MZ2*C0i(cc12, m32, m12, m22, 
    MA02, MHH2, MZ2)*Cos(alp - beta)*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2) - (0.009947183943243459*Alfa*(m12 + 3*m22 - m32)*MZ2*
   C0i(cc22, m32, m12, m22, MA02, MHH2, MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*C0i(cc1, m12, m22, m32, Mh02, MZ2, MZ2)*
   Cos(alp - beta)*((MA02 - Mh02)*Mh02 + (MA02 - Mh02 + m12 - m22 + m32)*
     MZ2 + (pow(MZ2, 2)))*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2) - (0.0625*Alfa2*MZ2*(-MA02 + MHH2 + MZ2)*v2*
   C0i(cc0, m12, m22, m32, MHH2, MZ2, MZ2)*Cos(alp - beta)*
   (pow(Csc(w), 2))*(pow(Sin(alp - beta), 2)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*(-(Mh02*MHH2) + MA02*(Mh02 + MZ2) + 
    MZ2*(-MHH2 + m12 - m22 + m32 + MZ2))*v2*C0i(cc1, m12, m22, m32, 
    MHH2, MZ2, MZ2)*Cos(alp - beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*v2*C0i(cc2, m32, m12, m22, MA02, MHH2, MZ2)*
   Cos(alp - beta)*(Mh02*MHH2 - MA02*(Mh02 + MHH2) + 
    MZ2*(m12 - m22 - m32 + MZ2) + (pow(MA02, 2)))*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-46*MW2*MZ2 + 20*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(A0i(aa0, MB2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-70*MW2*MZ2 + 56*(pow(MW2, 2)) + 23*(pow(MZ2, 2)))*
   Re(A0i(aa0, MC2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-46*MW2*MZ2 + 20*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(A0i(aa0, MD2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-22*MW2*MZ2 + 12*(pow(MW2, 2)) + 9*(pow(MZ2, 2)))*
   Re(A0i(aa0, ME2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*Re(A0i(aa0, MHp2)))/
  (CW*MW2*MZ2*SW*SW2) - (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-22*MW2*MZ2 + 12*(pow(MW2, 2)) + 9*(pow(MZ2, 2)))*
   Re(A0i(aa0, ML2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-22*MW2*MZ2 + 12*(pow(MW2, 2)) + 9*(pow(MZ2, 2)))*
   Re(A0i(aa0, MM2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-46*MW2*MZ2 + 20*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(A0i(aa0, MS2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-70*MW2*MZ2 + 56*(pow(MW2, 2)) + 23*(pow(MZ2, 2)))*
   Re(A0i(aa0, MT2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-70*MW2*MZ2 + 56*(pow(MW2, 2)) + 23*(pow(MZ2, 2)))*
   Re(A0i(aa0, MU2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(4*MW2 - 3*MZ2)*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(A0i(aa0, MW2)))/(CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*Re(A0i(aa0, MZ2)))/
  (CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) + 
 (0.17904931097838225*Alfa*Cos(alp - beta)*Re(B0i(bb0, 0, MW2, MW2)))/
  (CW*SW) - (0.0049735919716217296*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb0, MA02, Mh02, MZ2)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MA02, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.009947183943243459*Alfa*Cos(alp - beta)*
   Re(B0i(bb0, MA02, MHp2, MW2)))/(CW*SW*SW2) + 
 (0.0024867959858108648*Alfa*Cos(alp - beta)*(Mh02*MHH2 - 2*Mh02*MZ2 + 
    MHH2*MZ2 - MA02*(Mh02 + MHH2 + MZ2) + Cos(2*(alp - beta))*
     (MHH2*(-Mh02 + MZ2) + MA02*(Mh02 + MHH2 + MZ2) - (pow(MA02, 2))) + 
    (pow(MA02, 2)))*Re(B0i(bb0, Mh02, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.0049735919716217296*Alfa*Cos(alp - beta)*
   ((Mh02 - MHp2)*(MHH2 - MHp2) + (-2*Mh02 + MHH2 - MHp2)*MW2 + 
    ((Mh02 - MHp2)*(-MHH2 + MHp2) + (MHH2 + MHp2)*MW2)*Cos(2*(alp - beta)))*
   Re(B0i(bb0, Mh02, MHp2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*(Mh02*MHH2 + 2*MW2*(-MHH2 + 6*MW2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, MW2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.0024867959858108648*Alfa*(Mh02*MHH2 + 2*MZ2*(-MHH2 + 6*MZ2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*Cos(alp - beta)*(MHH2*(-Mh02 + MZ2) + 
    MA02*(Mh02 + MHH2 + MZ2) - (pow(MA02, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MHH2, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*((Mh02 - MHp2)*(MHH2 - MHp2) - (MHH2 + MHp2)*MW2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MHp2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*(Mh02*MHH2 + 2*MW2*(-MHH2 + 6*MW2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MW2, MW2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0024867959858108648*Alfa*(Mh02*MHH2 + 2*MZ2*(-MHH2 + 6*MZ2))*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.07957747154594767*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, 0, MW2)))/(CW*MZ2*SW*SW2) + 
 (0.05968310365946075*Alfa*(MT2 - MW2)*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, MB2, MT2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.05968310365946075*Alfa*MC2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, MC2, MS2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.05968310365946075*Alfa*(MU2 - MW2)*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MW2, MD2, MU2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MW2, Mh02, MW2)))/
  (CW*MZ2*SW*(pow(SW2, 2))) - (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   (pow(Cos(alp - beta), 3))*Re(B0i(bb0, MW2, MHH2, MW2)))/
  (CW*MZ2*SW*(pow(SW2, 2))) + (0.019894367886486918*Alfa*Cos(alp - beta)*
   (6*MZ2*(pow(MW2, 2)) + 4*(pow(MW2, 3)) - 6*MW2*(pow(MZ2, 2)) + 
    (pow(MZ2, 3)))*Re(B0i(bb0, MW2, MW2, MZ2)))/
  (CW*MW2*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MB2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MB2, MB2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MC2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MC2, MC2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MD2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MD2, MD2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*ME2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, ME2, ME2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb0, MZ2, Mh02, MZ2)))/
  (CW*MW2*SW*(pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   (pow(Cos(alp - beta), 3))*Re(B0i(bb0, MZ2, MHH2, MZ2)))/
  (CW*MW2*SW*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*ML2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, ML2, ML2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*MM2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MM2, MM2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MS2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MS2, MS2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MT2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MT2, MT2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*MU2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb0, MZ2, MU2, MU2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-8*MW2*MZ2 + 13*(pow(MW2, 2)) + 2*(pow(MZ2, 2)))*
   Re(B0i(bb0, MZ2, MW2, MW2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.026525823848649224*Alfa*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MB2, MB2)))/(CW*MW2*MZ2*SW) - 
 (0.05305164769729845*Alfa*(8*MW2 - 5*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MC2, MC2)))/(CW*MW2*MZ2*SW) - 
 (0.026525823848649224*Alfa*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MD2, MD2)))/(CW*MW2*MZ2*SW) - 
 (0.07957747154594767*Alfa*(4*MW2 - 3*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, ME2, ME2)))/(CW*MW2*MZ2*SW) + 
 (0.07957747154594767*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MHp2, MHp2)))/(CW*MW2*MZ2*SW) - 
 (0.07957747154594767*Alfa*(4*MW2 - 3*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, ML2, ML2)))/(CW*MW2*MZ2*SW) - 
 (0.07957747154594767*Alfa*(4*MW2 - 3*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MM2, MM2)))/(CW*MW2*MZ2*SW) - 
 (0.026525823848649224*Alfa*(4*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MS2, MS2)))/(CW*MW2*MZ2*SW) - 
 (0.05305164769729845*Alfa*(8*MW2 - 5*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MT2, MT2)))/(CW*MW2*MZ2*SW) - 
 (0.05305164769729845*Alfa*(8*MW2 - 5*MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MU2, MU2)))/(CW*MW2*MZ2*SW) + 
 (0.07957747154594767*Alfa*(6*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, 0, MW2, MW2)))/(CW*MW2*MZ2*SW) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, ME2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, ML2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, MM2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.15915494309189535*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, 0, MW2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MA02, MHp2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.1193662073189215*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MB2, MT2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.1193662073189215*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MC2, MS2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.1193662073189215*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MD2, MU2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MW2, Mh02, MHp2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MW2, Mh02, MW2)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MW2, MHH2, MHp2)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MW2, MHH2, MW2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(8*MW2 + MZ2)*Cos(alp - beta)*
   Re(B0i(bb00, MW2, MW2, MZ2)))/(CW*MW2*SW*(pow(MZ2, 2))*
   (pow(SW2, 2))) + (0.05968310365946075*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*Re(B0i(bb00, MZ2, 0, 0)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MZ2, MA02, Mh02)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MZ2, MA02, MHH2)))/
  (CW*MW2*MZ2*SW*(pow(SW2, 2))) + 
 (0.006631455962162306*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MB2, MB2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.006631455962162306*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 
    17*(pow(MZ2, 2)))*Re(B0i(bb00, MZ2, MC2, MC2)))/
  (CW*MW2*SW*(pow(MZ2, 3))*(pow(SW2, 2))) + 
 (0.006631455962162306*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MD2, MD2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, ME2, ME2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) - (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb00, MZ2, Mh02, MZ2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb00, MZ2, MHH2, MZ2)))/(CW*MW2*MZ2*SW*(pow(SW2, 2))) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*(pow(2*MW2 - MZ2, 3))*
   Re(B0i(bb00, MZ2, MHp2, MHp2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, ML2, ML2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MM2, MM2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.006631455962162306*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MS2, MS2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.006631455962162306*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 
    17*(pow(MZ2, 2)))*Re(B0i(bb00, MZ2, MT2, MT2)))/
  (CW*MW2*SW*(pow(MZ2, 3))*(pow(SW2, 2))) + 
 (0.006631455962162306*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MU2, MU2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) - (0.019894367886486918*Alfa*(2*MW2 - MZ2)*
   Cos(alp - beta)*(-4*MW2*MZ2 + 12*(pow(MW2, 2)) + (pow(MZ2, 2)))*
   Re(B0i(bb00, MZ2, MW2, MW2)))/(CW*MW2*SW*(pow(MZ2, 3))*
   (pow(SW2, 2))) + (0.026525823848649224*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MB2, MB2)))/(CW*SW) + 
 (0.1061032953945969*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MC2, MC2)))/
  (CW*SW) + (0.026525823848649224*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MD2, MD2)))/(CW*SW) + 
 (0.07957747154594767*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, ME2, ME2)))/
  (CW*SW) + (0.07957747154594767*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, ML2, ML2)))/(CW*SW) + 
 (0.07957747154594767*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MM2, MM2)))/
  (CW*SW) + (0.026525823848649224*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MS2, MS2)))/(CW*SW) + 
 (0.1061032953945969*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MT2, MT2)))/
  (CW*SW) + (0.1061032953945969*Alfa*Cos(alp - beta)*
   Re(B0i(bb1, 0, MU2, MU2)))/(CW*SW) + 
 (0.039788735772973836*Alfa*Cos(alp - beta)*Re(B0i(bb1, 0, MW2, MW2)))/
  (CW*SW) - (0.05968310365946075*Alfa*MB2*Cos(alp - beta)*
   (pow(Tan(beta), 2))*Re(B0i(bb1, MA02, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MC2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MD2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(bb1, MA02, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*ME2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, ME2, ME2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(bb1, MA02, Mh02, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MZ2*Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb1, MA02, MHH2, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*Re(B0i(bb1, MA02, MHp2, MW2)))/
  (CW*SW*SW2) - (0.019894367886486918*Alfa*ML2*Cos(alp - beta)*
   (pow(Cot(beta), 2))*Re(B0i(bb1, MA02, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MM2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MS2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(bb1, MA02, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MT2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MU2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb1, MA02, MU2, MU2)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*MZ2*Cos(alp - beta)*
   (-2*Mh02 + MHH2 + MHH2*Cos(2*(alp - beta)))*Re(B0i(bb1, Mh02, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-2*Mh02 + MHH2 + MHH2*Cos(2*(alp - beta)))*Re(B0i(bb1, Mh02, MHp2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) - (0.019894367886486918*Alfa*MHH2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, Mh02, MW2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) - (0.009947183943243459*Alfa*MHH2*MZ2*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb1, Mh02, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MHH2*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, MHH2, MA02, MZ2)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MHH2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, MHH2, MHp2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) + (0.019894367886486918*Alfa*MHH2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb1, MHH2, MW2, MW2)))/
  (CW*(Mh02 - MHH2)*SW*SW2) + (0.009947183943243459*Alfa*MHH2*MZ2*
   Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb1, MHH2, MZ2, MZ2)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, ME2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, ML2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.019894367886486918*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, MM2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.039788735772973836*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, 0, MW2)))/(CW*MZ2*SW*SW2) - 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MB2, MT2)))/(CW*MZ2*SW*(pow(SW2, 2))) + 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MC2, MS2)))/(CW*MZ2*SW*(pow(SW2, 2))) - 
 (0.05968310365946075*Alfa*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MD2, MU2)))/(CW*MZ2*SW*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*MW2*(2*MW2 - MZ2)*Cos(alp - beta)*
   Re(B0i(bb1, MW2, MW2, MZ2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.029841551829730376*Alfa*Cos(alp - beta)*Re(B0i(bb1, MZ2, 0, 0)))/
  (CW*SW*(pow(SW2, 2))) - (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MC2, MC2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, ME2, ME2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, ML2, ML2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MM2, MM2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MT2, MT2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(bb1, MZ2, MU2, MU2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) - 
 (0.039788735772973836*Alfa*Cos(alp - beta)*(pow(MW2, 2))*
   Re(B0i(bb1, MZ2, MW2, MW2)))/(CW*SW*(pow(MZ2, 2))*(pow(SW2, 2))) + 
 (0.0012433979929054324*Alfa*Cos(alp - beta)*
   (pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*Mh02)*
       Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, MA02, MA02, Mh02)))/(CW*MW2*SW*SW2) + 
 (0.0012433979929054324*Alfa*Cos(alp - beta)*(pow(Csc(2*beta), 2))*
   (pow((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
       Sin(alp + beta), 2))*Re(B0i(dbb0, MA02, MA02, MHH2)))/
  (CW*MW2*SW*SW2) - (0.0049735919716217296*Alfa*
   ((MA02 + Mh02)*MZ2 - (pow(MA02 - Mh02, 2)))*
   (pow(Cos(alp - beta), 3))*Re(B0i(dbb0, MA02, Mh02, MZ2)))/
  (CW*MW2*SW*SW2) - (0.0049735919716217296*Alfa*Cos(alp - beta)*
   ((MA02 + MHH2)*MZ2 - (pow(MA02 - MHH2, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, MA02, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) + (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-((MA02 + MHp2)*MW2) + (pow(MA02 - MHp2, 2)))*
   Re(B0i(dbb0, MA02, MHp2, MW2)))/(CW*MW2*SW*SW2) + 
 (0.0006216989964527162*Alfa*Cos(alp - beta)*
   (pow((-2*MA02 + Mh02)*Cos(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*Mh02)*
       Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh02, MA02, MA02)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*((MA02 + Mh02)*MZ2 - (pow(MA02 - Mh02, 2)))*
   (pow(Cos(alp - beta), 3))*Re(B0i(dbb0, Mh02, MA02, MZ2)))/
  (CW*MW2*SW*SW2) + (0.005595290968074445*Alfa*Cos(alp - beta)*
   (pow(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
      (2*M2 - 3*Mh02)*Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh02, Mh02, Mh02)))/(CW*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*(pow(Cos(alp - beta), 3))*
   (pow(M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*Sin(2*alp), 2))*
   Re(B0i(dbb0, Mh02, Mh02, MHH2)))/(CW*MW2*SW*SW2) + 
 (0.0024867959858108648*Alfa*Cos(alp - beta)*(pow(Csc(2*beta), 2))*
   (pow(Sin(alp - beta), 2))*
   (pow((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta), 2))*
   Re(B0i(dbb0, Mh02, MHH2, MHH2)))/(CW*MW2*SW*SW2) + 
 (0.0012433979929054324*Alfa*Cos(alp - beta)*
   (pow((Mh02 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh02 + 2*MHp2)*
       Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh02, MHp2, MHp2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*(-((Mh02 + MHp2)*MW2) + 
    (pow(Mh02 - MHp2, 2)))*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb0, Mh02, MHp2, MW2)))/(CW*MW2*SW*SW2) + 
 (0.0049735919716217296*Alfa*Cos(alp - beta)*
   (-2*Mh02*MW2 + (pow(Mh02, 2)) + 12*(pow(MW2, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, Mh02, MW2, MW2)))/
  (CW*MW2*SW*SW2) + (0.0024867959858108648*Alfa*Cos(alp - beta)*
   (-2*Mh02*MZ2 + (pow(Mh02, 2)) + 12*(pow(MZ2, 2)))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, Mh02, MZ2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.029841551829730376*Alfa*MB2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MC2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MD2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*ME2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, ME2, ME2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*(pow(MZ2, 2))*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb0, MZ2, Mh02, MZ2)))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*(pow(MZ2, 2))*
   (pow(Cos(alp - beta), 3))*Re(B0i(dbb0, MZ2, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.009947183943243459*Alfa*ML2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*MM2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MS2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MT2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MU2*MZ2*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MU2, MU2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*(9*MW2 - 2*MZ2)*Cos(alp - beta)*
   Re(B0i(dbb0, MZ2, MW2, MW2)))/(CW*SW*SW2) - 
 (0.05305164769729845*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MB2, MB2)))/
  (CW*SW) - (0.2122065907891938*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, MC2, MC2)))/(CW*SW) - 
 (0.05305164769729845*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MD2, MD2)))/
  (CW*SW) - (0.15915494309189535*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, ME2, ME2)))/(CW*SW) + 
 (0.07957747154594767*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MHp2, MHp2)))/
  (CW*SW) - (0.15915494309189535*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, ML2, ML2)))/(CW*SW) - 
 (0.15915494309189535*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MM2, MM2)))/
  (CW*SW) - (0.05305164769729845*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, MS2, MS2)))/(CW*SW) - 
 (0.2122065907891938*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MT2, MT2)))/
  (CW*SW) - (0.2122065907891938*Alfa*Cos(alp - beta)*
   Re(B0i(dbb00, 0, MU2, MU2)))/(CW*SW) + 
 (0.238732414637843*Alfa*Cos(alp - beta)*Re(B0i(dbb00, 0, MW2, MW2)))/
  (CW*SW) + (0.05968310365946075*Alfa*MZ2*Cos(alp - beta)*
   Re(B0i(dbb00, MZ2, 0, 0)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb00, MZ2, MA02, Mh02)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(dbb00, MZ2, MA02, MHH2)))/(CW*MW2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MB2, MB2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MC2, MC2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MD2, MD2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, ME2, ME2)))/(CW*MW2*MZ2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*Cos(alp - beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(dbb00, MZ2, Mh02, MZ2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb00, MZ2, MHH2, MZ2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*(pow(-2*MW2 + MZ2, 2))*
   Re(B0i(dbb00, MZ2, MHp2, MHp2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, ML2, ML2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MM2, MM2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MS2, MS2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MT2, MT2)))/(CW*MW2*MZ2*SW*SW2) + 
 (0.006631455962162306*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MU2, MU2)))/(CW*MW2*MZ2*SW*SW2) - 
 (0.019894367886486918*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 12*(pow(MW2, 2)) + (pow(MZ2, 2)))*
   Re(B0i(dbb00, MZ2, MW2, MW2)))/(CW*MW2*MZ2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MB2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(dbb1, MA02, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MC2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MD2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(dbb1, MA02, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MA02*ME2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, ME2, ME2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MA02*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb1, MA02, Mh02, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MA02*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb1, MA02, MHH2, MZ2)))/
  (CW*MW2*SW*SW2) + (0.019894367886486918*Alfa*MA02*Cos(alp - beta)*
   Re(B0i(dbb1, MA02, MHp2, MW2)))/(CW*SW*SW2) - 
 (0.019894367886486918*Alfa*MA02*ML2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*MA02*MM2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MS2*Cos(alp - beta)*(pow(Tan(beta), 2))*
   Re(B0i(dbb1, MA02, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MT2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MA02*MU2*Cos(alp - beta)*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MA02, MU2, MU2)))/(CW*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*Mh02*MZ2*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb1, Mh02, MA02, MZ2)))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*Mh02*(pow(Cos(alp - beta), 3))*
   Re(B0i(dbb1, Mh02, MHp2, MW2)))/(CW*SW*SW2) + 
 (0.019894367886486918*Alfa*Mh02*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb1, Mh02, MW2, MW2)))/(CW*SW*SW2) + 
 (0.009947183943243459*Alfa*Mh02*MZ2*Cos(alp - beta)*
   (pow(Sin(alp - beta), 2))*Re(B0i(dbb1, Mh02, MZ2, MZ2)))/
  (CW*MW2*SW*SW2) - (0.029841551829730376*Alfa*Cos(alp - beta)*
   (pow(MZ2, 2))*Re(B0i(dbb1, MZ2, 0, 0)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MB2, MB2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MC2, MC2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MD2, MD2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, ME2, ME2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, ML2, ML2)))/(CW*MW2*SW*SW2) - 
 (0.009947183943243459*Alfa*Cos(alp - beta)*
   (-12*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MM2, MM2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-4*MW2*MZ2 + 8*(pow(MW2, 2)) + 5*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MS2, MS2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MT2, MT2)))/(CW*MW2*SW*SW2) - 
 (0.003315727981081153*Alfa*Cos(alp - beta)*
   (-40*MW2*MZ2 + 32*(pow(MW2, 2)) + 17*(pow(MZ2, 2)))*
   Re(B0i(dbb1, MZ2, MU2, MU2)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*MW2*Cos(alp - beta)*Re(B0i(dbb1, MZ2, MW2, MW2)))/
  (CW*SW*SW2) + (0.05968310365946075*Alfa*MC2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MC2, MC2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*ME2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, ME2, ME2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*ML2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, ML2, ML2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*MM2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MM2, MM2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MT2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MT2, MT2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MU2*Cos(alp)*Csc(beta)*
   Re(B0i(bb1, Mh02, MU2, MU2))*
   (Cot(beta)*(-Mh02 + MHH2*(pow(Cos(alp), 2))) + 
    MHH2*Cos(alp)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*(pow(MC2, 2))*
   Re(B0i(bb0, Mh02, MC2, MC2))*Sin(alp)*(Cos(alp) - Cot(beta)*Sin(alp)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*
   (pow(MT2, 2))*Re(B0i(bb0, Mh02, MT2, MT2))*Sin(alp)*
   (Cos(alp) - Cot(beta)*Sin(alp)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.1193662073189215*Alfa*Cos(alp)*Csc(beta)*(pow(MU2, 2))*
   Re(B0i(bb0, Mh02, MU2, MU2))*Sin(alp)*(Cos(alp) - Cot(beta)*Sin(alp)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (1.*Re(B0i(bb0, Mh02, ME2, ME2))*
   (-(Alfa*Cos(alp)*Cot(beta)*Csc(beta)*(pow(ME2, 2))*
      (pow(Sin(alp), 2))) + Alfa*Csc(beta)*(pow(ME2, 2))*
     (pow(Cos(alp), 2))*Sin(alp)))/(8*CW*Mh02*MW2*Pi*SW*SW2 - 
   8*CW*MHH2*MW2*Pi*SW*SW2) + (1.*Re(B0i(bb0, Mh02, ML2, ML2))*
   (-(Alfa*Cos(alp)*Cot(beta)*Csc(beta)*(pow(ML2, 2))*
      (pow(Sin(alp), 2))) + Alfa*Csc(beta)*(pow(ML2, 2))*
     (pow(Cos(alp), 2))*Sin(alp)))/(8*CW*Mh02*MW2*Pi*SW*SW2 - 
   8*CW*MHH2*MW2*Pi*SW*SW2) + (1.*Re(B0i(bb0, Mh02, MM2, MM2))*
   (-(Alfa*Cos(alp)*Cot(beta)*Csc(beta)*(pow(MM2, 2))*
      (pow(Sin(alp), 2))) + Alfa*Csc(beta)*(pow(MM2, 2))*
     (pow(Cos(alp), 2))*Sin(alp)))/(8*CW*Mh02*MW2*Pi*SW*SW2 - 
   8*CW*MHH2*MW2*Pi*SW*SW2) + (0.03125*Alfa2*(-MA02 + Mh02 + MZ2)*v2*
   C0i(cc0, m22, m32, m12, Mh02, MHH2, MZ2)*(pow(Cos(alp - beta), 3))*
   (pow(Csc(w), 2))*(M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*
     Sin(2*alp)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.03125*Alfa2*(MA02 - Mh02 + MZ2)*v2*C0i(cc1, m22, m32, m12, Mh02, MHH2, 
    MZ2)*(pow(Cos(alp - beta), 3))*(pow(Csc(w), 2))*
   (M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*Sin(2*alp)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.03125*Alfa2*(MA02 - Mh02 + MZ2)*v2*C0i(cc2, m22, m32, m12, Mh02, MHH2, 
    MZ2)*(pow(Cos(alp - beta), 3))*(pow(Csc(w), 2))*
   (M2 + (-3*M2 + 2*Mh02 + MHH2)*Csc(2*beta)*Sin(2*alp)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.05968310365946075*Alfa*MC2*B0i(bb1, MA02, MC2, MC2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*ME2*B0i(bb1, MA02, ME2, ME2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*ML2*B0i(bb1, MA02, ML2, ML2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*MM2*B0i(bb1, MA02, MM2, MM2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MT2*B0i(bb1, MA02, MT2, MT2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MU2*B0i(bb1, MA02, MU2, MU2)*Cot(beta)*
   Sin(alp - beta))/(CW*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*(pow(MB2, 2))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, Mh02, MB2, MB2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.05968310365946075*Alfa*(pow(MD2, 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh02, MD2, MD2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*(pow(MS2, 2))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, Mh02, MS2, MS2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - (0.05968310365946075*Alfa*(pow(MB2, 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MB2, MB2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*(pow(MC2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MC2, MC2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - (0.05968310365946075*Alfa*(pow(MD2, 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MD2, MD2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(pow(ME2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, ME2, ME2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(pow(ML2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, ML2, ML2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.019894367886486918*Alfa*(pow(MM2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MM2, MM2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - (0.05968310365946075*Alfa*(pow(MS2, 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MS2, MS2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*(pow(MT2, 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MT2, MT2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.05968310365946075*Alfa*(pow(MU2, 2))*
   (pow(Csc(beta), 2))*Re(B0i(bb0, MHH2, MU2, MU2))*Sin(2*alp)*
   Sin(alp - beta))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MB2*MHH2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, MHH2, MB2, MB2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MC2*MHH2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MC2, MC2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MD2*MHH2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, MHH2, MD2, MD2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*ME2*MHH2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, ME2, ME2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MHH2*ML2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, ML2, ML2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.009947183943243459*Alfa*MHH2*MM2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MM2, MM2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.029841551829730376*Alfa*MHH2*MS2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, MHH2, MS2, MS2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MHH2*MT2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MT2, MT2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.029841551829730376*Alfa*MHH2*MU2*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MU2, MU2))*Sin(2*alp)*Sin(alp - beta))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.0012433979929054324*Alfa*(MA02 - Mh02)*
   B0i(bb0, 0, MA02, Mh02)*((2*MA02 - Mh02)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh02)*Cos(alp + beta))*Csc(2*beta)*
   Sin(2*(alp - beta)))/(CW*MA02*MW2*SW*SW2) - 
 (0.0012433979929054324*Alfa*(MA02 - Mh02)*B0i(bb0, MA02, MA02, Mh02)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(2*beta)*Sin(2*(alp - beta)))/(CW*MA02*MW2*SW*SW2) + 
 (0.003730193978716297*Alfa*(-MA02 + Mh02 + MZ2)*C0i(cc0, m22, m32, m12, 
    Mh02, Mh02, MZ2)*(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh02)*Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) - (0.003730193978716297*Alfa*(MA02 - Mh02 + MZ2)*
   C0i(cc1, m22, m32, m12, Mh02, Mh02, MZ2)*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh02)*Cos(3*alp - beta) + (2*M2 - 3*Mh02)*Cos(alp + beta))*
   Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/(CW*MW2*SW*SW2) - 
 (0.003730193978716297*Alfa*(MA02 - Mh02 + MZ2)*C0i(cc2, m22, m32, m12, 
    Mh02, Mh02, MZ2)*(M2*Cos(alp - 3*beta) + (M2 - Mh02)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh02)*Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) + (0.0012433979929054324*Alfa*(MA02 - Mh02 + MZ2)*
   C0i(cc0, m12, m32, m22, MA02, Mh02, MZ2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) - (0.0012433979929054324*Alfa*(-MA02 + Mh02 + MZ2)*
   C0i(cc1, m12, m32, m22, MA02, Mh02, MZ2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) - (0.0012433979929054324*Alfa*(-MA02 + Mh02 + MZ2)*
   C0i(cc2, m12, m32, m22, MA02, Mh02, MZ2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*Csc(beta)*Sec(beta)*Sin(2*(alp - beta)))/
  (CW*MW2*SW*SW2) + (0.05968310365946075*Alfa*MB2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, Mh02, MB2, MB2))*Sin(alp)*(MHH2*Cos(alp - beta)*Sin(alp) - 
    Mh02*Sin(beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MD2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, Mh02, MD2, MD2))*Sin(alp)*(MHH2*Cos(alp - beta)*Sin(alp) - 
    Mh02*Sin(beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.05968310365946075*Alfa*MS2*(pow(Sec(beta), 2))*
   Re(B0i(bb1, Mh02, MS2, MS2))*Sin(alp)*(MHH2*Cos(alp - beta)*Sin(alp) - 
    Mh02*Sin(beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.03125*Alfa2*(-MA02 + MHH2 + MZ2)*v2*C0i(cc0, m22, m32, m12, MHH2, 
    MHH2, MZ2)*Cos(alp - beta)*Csc(2*beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2))*((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*(MA02 - MHH2 + MZ2)*v2*C0i(cc1, m22, m32, m12, MHH2, MHH2, 
    MZ2)*Cos(alp - beta)*Csc(2*beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2))*((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.03125*Alfa2*(MA02 - MHH2 + MZ2)*v2*C0i(cc2, m22, m32, m12, MHH2, MHH2, 
    MZ2)*Cos(alp - beta)*Csc(2*beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2))*((3*M2 - Mh02 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.009947183943243459*Alfa*MZ2*C0i(cc0, m22, m12, m32, Mh02, MHH2, MZ2)*
   Cos(alp - beta)*Csc(beta)*(pow(Sin(alp - beta), 2))*Sec(beta)*
   ((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(CW*MW2*SW*SW2) - 
 (0.0049735919716217296*Alfa*(MA02 - MHH2 + MZ2)*C0i(cc1, m22, m12, m32, 
    Mh02, MHH2, MZ2)*Cos(alp - beta)*Csc(beta)*(pow(Sin(alp - beta), 2))*
   Sec(beta)*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*MW2*SW*SW2) + (0.0078125*Alfa2*v2*C0i(cc1, m12, m22, m32, MA02, 
    Mh02, MHH2)*((2*MA02 - Mh02)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh02)*Cos(alp + beta))*(pow(Csc(2*beta), 2))*
   (pow(Csc(w), 2))*Sin(2*(alp - beta))*
   ((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0018650969893581485*Alfa*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh02)*Cos(3*alp - beta) + (2*M2 - 3*Mh02)*Cos(alp + beta))*
   (pow(Csc(2*beta), 2))*Re(B0i(bb0, Mh02, Mh02, Mh02))*
   Sin(2*(alp - beta))*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0018650969893581485*Alfa*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh02)*Cos(3*alp - beta) + (2*M2 - 3*Mh02)*Cos(alp + beta))*
   (pow(Csc(2*beta), 2))*Re(B0i(bb0, MHH2, Mh02, Mh02))*
   Sin(2*(alp - beta))*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) - (0.0024867959858108648*Alfa*Cos(alp - beta)*
   (pow(Csc(2*beta), 2))*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, Mh02, MHH2))*(9*M2*(Mh02 + MHH2) - 
    (2*Mh02 + MHH2)*(Mh02 + 2*MHH2) + (3*M2 - Mh02 - 2*MHH2)*
     (3*M2 - 2*Mh02 - MHH2)*Cos(4*alp) - 8*(pow(M2, 2)) + 
    M2*(-(M2*Cos(4*beta)) + 2*(Mh02 - MHH2)*Sin(2*alp)*Sin(2*beta))))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + (0.0024867959858108648*Alfa*Cos(alp - beta)*
   (pow(Csc(2*beta), 2))*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, Mh02, MHH2))*(9*M2*(Mh02 + MHH2) - 
    (2*Mh02 + MHH2)*(Mh02 + 2*MHH2) + (3*M2 - Mh02 - 2*MHH2)*
     (3*M2 - 2*Mh02 - MHH2)*Cos(4*alp) - 8*(pow(M2, 2)) + 
    M2*(-(M2*Cos(4*beta)) + 2*(Mh02 - MHH2)*Sin(2*alp)*Sin(2*beta))))/
  (CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.003730193978716297*Alfa*Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, Mh02, MHH2, MHH2))*(M2 + (3*M2 - Mh02 - 2*MHH2)*Csc(2*beta)*
     Sin(2*alp))*(M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.003730193978716297*Alfa*Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb0, MHH2, MHH2, MHH2))*(M2 + (3*M2 - Mh02 - 2*MHH2)*Csc(2*beta)*
     Sin(2*alp))*(M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.015625*Alfa2*v2*C0i(cc1, m12, m22, m32, MA02, MHH2, MHH2)*Csc(2*beta)*
   (pow(Csc(w), 2))*(pow(Sin(alp - beta), 2))*
   (M2 + (3*M2 - Mh02 - 2*MHH2)*Csc(2*beta)*Sin(2*alp))*
   ((2*MA02 - MHH2)*Sin(alp - 3*beta) + (4*M2 - 2*MA02 - 3*MHH2)*
     Sin(alp + beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.015625*Alfa2*v2*C0i(cc2, m32, m22, m12, MA02, Mh02, MHH2)*
   (pow(Cos(alp - beta), 2))*(pow(Csc(2*beta), 2))*
   (pow(Csc(w), 2))*((-3*M2 + 2*Mh02 + MHH2)*Sin(2*alp) + M2*Sin(2*beta))*
   ((2*MA02 - MHH2)*Sin(alp - 3*beta) + (4*M2 - 2*MA02 - 3*MHH2)*
     Sin(alp + beta)))/(CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0024867959858108648*Alfa*(MA02 - MHH2)*B0i(bb0, 0, MA02, MHH2)*
   Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*MA02*MW2*SW*SW2) + 
 (0.0024867959858108648*Alfa*(MA02 - MHH2)*B0i(bb0, MA02, MA02, MHH2)*
   Csc(2*beta)*(pow(Sin(alp - beta), 2))*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*MA02*MW2*SW*SW2) + 
 (0.0078125*Alfa2*(MA02 - Mh02 + MZ2)*v2*C0i(cc0, m12, m32, m22, MA02, 
    MHH2, MZ2)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Sec(beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*(-MA02 + Mh02 + MZ2)*v2*C0i(cc1, m12, m32, m22, MA02, 
    MHH2, MZ2)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Sec(beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) - 
 (0.0078125*Alfa2*(-MA02 + Mh02 + MZ2)*v2*C0i(cc2, m12, m32, m22, MA02, 
    MHH2, MZ2)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Sec(beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0078125*Alfa2*v2*C0i(cc0, m22, m32, m12, MA02, MA02, MHH2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*(pow(Csc(2*beta), 2))*(pow(Csc(w), 2))*
   Sin(alp - beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0078125*Alfa2*v2*C0i(cc1, m22, m32, m12, MA02, MA02, MHH2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*(pow(Csc(2*beta), 2))*(pow(Csc(w), 2))*
   Sin(alp - beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.0078125*Alfa2*v2*C0i(cc2, m22, m32, m12, MA02, MA02, MHH2)*
   ((2*MA02 - Mh02)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh02)*
     Cos(alp + beta))*(pow(Csc(2*beta), 2))*(pow(Csc(w), 2))*
   Sin(alp - beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/
  (CW*MW2*SW*SW2*(pow(MW, 2))) + 
 (0.00015542474911317905*Alfa*((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh02, MA02, MA02))*Sin(alp - beta)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.00015542474911317905*Alfa*((-2*MA02 + Mh02)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh02)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MA02, MA02))*Sin(alp - beta)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) + 
 (0.0003108494982263581*Alfa*((Mh02 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh02 + 2*MHp2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh02, MHp2, MHp2))*Sin(alp - beta)*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.0003108494982263581*Alfa*((Mh02 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh02 + 2*MHp2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MHp2, MHp2))*Sin(alp - beta)*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/(CW*(Mh02 - MHH2)*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MC2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MC2, MC2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*Csc(beta)*(pow(ME2, 2))*
   (pow(Cos(alp), 3))*Re(B0i(dbb0, Mh02, ME2, ME2))*
   (Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*Csc(beta)*(pow(ML2, 2))*
   (pow(Cos(alp), 3))*Re(B0i(dbb0, Mh02, ML2, ML2))*
   (Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.039788735772973836*Alfa*Csc(beta)*(pow(MM2, 2))*
   (pow(Cos(alp), 3))*Re(B0i(dbb0, Mh02, MM2, MM2))*
   (Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MT2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MT2, MT2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*Csc(beta)*(pow(MU2, 2))*(pow(Cos(alp), 3))*
   Re(B0i(dbb0, Mh02, MU2, MU2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MC2*Mh02*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MC2, MC2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*ME2*Mh02*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, ME2, ME2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*Mh02*ML2*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, ML2, ML2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.019894367886486918*Alfa*Mh02*MM2*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MM2, MM2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*Mh02*MT2*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MT2, MT2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*Mh02*MU2*Csc(beta)*(pow(Cos(alp), 3))*
   Re(B0i(dbb1, Mh02, MU2, MU2))*(Cot(beta) + Tan(alp)))/(CW*MW2*SW*SW2) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MB2, MB2, MB2)*(pow(MB2, 2))*
   Sec(beta)*Sin(alp)*Tan(beta))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MD2, MD2, MD2)*(pow(MD2, 2))*
   Sec(beta)*Sin(alp)*Tan(beta))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MZ*MZ2*C0i(cc0, m22, m32, m12, MS2, MS2, MS2)*(pow(MS2, 2))*
   Sec(beta)*Sin(alp)*Tan(beta))/(-8*MW*MW2*MZ2*Pi*SW + 
   8*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MB2*MZ*MZ2*B0i(bb0, m32, MB2, MB2)*Sec(beta)*Sin(alp)*Tan(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MD2*MZ*MZ2*B0i(bb0, m32, MD2, MD2)*Sec(beta)*Sin(alp)*Tan(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MS2*MZ*MZ2*B0i(bb0, m32, MS2, MS2)*Sec(beta)*Sin(alp)*Tan(beta))/
  (-16*MW*MW2*MZ2*Pi*SW + 16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MB2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MB2, MB2, MB2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MD2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MD2, MD2, MD2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MS2*m22*MZ*MZ2*C0i(cc1, m22, m32, m12, MS2, MS2, MS2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MB2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MB2, MB2, MB2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MD2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MD2, MD2, MD2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) + 
 (3.*Alfa*MS2*m12*MZ*MZ2*C0i(cc2, m22, m32, m12, MS2, MS2, MS2)*Sec(beta)*
   Sin(alp)*Tan(beta))/(-16*MW*MW2*MZ2*Pi*SW + 
   16*MW*Pi*SW*(pow(MW2, 2))) - 
 (0.05968310365946075*Alfa*MB2*B0i(bb1, MA02, MB2, MB2)*Sin(alp - beta)*
   Tan(beta))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MD2*B0i(bb1, MA02, MD2, MD2)*Sin(alp - beta)*
   Tan(beta))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MS2*B0i(bb1, MA02, MS2, MS2)*Sin(alp - beta)*
   Tan(beta))/(CW*MW2*SW*SW2) - (0.1193662073189215*Alfa*(pow(MB2, 2))*
   (pow(Sin(alp), 3))*Re(B0i(dbb0, Mh02, MB2, MB2))*Sec(beta)*
   (Cot(alp) + Tan(beta)))/(CW*MW2*SW*SW2) - 
 (0.1193662073189215*Alfa*(pow(MD2, 2))*(pow(Sin(alp), 3))*
   Re(B0i(dbb0, Mh02, MD2, MD2))*Sec(beta)*(Cot(alp) + Tan(beta)))/
  (CW*MW2*SW*SW2) - (0.1193662073189215*Alfa*(pow(MS2, 2))*
   (pow(Sin(alp), 3))*Re(B0i(dbb0, Mh02, MS2, MS2))*Sec(beta)*
   (Cot(alp) + Tan(beta)))/(CW*MW2*SW*SW2) - 
 (0.05968310365946075*Alfa*MB2*Mh02*(pow(Sin(alp), 3))*
   Re(B0i(dbb1, Mh02, MB2, MB2))*Sec(beta)*(Cot(alp) + Tan(beta)))/
  (CW*MW2*SW*SW2) - (0.05968310365946075*Alfa*MD2*Mh02*(pow(Sin(alp), 3))*
   Re(B0i(dbb1, Mh02, MD2, MD2))*Sec(beta)*(Cot(alp) + Tan(beta)))/
  (CW*MW2*SW*SW2) - (0.05968310365946075*Alfa*Mh02*MS2*(pow(Sin(alp), 3))*
   Re(B0i(dbb1, Mh02, MS2, MS2))*Sec(beta)*(Cot(alp) + Tan(beta)))/
  (CW*MW2*SW*SW2)
;
}

} //end namespace TypeFL

