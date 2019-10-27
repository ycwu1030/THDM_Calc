#include "CouplingFunctionTypeLS.h"

namespace TypeLS{
ComplexType dKappahWW(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32)
{
return 0. + (0.75*MT2*A0i(aa0, MT2)*Cos(alp - beta)*Cot(beta))/(MA02*Pi2*v2) + 
 (0.375*MT2*B0i(bb1, MA02, MT2, MT2)*Cos(alp - beta)*Cot(beta))/(Pi2*v2) + 
 (0.75*MT2*(-m22 + MT2)*C0i(cc0, m22, m12, m32, MB2, MT2, MT2)*
   (-1 + Cos(alp)*Csc(beta)))/(Pi2*v2) - 
 (1.5*MT2*C0i(cc00, m22, m12, m32, MB2, MT2, MT2)*(-1 + Cos(alp)*Csc(beta)))/
  (Pi2*v2) - (0.1875*(m12 + 5*m22 - m32)*MT2*C0i(cc1, m22, m12, m32, MB2, 
    MT2, MT2)*(-1 + Cos(alp)*Csc(beta)))/(Pi2*v2) + 
 (0.1875*(3*m12 - 3*m22 - m32)*MT2*C0i(cc2, m22, m12, m32, MB2, MT2, MT2)*
   (-1 + Cos(alp)*Csc(beta)))/(Pi2*v2) - 
 (0.0234375*A0i(aa0, MA02)*Cos(alp - beta)*(3*(Mh2 - MHH2)*Cos(2*alp) + 
    (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 4*(-2*M2 + Mh2 + MHH2)*Cos(2*beta))*
   Csc(2*beta))/(MA02*Pi2*v2) - (0.015625*A0i(aa0, MHp2)*Cos(alp - beta)*
   (3*(Mh2 - MHH2)*Cos(2*alp) + (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 
    4*(-2*M2 + Mh2 + MHH2)*Cos(2*beta))*Csc(2*beta))/(MA02*Pi2*v2) + 
 (0.015625*B0i(bb0, m12, MA02, MA02)*((2*MA02 - Mh2)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh2)*Cos(alp + beta))*Csc(2*beta))/(Pi2*v2) + 
 (0.0625*C0i(cc00, m12, m32, m22, MA02, MA02, MHp2)*
   ((-2*MA02 + Mh2)*Cos(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*Mh2)*
     Cos(alp + beta))*Csc(2*beta))/(Pi2*v2) - 
 (0.03125*B0i(bb0, m12, MHp2, MHp2)*((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh2 + 2*MHp2)*Cos(alp + beta))*Csc(2*beta))/(Pi2*v2) + 
 (0.0625*C0i(cc00, m22, m12, m32, MA02, MHp2, MHp2)*
   ((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*Csc(2*beta))/(Pi2*v2) + 
 (0.75*A0i(aa0, MB2)*Cos(alp - beta)*Cot(beta)*(pow(MB, 2)))/
  (MA02*Pi2*v2) + (0.375*B0i(bb1, MA02, MB2, MB2)*Cos(alp - beta)*Cot(beta)*
   (pow(MB, 2)))/(Pi2*v2) - (1.5*C0i(cc00, m12, m32, m22, MB2, MB2, MT2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MB, 2)))/(Pi2*v2) + 
 (0.375*(2*m12 + m22 - m32)*C0i(cc1, m12, m32, m22, MB2, MB2, MT2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MB, 2)))/(Pi2*v2) + 
 (0.1875*(m12 + 5*m22 - m32)*C0i(cc2, m12, m32, m22, MB2, MB2, MT2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MB, 2)))/(Pi2*v2) + 
 (0.75*B0i(bb0, m32, MB2, MT2)*(-1 + Cos(alp)*Csc(beta))*
   (MT2 + (pow(MB, 2))))/(Pi2*v2) + 
 (0.1875*C0i(cc0, m12, m32, m22, MB2, MB2, MT2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MB, 2))*(m12 + m22 - m32 + 4*(pow(MB, 2))))/(Pi2*v2) + 
 (0.75*A0i(aa0, MC2)*Cos(alp - beta)*Cot(beta)*(pow(MC, 2)))/
  (MA02*Pi2*v2) + (0.375*B0i(bb1, MA02, MC2, MC2)*Cos(alp - beta)*Cot(beta)*
   (pow(MC, 2)))/(Pi2*v2) - (1.5*C0i(cc00, m12, m32, m22, MC2, MC2, MS2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MC, 2)))/(Pi2*v2) + 
 (0.375*(2*m12 + m22 - m32)*C0i(cc1, m12, m32, m22, MC2, MC2, MS2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MC, 2)))/(Pi2*v2) + 
 (0.1875*(m12 + 5*m22 - m32)*C0i(cc2, m12, m32, m22, MC2, MC2, MS2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MC, 2)))/(Pi2*v2) + 
 (0.1875*C0i(cc0, m12, m32, m22, MC2, MC2, MS2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MC, 2))*(m12 + m22 - m32 + 4*(pow(MC, 2))))/(Pi2*v2) + 
 (0.75*A0i(aa0, MD2)*Cos(alp - beta)*Cot(beta)*(pow(MD, 2)))/
  (MA02*Pi2*v2) + (0.375*B0i(bb1, MA02, MD2, MD2)*Cos(alp - beta)*Cot(beta)*
   (pow(MD, 2)))/(Pi2*v2) - (1.5*C0i(cc00, m12, m32, m22, MD2, MD2, MU2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MD, 2)))/(Pi2*v2) + 
 (0.375*(2*m12 + m22 - m32)*C0i(cc1, m12, m32, m22, MD2, MD2, MU2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MD, 2)))/(Pi2*v2) + 
 (0.1875*(m12 + 5*m22 - m32)*C0i(cc2, m12, m32, m22, MD2, MD2, MU2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MD, 2)))/(Pi2*v2) + 
 (0.1875*C0i(cc0, m12, m32, m22, MD2, MD2, MU2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MD, 2))*(m12 + m22 - m32 + 4*(pow(MD, 2))))/(Pi2*v2) + 
 (0.75*A0i(aa0, MS2)*Cos(alp - beta)*Cot(beta)*(pow(MS, 2)))/
  (MA02*Pi2*v2) + (0.375*B0i(bb1, MA02, MS2, MS2)*Cos(alp - beta)*Cot(beta)*
   (pow(MS, 2)))/(Pi2*v2) - (1.5*C0i(cc00, m22, m12, m32, MC2, MS2, MS2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MS, 2)))/(Pi2*v2) - 
 (0.1875*(m12 + 5*m22 - m32)*C0i(cc1, m22, m12, m32, MC2, MS2, MS2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MS, 2)))/(Pi2*v2) + 
 (0.1875*(3*m12 - 3*m22 - m32)*C0i(cc2, m22, m12, m32, MC2, MS2, MS2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MS, 2)))/(Pi2*v2) + 
 (0.75*C0i(cc0, m22, m12, m32, MC2, MS2, MS2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MS, 2))*(-m22 + (pow(MS, 2))))/(Pi2*v2) + 
 (0.75*B0i(bb0, m32, MC2, MS2)*(-1 + Cos(alp)*Csc(beta))*
   ((pow(MC, 2)) + (pow(MS, 2))))/(Pi2*v2) + 
 (0.75*A0i(aa0, MU2)*Cos(alp - beta)*Cot(beta)*(pow(MU, 2)))/
  (MA02*Pi2*v2) + (0.375*B0i(bb1, MA02, MU2, MU2)*Cos(alp - beta)*Cot(beta)*
   (pow(MU, 2)))/(Pi2*v2) - (1.5*C0i(cc00, m22, m12, m32, MD2, MU2, MU2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MU, 2)))/(Pi2*v2) - 
 (0.1875*(m12 + 5*m22 - m32)*C0i(cc1, m22, m12, m32, MD2, MU2, MU2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MU, 2)))/(Pi2*v2) + 
 (0.1875*(3*m12 - 3*m22 - m32)*C0i(cc2, m22, m12, m32, MD2, MU2, MU2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MU, 2)))/(Pi2*v2) + 
 (0.75*C0i(cc0, m22, m12, m32, MD2, MU2, MU2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MU, 2))*(-m22 + (pow(MU, 2))))/(Pi2*v2) + 
 (0.75*B0i(bb0, m32, MD2, MU2)*(-1 + Cos(alp)*Csc(beta))*
   ((pow(MD, 2)) + (pow(MU, 2))))/(Pi2*v2) - 
 (0.015625*(MA02 - Mh2)*B0i(bb0, 0, MA02, Mh2)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*Csc(2*beta)*(pow(Cos(alp - beta), 2)))/
  (MA02*Pi2*v2) - (0.015625*(MA02 - Mh2)*B0i(bb0, MA02, MA02, Mh2)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*Csc(2*beta)*(pow(Cos(alp - beta), 2)))/
  (MA02*Pi2*v2) + (0.015625*A0i(aa0, Mh2)*(2*(M2 - MA02)*Cos(alp - 3*beta) + 
    (-Mh2 + MHH2)*Cos(3*alp - beta) + (2*M2 + 2*MA02 - 3*Mh2 - MHH2)*
     Cos(alp + beta))*Csc(2*beta)*(pow(Cos(alp - beta), 2)))/
  (MA02*Pi2*v2) + (0.0625*C0i(cc00, m22, m12, m32, Mh2, MHp2, MHp2)*
   ((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*Csc(2*beta)*(pow(Cos(alp - beta), 2)))/(Pi2*v2) - 
 1.*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)) + 
 (0.125*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, MW2, MZ2, MZ2)*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.125*MW2*C0i(cc1, m12, m32, m22, MW2, MW2, MZ2)*
   (m12 + (3*m12 + 2*m22 - 2*m32)*Cos(2*w))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.125*C0i(cc00, m12, m32, m22, MW2, MW2, MZ2)*
   (Mh2 + 10*MW2 + 8*MW2*Cos(2*w))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.0625*C0i(cc00, m22, m12, m32, MW2, MZ2, MZ2)*
   (Mh2 + 20*MW2 + (Mh2 + 16*MW2)*Cos(2*w))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.015625*MW2*B0i(bb0, m32, MW2, MZ2)*(23 + 36*Cos(2*w) + 5*Cos(4*w))*
   (pow(Sec(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2)))/(Pi2*v2) - (0.015625*MW2*C0i(cc0, m12, m32, m22, MW2, MW2, MZ2)*
   (-14*m12 + 17*m22 + 13*m32 - 3*Mh2 + 10*MW2 + 
    4*(-4*m12 + 6*m22 + 4*m32 + Mh2 + 4*MW2)*Cos(2*w) + 
    (-2*m12 + 7*m22 + 3*m32 - Mh2 + 6*MW2)*Cos(4*w))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.015625*MW2*C0i(cc0, m22, m12, m32, MW2, MZ2, MZ2)*
   (-15*m12 + 3*m22 + 15*m32 + 10*MW2 + 4*(-5*m12 + m22 + 5*m32 + 6*MW2)*
     Cos(2*w) + (-5*m12 + m22 + 5*m32 - 2*MW2)*Cos(4*w))*(pow(Sec(w), 4))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.125*(Mh2 + 2*MW2)*C0i(cc00, m22, m12, m32, Mh2, MW2, MW2)*
   (1 + (pow(Sin(alp - beta), 3))))/(Pi2*v2) + 
 (0.25*C0i(cc0, m22, m12, m32, Mh2, MW2, MW2)*(pow(MW, 4))*
   (1 + (pow(Sin(alp - beta), 3))))/(Pi2*v2) - 
 (0.125*MW2*(-3*m12 - 2*m22 + 4*m32 - Mh2 + 6*MW2)*
   C0i(cc0, m22, m12, m32, 0, MW2, MW2)*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   (pow(Sin(w), 2)))/(Pi2*v2) + 
 (0.125*(m12 + 9*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, 0, MW2, MW2)*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   (pow(Sin(w), 2)))/(Pi2*v2) + 
 (0.00390625*((Mh2 - MHH2)*Cos(3*alp - 5*beta) + 
    (-8*M2 - 4*MA02 + 5*Mh2 + 3*MHH2)*Cos(alp - 3*beta) + 
    3*Mh2*Cos(3*alp - beta) - 3*MHH2*Cos(3*alp - beta) - 
    8*M2*Cos(alp + beta) + 4*MA02*Cos(alp + beta) + 7*Mh2*Cos(alp + beta) + 
    MHH2*Cos(alp + beta))*Csc(2*beta)*Re(A0i(aa0, MA02)))/(MA02*Pi2*v2) + 
 (0.0078125*((Mh2 - MHH2)*Cos(3*alp - 5*beta) + 
    (-8*M2 - 6*MA02 + 5*Mh2 + 3*MHH2)*Cos(alp - 3*beta) + 
    3*Mh2*Cos(3*alp - beta) - 3*MHH2*Cos(3*alp - beta) - 
    8*M2*Cos(alp + beta) + 6*MA02*Cos(alp + beta) + 7*Mh2*Cos(alp + beta) + 
    MHH2*Cos(alp + beta) + MA02*Cos(alp - 3*beta - 4*w) - 
    MA02*Cos(alp + beta - 4*w) - 8*MA02*Cos(alp - 3*beta - 2*w) + 
    8*MA02*Cos(alp + beta - 2*w) - 8*MA02*Cos(alp - 3*beta + 2*w) + 
    8*MA02*Cos(alp + beta + 2*w) + MA02*Cos(alp - 3*beta + 4*w) - 
    MA02*Cos(alp + beta + 4*w))*Csc(2*beta)*Re(A0i(aa0, MHp2)))/
  (MA02*Pi2*v2) + (0.125*MW2*(-1 + 5*Cos(2*w))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.375*(MT2 - MW2*Cos(2*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MW2, MB2, MT2)))/(Pi2*v2) + 
 (0.375*(-(MW2*Cos(2*w)) + (pow(MU, 2)))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MW2, MD2, MU2)))/(Pi2*v2) - 
 (0.125*MW2*(pow(Csc(w), 2))*(1 + (pow(Sin(alp - beta), 3)))*
   Re(B0i(bb0, MW2, Mh2, MW2)))/(Pi2*v2) + 
 (0.0078125*MW2*(22 + 45*Cos(2*w) + 10*Cos(4*w) + 3*Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.125*MW2*(pow(Csc(w), 2))*(1 + (pow(Sin(alp - beta), 3)))*
   Re(B0i(bb0, MZ2, Mh2, MZ2)))/(Pi2*v2) - 
 (0.5*(-1 + 2*Cos(2*w))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, ME2, ME2)))/(Pi2*v2) - 
 (0.5*(-1 + 2*Cos(2*w))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, ML2, ML2)))/(Pi2*v2) - 
 (0.5*(-1 + 2*Cos(2*w))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, MM2, MM2)))/(Pi2*v2) + 
 (0.125*(pow(Csc(w), 2))*(1 + (pow(Sin(alp - beta), 3)))*
   Re(B0i(bb00, MW2, Mh2, MW2)))/(Pi2*v2) + 
 (0.125*(pow(Csc(w), 2))*(pow(Sin(alp - beta), 3))*
   Re(B0i(bb00, MW2, MHH2, MHp2)))/(Pi2*v2) - 
 (0.125*(pow(Cot(w), 2))*(pow(Sin(alp - beta), 3))*
   Re(B0i(bb00, MZ2, MA02, MHH2)))/(Pi2*v2) + 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, ME2, ME2)))/(Pi2*v2) - 
 (0.125*(pow(Cot(w), 2))*(1 + (pow(Sin(alp - beta), 3)))*
   Re(B0i(bb00, MZ2, Mh2, MZ2)))/(Pi2*v2) + 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, MU2, MU2)))/(Pi2*v2) - 
 (0.0625*(7 + 8*Cos(2*w) + 3*Cos(4*w))*(pow(Cot(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, MW2, MW2)))/(Pi2*v2) - 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (pow(Csc(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(bb1, MZ2, MC2, MC2)))/(Pi2*v2) - 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(Pi2*v2) - 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb1, MZ2, ME2, ME2)))/(Pi2*v2) - 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb1, MZ2, ML2, ML2)))/(Pi2*v2) - 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb1, MZ2, MM2, MM2)))/(Pi2*v2) - 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (pow(Csc(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(bb1, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (pow(Csc(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(bb1, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.015625*(pow(Csc(2*beta), 2))*(pow(Sin(alp - beta), 3))*
   (pow((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta), 2))*
   Re(B0i(dbb0, Mh2, MHH2, MHH2)))/(Pi2*v2) + 
 (0.03125*(-2*Mh2*MW2 + (pow(Mh, 4)) + 12*(pow(MW, 4)))*
   (1 + (pow(Sin(alp - beta), 3)))*Re(B0i(dbb0, Mh2, MW2, MW2)))/
  (Pi2*v2) + (0.75*MW2*(-MT2 + MW2)*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb0, MW2, MB2, MT2)))/(Pi2*v2) + 
 (0.75*MW2*(MW2 - (pow(MU, 2)))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb0, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.25*(pow(MW, 4))*(1 + (pow(Sin(alp - beta), 3)))*
   Re(B0i(dbb0, MW2, Mh2, MW2)))/(Pi2*v2) - 
 (0.0625*(7 + 12*Cos(2*w) + Cos(4*w))*(pow(MW, 4))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb0, MW2, MW2, MZ2)))/(Pi2*v2) - 
 (0.25*MW2*(1 + (pow(Sin(alp - beta), 3)))*Re(B0i(dbb00, MW2, Mh2, MW2)))/
  (Pi2*v2) - (0.25*MW2*(pow(Sin(alp - beta), 3))*
   Re(B0i(dbb00, MW2, MHH2, MHp2)))/(Pi2*v2) + 
 (0.125*Mh2*MW2*(1 + (pow(Sin(alp - beta), 3)))*
   Re(B0i(dbb1, Mh2, MW2, MW2)))/(Pi2*v2) + 
 (0.0625*Mh2*MW2*(pow(Sec(w), 2))*(1 + (pow(Sin(alp - beta), 3)))*
   Re(B0i(dbb1, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (0.09375*C0i(cc00, m12, m32, m22, Mh2, Mh2, MHp2)*
   (M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh2)*Cos(alp + beta))*Csc(beta)*(pow(Cos(alp - beta), 2))*
   Sec(beta))/(Pi2*v2) + (0.03125*C0i(cc00, m22, m12, m32, MHH2, MHp2, MHp2)*
   ((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*Csc(beta)*(pow(Sin(alp - beta), 2))*Sec(beta))/
  (Pi2*v2) + (0.001953125*Re(A0i(aa0, Mh2))*
   (-16 - (2*M2*Cos(3*alp - 5*beta) + (14*M2 + 4*MA02 - 11*Mh2 - MHH2)*
       Cos(alp - 3*beta) - Mh2*Cos(5*alp - 3*beta) + 
      MHH2*Cos(5*alp - 3*beta) + 10*M2*Cos(3*alp - beta) - 
      13*Mh2*Cos(3*alp - beta) + MHH2*Cos(3*alp - beta) + 
      22*M2*Cos(alp + beta) - 4*MA02*Cos(alp + beta) - 
      23*Mh2*Cos(alp + beta) - MHH2*Cos(alp + beta))*Csc(beta)*
     (pow(MA02, -1))*Sec(beta)))/(Pi2*v2) + 
 (0.09375*C0i(cc00, m12, m32, m22, Mh2, Mh2, MW2)*
   (-4*Mh2 - (M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
      (2*M2 - 3*Mh2)*Cos(alp + beta))*Csc(beta)*(pow(Sin(alp - beta), 2))*
     Sec(beta)))/(Pi2*v2) + 
 (0.09375*MW2*C0i(cc0, m12, m32, m22, Mh2, Mh2, MW2)*
   (4*Mh2 + (M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
      (2*M2 - 3*Mh2)*Cos(alp + beta))*Csc(beta)*(pow(Sin(alp - beta), 2))*
     Sec(beta)))/(Pi2*v2) - (0.25*B0i(bb0, m32, 0, ME2)*(pow(ME, 2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.5*C0i(cc00, m22, m12, m32, 0, ME2, ME2)*(pow(ME, 2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.0625*(m12 + 5*m22 - m32)*C0i(cc1, m22, m12, m32, 0, ME2, ME2)*
   (pow(ME, 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*(3*m12 - 3*m22 - m32)*C0i(cc2, m22, m12, m32, 0, ME2, ME2)*
   (pow(ME, 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.25*C0i(cc0, m22, m12, m32, 0, ME2, ME2)*(pow(ME, 2))*
   (-m22 + (pow(ME, 2)))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.25*B0i(bb0, m32, 0, ML2)*(pow(ML, 2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.5*C0i(cc00, m22, m12, m32, 0, ML2, ML2)*(pow(ML, 2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.0625*(m12 + 5*m22 - m32)*C0i(cc1, m22, m12, m32, 0, ML2, ML2)*
   (pow(ML, 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*(3*m12 - 3*m22 - m32)*C0i(cc2, m22, m12, m32, 0, ML2, ML2)*
   (pow(ML, 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.25*C0i(cc0, m22, m12, m32, 0, ML2, ML2)*(pow(ML, 2))*
   (-m22 + (pow(ML, 2)))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.25*B0i(bb0, m32, 0, MM2)*(pow(MM, 2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.5*C0i(cc00, m22, m12, m32, 0, MM2, MM2)*(pow(MM, 2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.0625*(m12 + 5*m22 - m32)*C0i(cc1, m22, m12, m32, 0, MM2, MM2)*
   (pow(MM, 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*(3*m12 - 3*m22 - m32)*C0i(cc2, m22, m12, m32, 0, MM2, MM2)*
   (pow(MM, 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.25*C0i(cc0, m22, m12, m32, 0, MM2, MM2)*(pow(MM, 2))*
   (-m22 + (pow(MM, 2)))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.75*Cos(alp)*Csc(beta)*(pow(MB, 4))*Re(B0i(bb0, Mh2, MB2, MB2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*Csc(beta)*(pow(MC, 4))*Re(B0i(bb0, Mh2, MC2, MC2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*Csc(beta)*(pow(MD, 4))*Re(B0i(bb0, Mh2, MD2, MD2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*Csc(beta)*(pow(MS, 4))*Re(B0i(bb0, Mh2, MS2, MS2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*Csc(beta)*(pow(MT, 4))*Re(B0i(bb0, Mh2, MT2, MT2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*Csc(beta)*(pow(MU, 4))*Re(B0i(bb0, Mh2, MU2, MU2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*Csc(beta)*(pow(MB, 4))*Re(B0i(bb0, MHH2, MB2, MB2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*Csc(beta)*(pow(MC, 4))*Re(B0i(bb0, MHH2, MC2, MC2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*Csc(beta)*(pow(MD, 4))*Re(B0i(bb0, MHH2, MD2, MD2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*Csc(beta)*(pow(MS, 4))*Re(B0i(bb0, MHH2, MS2, MS2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*Csc(beta)*(pow(MT, 4))*Re(B0i(bb0, MHH2, MT2, MT2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*Csc(beta)*(pow(MU, 4))*Re(B0i(bb0, MHH2, MU2, MU2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*Csc(beta)*(pow(MB, 2))*Re(B0i(bb1, MHH2, MB2, MB2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*Csc(beta)*(pow(MC, 2))*Re(B0i(bb1, MHH2, MC2, MC2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*Csc(beta)*(pow(MD, 2))*Re(B0i(bb1, MHH2, MD2, MD2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*Csc(beta)*(pow(MS, 2))*Re(B0i(bb1, MHH2, MS2, MS2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*MT2*Cos(alp)*Csc(beta)*Re(B0i(bb1, MHH2, MT2, MT2))*Sin(alp)*
   (Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*Csc(beta)*(pow(MU, 2))*Re(B0i(bb1, MHH2, MU2, MU2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*(pow(ME, 2))*Re(B0i(bb1, Mh2, ME2, ME2))*
   (2*(Mh2 - MHH2) + MHH2*Cos(2*alp - beta)*(pow(Sec(beta), 2))*
     Sin(alp) + (2*Mh2 - MHH2)*Sec(beta)*Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*(pow(ML, 2))*Re(B0i(bb1, Mh2, ML2, ML2))*
   (2*(Mh2 - MHH2) + MHH2*Cos(2*alp - beta)*(pow(Sec(beta), 2))*
     Sin(alp) + (2*Mh2 - MHH2)*Sec(beta)*Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*(pow(MM, 2))*Re(B0i(bb1, Mh2, MM2, MM2))*
   (2*(Mh2 - MHH2) + MHH2*Cos(2*alp - beta)*(pow(Sec(beta), 2))*
     Sin(alp) + (2*Mh2 - MHH2)*Sec(beta)*Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0625*(Mh2 - MHH2)*A0i(aa0, MW2)*(pow(Cos(alp - beta), 2))*
   Sin(alp - beta))/(MA02*Pi2*v2) + 
 (0.09375*(Mh2 - MHH2)*A0i(aa0, MZ2)*(pow(Cos(alp - beta), 2))*
   Sin(alp - beta))/(MA02*Pi2*v2) - 
 (0.125*(Mh2 - MHp2 + MW2)*C0i(cc00, m22, m12, m32, Mh2, MHp2, MW2)*
   (pow(Cos(alp - beta), 2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.125*(Mh2 - MHp2 + MW2)*C0i(cc00, m22, m12, m32, MHH2, MHp2, MW2)*
   (pow(Cos(alp - beta), 2))*Sin(alp - beta))/(Pi2*v2) - 
 (0.125*(Mh2 + 2*MW2)*C0i(cc00, m22, m12, m32, MHH2, MW2, MW2)*
   (pow(Cos(alp - beta), 2))*Sin(alp - beta))/(Pi2*v2) - 
 (0.125*(Mh2 - MHp2 + MW2)*C0i(cc00, m32, m12, m22, Mh2, MHp2, MW2)*
   (pow(Cos(alp - beta), 2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.125*(Mh2 - MHp2 + MW2)*C0i(cc00, m32, m12, m22, MHH2, MHp2, MW2)*
   (pow(Cos(alp - beta), 2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.25*C0i(cc0, m22, m12, m32, MHH2, MW2, MW2)*(pow(MW, 4))*
   (pow(Cos(alp - beta), 2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.0625*MW2*B0i(bb1, MA02, Mh2, MZ2)*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Sin(alp - beta))/(Pi2*v2) - 
 (0.0625*MW2*B0i(bb1, MA02, MHH2, MZ2)*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.015625*Mh2*B0i(bb0, 0, Mh2, MZ2)*(-MA02 + Mh2 - 2*MW2 + 
    (-MA02 + Mh2)*Cos(2*w))*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Sin(alp - beta))/(MA02*Pi2*v2) + 
 (0.015625*MHH2*B0i(bb0, 0, MHH2, MZ2)*(MA02 - MHH2 + 2*MW2 + 
    (MA02 - MHH2)*Cos(2*w))*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Sin(alp - beta))/(MA02*Pi2*v2) + 
 (0.015625*B0i(bb0, MA02, MHH2, MZ2)*(-(MHH2*(MHH2 - 2*MW2)) + 
    MA02*(MHH2 + 2*MW2) + (MA02 - MHH2)*MHH2*Cos(2*w))*
   (pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*Sin(alp - beta))/
  (MA02*Pi2*v2) - (0.015625*(MHH2*(Mh2 - 2*MW2) - MA02*(Mh2 + MHH2 + 2*MW2) + 
    (MA02 - Mh2)*(MA02 - MHH2)*Cos(2*w) + (pow(MA02, 2)))*
   (pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*
   Re(B0i(bb0, Mh2, MA02, MZ2))*Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*(Mh2*(MHH2 - MHp2) + MHp2*(MHp2 - MW2) - MHH2*(MHp2 + MW2))*
   (pow(Cos(alp - beta), 2))*Re(B0i(bb0, Mh2, MHp2, MW2))*
   Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0625*(Mh2*(MHH2 - MHp2) + MHp2*(MHp2 - MW2) - MHH2*(MHp2 + MW2))*
   (pow(Cos(alp - beta), 2))*Re(B0i(bb0, MHH2, MHp2, MW2))*
   Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.03125*(Mh2*MHH2 - 2*MHH2*MW2 + 12*(pow(MW, 4)))*
   (pow(Cos(alp - beta), 2))*Re(B0i(bb0, MHH2, MW2, MW2))*
   Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.001953125*(3*Mh2*MHH2 - 8*MHH2*MW2 + 4*MHH2*(Mh2 - 2*MW2)*Cos(2*w) + 
    Mh2*MHH2*Cos(4*w) + 96*(pow(MW, 4)))*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 4))*Re(B0i(bb0, MHH2, MZ2, MZ2))*Sin(alp - beta))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.125*MW2*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb0, MW2, MHH2, MW2))*Sin(alp - beta))/
  (Pi2*v2) + (0.125*MW2*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Re(B0i(bb0, MZ2, MHH2, MZ2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.5*Cos(2*w)*(pow(Sin(w), 2))*Re(B0i(bb00, 0, MHp2, MHp2))*
   Sin(alp - beta))/(Pi2*v2) + 
 (0.125*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, MA02, MHp2))*Sin(alp - beta))/
  (Pi2*v2) + (0.125*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MW2, Mh2, MHp2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.125*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MW2, MHH2, MW2))*Sin(alp - beta))/(Pi2*v2) - 
 (0.125*(pow(Cos(alp - beta), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MA02, Mh2))*Sin(alp - beta))/(Pi2*v2) - 
 (0.125*(pow(Cos(alp - beta), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MHH2, MZ2))*Sin(alp - beta))/(Pi2*v2) - 
 (0.0078125*(pow(Csc(w), 4))*(pow(Sin(4*w), 2))*
   Re(B0i(bb00, MZ2, MHp2, MHp2))*Sin(alp - beta))/(Pi2*v2) - 
 (0.0625*MHH2*MW2*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*
   Re(B0i(bb1, Mh2, MA02, MZ2))*Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.125*MHH2*MW2*(pow(Cos(alp - beta), 2))*Re(B0i(bb1, Mh2, MHp2, MW2))*
   Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0625*MHH2*MW2*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*
   Re(B0i(bb1, MHH2, MA02, MZ2))*Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.125*MHH2*MW2*(pow(Cos(alp - beta), 2))*Re(B0i(bb1, MHH2, MHp2, MW2))*
   Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.125*MHH2*MW2*(pow(Cos(alp - beta), 2))*Re(B0i(bb1, MHH2, MW2, MW2))*
   Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*MHH2*MW2*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*
   Re(B0i(bb1, MHH2, MZ2, MZ2))*Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.00390625*(pow((2*MA02 - Mh2)*Cos(alp - 3*beta) + 
      (4*M2 - 2*MA02 - 3*Mh2)*Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh2, MA02, MA02))*Sin(alp - beta))/(Pi2*v2) + 
 (0.015625*(-2*Mh2*MW2 - 2*MA02*(Mh2 + MW2) + (pow(MA02, 2)) + 
    (pow(Mh, 4)) + Cos(2*w)*(pow(MA02 - Mh2, 2)))*
   (pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, Mh2, MA02, MZ2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.03125*(pow(Cos(alp - beta), 2))*(pow(Csc(2*beta), 2))*
   (pow((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta), 2))*
   Re(B0i(dbb0, Mh2, Mh2, MHH2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.0078125*(pow((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + 
      (-4*M2 + 3*Mh2 + 2*MHp2)*Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh2, MHp2, MHp2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.0625*(MHp2*(MHp2 - MW2) - Mh2*(2*MHp2 + MW2) + (pow(Mh, 4)))*
   (pow(Cos(alp - beta), 2))*Re(B0i(dbb0, Mh2, MHp2, MW2))*
   Sin(alp - beta))/(Pi2*v2) + 
 (0.25*(pow(MW, 4))*(pow(Cos(alp - beta), 2))*
   Re(B0i(dbb0, MW2, MHH2, MW2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MHp2, MHp2))*Sin(alp - beta))/
  (Pi2*v2) - (0.25*MW2*Re(B0i(dbb00, MW2, MA02, MHp2))*Sin(alp - beta))/
  (Pi2*v2) - (0.25*MW2*(pow(Cos(alp - beta), 2))*
   Re(B0i(dbb00, MW2, Mh2, MHp2))*Sin(alp - beta))/(Pi2*v2) - 
 (0.25*MW2*(pow(Cos(alp - beta), 2))*Re(B0i(dbb00, MW2, MHH2, MW2))*
   Sin(alp - beta))/(Pi2*v2) + (0.0625*Mh2*MW2*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Re(B0i(dbb1, Mh2, MA02, MZ2))*Sin(alp - beta))/
  (Pi2*v2) + (0.125*Mh2*MW2*(pow(Cos(alp - beta), 2))*
   Re(B0i(dbb1, Mh2, MHp2, MW2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.0625*(Mh2 + 12*MW2)*B0i(bb0, m12, MW2, MW2)*(1 + Sin(alp - beta)))/
  (Pi2*v2) + (0.03125*(Mh2 + 24*MW2)*B0i(bb0, m12, MZ2, MZ2)*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.125*MW2*B0i(bb0, m22, Mh2, MW2)*(1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.125*MW2*B0i(bb0, m32, Mh2, MW2)*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.125*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, MW2, MZ2, MZ2)*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.0625*MW2*C0i(cc2, m12, m32, m22, MW2, MW2, MZ2)*
   (m12 + m22 - m32 + (m12 + 9*m22 - m32)*Cos(2*w))*(1 + Sin(alp - beta)))/
  (Pi2*v2) + (0.125*MW2*B0i(bb0, m22, 0, MW2)*(pow(Sin(w), 2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.625*MW2*B0i(bb0, m32, 0, MW2)*(pow(Sin(w), 2))*(1 + Sin(alp - beta)))/
  (Pi2*v2) - (2.*MW2*C0i(cc00, m22, m12, m32, 0, MW2, MW2)*
   (pow(Sin(w), 2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.125*(5*m12 - 5*m22 - 3*m32)*MW2*C0i(cc2, m22, m12, m32, 0, MW2, MW2)*
   (pow(Sin(w), 2))*(1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.125*MW2*B0i(bb0, m22, MW2, MZ2)*(pow(Sin(w), 2))*(pow(Tan(w), 2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (1.125*MW2*(pow(Sin(w), 2))*Re(B0i(bb0, 0, MW2, MW2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.375*(pow(MC, 2))*(pow(Csc(w), 2))*Re(B0i(bb0, MW2, MC2, MS2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.1875*(pow(MB, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MB2, MB2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.1875*(pow(MC, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MC2, MC2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.1875*(pow(MD, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MD2, MD2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.0625*(pow(ME, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, ME2, ME2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.0625*(pow(ML, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, ML2, ML2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.0625*(pow(MM, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MM2, MM2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.1875*(pow(MS, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MS2, MS2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.1875*MT2*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MT2, MT2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.1875*(pow(MU, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MU2, MU2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.0625*MW2*(5 + 9*Cos(2*w))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MW2, MW2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MB2, MB2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MC2, MC2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MD2, MD2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MS2, MS2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MT2, MT2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MU2, MU2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.25*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, 0, ME2))*(1 + Sin(alp - beta)))/
  (Pi2*v2) - (0.25*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, 0, ML2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.25*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, 0, MM2))*(1 + Sin(alp - beta)))/
  (Pi2*v2) + (1.*Re(B0i(bb00, MW2, 0, MW2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.75*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, MB2, MT2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.75*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, MC2, MS2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.75*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, MD2, MU2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.125*(5 + 4*Cos(2*w))*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, MW2, MZ2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.375*(pow(Cot(w), 2))*Re(B0i(bb00, MZ2, 0, 0))*(1 + Sin(alp - beta)))/
  (Pi2*v2) + (0.16666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MB2, MB2))*(1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MC2, MC2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MD2, MD2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ME2, ME2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ML2, ML2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MM2, MM2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MS2, MS2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MT2, MT2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MU2, MU2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.25*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MW2, MW2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.125*MW2*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb1, MW2, 0, ME2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.125*MW2*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb1, MW2, 0, ML2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.125*MW2*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb1, MW2, 0, MM2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.25*MW2*Cos(2*w)*Re(B0i(bb1, MW2, 0, MW2))*(1 + Sin(alp - beta)))/
  (Pi2*v2) - (0.375*MW2*Cos(2*w)*(pow(Csc(w), 2))*
   Re(B0i(bb1, MW2, MB2, MT2))*(1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.375*MW2*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb1, MW2, MC2, MS2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.375*MW2*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb1, MW2, MD2, MU2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.25*MW2*Cos(2*w)*(pow(Cot(w), 2))*Re(B0i(bb1, MW2, MW2, MZ2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.1875*MW2*(pow(Csc(w), 2))*Re(B0i(bb1, MZ2, 0, 0))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.25*MW2*(pow(Cos(w), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb1, MZ2, MW2, MW2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (1.*(pow(MW, 4))*(pow(Sin(w), 2))*Re(B0i(dbb0, MW2, 0, MW2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.75*MW2*(pow(MC, 2))*Re(B0i(dbb0, MW2, MC2, MS2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MB2, MB2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MC2, MC2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MD2, MD2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ME2, ME2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ML2, ML2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MM2, MM2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MS2, MS2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MT2, MT2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MU2, MU2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (1.5*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MW2, MW2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.5*MW2*Re(B0i(dbb00, MW2, 0, ME2))*(1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.5*MW2*Re(B0i(dbb00, MW2, 0, ML2))*(1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.5*MW2*Re(B0i(dbb00, MW2, 0, MM2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (2.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, MW2, 0, MW2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (1.5*MW2*Re(B0i(dbb00, MW2, MB2, MT2))*(1 + Sin(alp - beta)))/(Pi2*v2) + 
 (1.5*MW2*Re(B0i(dbb00, MW2, MC2, MS2))*(1 + Sin(alp - beta)))/(Pi2*v2) + 
 (1.5*MW2*Re(B0i(dbb00, MW2, MD2, MU2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.25*MW2*(5 + 4*Cos(2*w))*Re(B0i(dbb00, MW2, MW2, MZ2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.25*(pow(MW, 4))*Re(B0i(dbb1, MW2, 0, ME2))*(1 + Sin(alp - beta)))/
  (Pi2*v2) - (0.25*(pow(MW, 4))*Re(B0i(dbb1, MW2, 0, ML2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.25*(pow(MW, 4))*Re(B0i(dbb1, MW2, 0, MM2))*(1 + Sin(alp - beta)))/
  (Pi2*v2) - (0.5*(pow(MW, 4))*(pow(Sin(w), 2))*
   Re(B0i(dbb1, MW2, 0, MW2))*(1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.75*(pow(MW, 4))*Re(B0i(dbb1, MW2, MB2, MT2))*(1 + Sin(alp - beta)))/
  (Pi2*v2) - (0.75*(pow(MW, 4))*Re(B0i(dbb1, MW2, MC2, MS2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.75*(pow(MW, 4))*Re(B0i(dbb1, MW2, MD2, MU2))*(1 + Sin(alp - beta)))/
  (Pi2*v2) + (0.5*(pow(MW, 4))*(pow(Cos(w), 2))*
   Re(B0i(dbb1, MW2, MW2, MZ2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.0009765625*(-8*Mh2*MW2 + 4*Mh2*(Mh2 - 2*MW2)*Cos(2*w) + 
    3*(pow(Mh, 4)) + Cos(4*w)*(pow(Mh, 4)) + 96*(pow(MW, 4)))*
   (pow(Sec(w), 4))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(dbb0, Mh2, MZ2, MZ2))*(-3 + Cos(2*(alp - beta)) + 
    2*Sin(alp - beta)))/(Pi2*v2) + 
 (0.0625*MW2*Re(B0i(bb1, Mh2, MW2, MW2))*(2*(Mh2 - MHH2) + 
    (2*Mh2 - MHH2 + MHH2*Cos(2*(alp - beta)))*Sin(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.03125*MW2*(pow(Sec(w), 2))*
   Re(B0i(bb1, Mh2, MZ2, MZ2))*(2*(Mh2 - MHH2) + 
    (2*Mh2 - MHH2 + MHH2*Cos(2*(alp - beta)))*Sin(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.0078125*(pow(Csc(w), 2))*Re(A0i(aa0, MW2))*
   (-2 + 23*Cos(2*w) + 30*Cos(4*w) - 3*Cos(6*w) + 
    (-2*MA02 - 2*Mh2 + 2*MHH2 - 2*(Mh2 - MHH2)*Cos(2*(alp - beta)) + 
      (Mh2 - MHH2)*Cos(2*(alp - beta - w)) + 23*MA02*Cos(2*w) + 
      2*Mh2*Cos(2*w) - 2*MHH2*Cos(2*w) + 30*MA02*Cos(4*w) - 3*MA02*Cos(6*w) + 
      Mh2*Cos(2*(alp - beta + w)) - MHH2*Cos(2*(alp - beta + w)))*
     (pow(MA02, -1))*Sin(alp - beta)))/(Pi2*v2) + 
 (0.015625*Re(B0i(bb0, Mh2, MW2, MW2))*(4*(-Mh2 + MHH2)*MW2 + 
    (Mh2*(MHH2 - 4*MW2) + 2*MW2*(MHH2 + 6*MW2) + Cos(2*(alp - beta))*
       (Mh2*MHH2 - 2*MHH2*MW2 + 12*(pow(MW, 4))))*Sin(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.75*(pow(MB, 4))*Re(B0i(dbb0, Mh2, MB2, MB2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.75*(pow(MC, 4))*Re(B0i(dbb0, Mh2, MC2, MC2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.75*(pow(MD, 4))*Re(B0i(dbb0, Mh2, MD2, MD2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.75*(pow(MS, 4))*Re(B0i(dbb0, Mh2, MS2, MS2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.75*(pow(MT, 4))*Re(B0i(dbb0, Mh2, MT2, MT2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.75*(pow(MU, 4))*Re(B0i(dbb0, Mh2, MU2, MU2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.375*Mh2*(pow(MB, 2))*Re(B0i(dbb1, Mh2, MB2, MB2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.375*Mh2*(pow(MC, 2))*Re(B0i(dbb1, Mh2, MC2, MC2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.375*Mh2*(pow(MD, 2))*Re(B0i(dbb1, Mh2, MD2, MD2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.375*Mh2*(pow(MS, 2))*Re(B0i(dbb1, Mh2, MS2, MS2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.375*Mh2*MT2*Re(B0i(dbb1, Mh2, MT2, MT2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.375*Mh2*(pow(MU, 2))*Re(B0i(dbb1, Mh2, MU2, MU2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) + (0.03515625*Re(B0i(dbb0, Mh2, Mh2, Mh2))*
   (4*(pow(Mh, 4)) + 
    (pow(M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
        (2*M2 - 3*Mh2)*Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
     Sin(alp - beta)))/(Pi2*v2) + (0.00048828125*Re(B0i(bb0, Mh2, MZ2, MZ2))*
   (-256*MW2*(pow(Csc(2*w), 2)) + (pow(Mh2 - MHH2, -1))*
     (6*Mh2*MHH2 - 32*Mh2*MW2 + 16*MHH2*MW2 + 
      Mh2*MHH2*Cos(2*(alp - beta - 2*w)) + 4*Mh2*MHH2*
       Cos(2*(alp - beta - w)) - 8*MHH2*MW2*Cos(2*(alp - beta - w)) + 
      8*Mh2*MHH2*Cos(2*w) - 32*Mh2*MW2*Cos(2*w) + 16*MHH2*MW2*Cos(2*w) + 
      2*Mh2*MHH2*Cos(4*w) + 4*Mh2*MHH2*Cos(2*(alp - beta + w)) - 
      8*MHH2*MW2*Cos(2*(alp - beta + w)) + 
      Mh2*MHH2*Cos(2*(alp - beta + 2*w)) + 192*(pow(MW, 4)) + 
      2*Cos(2*(alp - beta))*(3*Mh2*MHH2 - 8*MHH2*MW2 + 96*(pow(MW, 4))))*
     (pow(Csc(w), 2))*(pow(Sec(w), 4))*Sin(alp - beta)))/
  (Pi2*v2*(pow(Csc(w), 2))) - 
 (0.25*(pow(ME, 4))*Re(B0i(dbb0, Mh2, ME2, ME2))*
   (1 + (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.25*(pow(ML, 4))*Re(B0i(dbb0, Mh2, ML2, ML2))*
   (1 + (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.25*(pow(MM, 4))*Re(B0i(dbb0, Mh2, MM2, MM2))*
   (1 + (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.125*Mh2*(pow(ME, 2))*Re(B0i(dbb1, Mh2, ME2, ME2))*
   (1 + (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.125*Mh2*(pow(ML, 2))*Re(B0i(dbb1, Mh2, ML2, ML2))*
   (1 + (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.125*Mh2*(pow(MM, 2))*Re(B0i(dbb1, Mh2, MM2, MM2))*
   (1 + (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.00390625*(pow(Csc(w), 2))*Re(A0i(aa0, MZ2))*
   (52*MA02 + (52*MA02 + 2*Mh2 - 2*MHH2 + 2*(Mh2 - MHH2)*
       Cos(2*(alp - beta)) + (-Mh2 + MHH2)*Cos(2*(alp - beta - w)) - 
      Mh2*Cos(2*(alp - beta + w)) + MHH2*Cos(2*(alp - beta + w)))*
     Sin(alp - beta) + 2*Cos(2*w)*(22*MA02 + (22*MA02 - Mh2 + MHH2)*
       Sin(alp - beta))))/(MA02*Pi2*v2) - 
 (0.015625*B0i(bb0, MA02, Mh2, MZ2)*Cos(alp - beta)*
   ((MA02 - Mh2)*Mh2 + (MA02 + Mh2)*MW2*(pow(Sec(w), 2)))*
   Sin(2*(alp - beta)))/(MA02*Pi2*v2) + 
 (0.015625*Cos(alp - beta)*((MA02 - Mh2)*(MA02 - MHH2) - 
    (MA02 + MHH2)*MW2*(pow(Sec(w), 2)))*Re(B0i(bb0, MHH2, MA02, MZ2))*
   Sin(2*(alp - beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MB, 2))*Re(B0i(bb1, Mh2, MB2, MB2))*
   (-2*Mh2 + 2*MHH2 + Cos(alp)*Csc(beta)*(2*Mh2 - MHH2 + 
      MHH2*Csc(beta)*Sin(2*alp - beta))))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MC, 2))*Re(B0i(bb1, Mh2, MC2, MC2))*
   (-2*Mh2 + 2*MHH2 + Cos(alp)*Csc(beta)*(2*Mh2 - MHH2 + 
      MHH2*Csc(beta)*Sin(2*alp - beta))))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MD, 2))*Re(B0i(bb1, Mh2, MD2, MD2))*
   (-2*Mh2 + 2*MHH2 + Cos(alp)*Csc(beta)*(2*Mh2 - MHH2 + 
      MHH2*Csc(beta)*Sin(2*alp - beta))))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MS, 2))*Re(B0i(bb1, Mh2, MS2, MS2))*
   (-2*Mh2 + 2*MHH2 + Cos(alp)*Csc(beta)*(2*Mh2 - MHH2 + 
      MHH2*Csc(beta)*Sin(2*alp - beta))))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*MT2*Re(B0i(bb1, Mh2, MT2, MT2))*(-2*Mh2 + 2*MHH2 + 
    Cos(alp)*Csc(beta)*(2*Mh2 - MHH2 + MHH2*Csc(beta)*Sin(2*alp - beta))))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.1875*(pow(MU, 2))*
   Re(B0i(bb1, Mh2, MU2, MU2))*(-2*Mh2 + 2*MHH2 + 
    Cos(alp)*Csc(beta)*(2*Mh2 - MHH2 + MHH2*Csc(beta)*Sin(2*alp - beta))))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.010416666666666666*Re(A0i(aa0, MT2))*
   (-29*MA02 - 72*MT2*Cos(beta)*Cot(beta)*(pow(Cos(alp), 3)) - 
    29*MA02*(pow(Csc(w), 2)) - 20*MA02*Cos(4*w)*(pow(Csc(w), 2)) + 
    2*MA02*Cos(6*w)*(pow(Csc(w), 2)) - 72*MT2*Cos(beta)*
     (pow(Sin(alp), 3)) - 72*MT2*Cos(beta)*(pow(Cos(alp), 2))*
     Sin(alp) + 48*MA02*Cos(beta)*(pow(Sin(w), 2))*Sin(alp) - 
    128*MA02*Cos(beta)*(pow(Sin(w), 4))*Sin(alp) + 
    MA02*(pow(Cot(w), 2))*(29 - 2*Cos(beta)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*
       Sin(alp)) + Cos(alp)*(-72*MT2*Cos(beta)*Cot(beta)*
       (pow(Sin(alp), 2)) + MA02*(29 - 29*Cos(2*w) + 20*Cos(4*w) - 
        2*Cos(6*w))*(pow(Csc(w), 2))*Sin(beta))))/(MA02*Pi2*v2) + 
 (0.010416666666666666*Re(A0i(aa0, MC2))*
   (-29*MA02 - 72*Cos(beta)*Cot(beta)*(pow(MC, 2))*
     (pow(Cos(alp), 3)) - 29*MA02*(pow(Csc(w), 2)) - 
    20*MA02*Cos(4*w)*(pow(Csc(w), 2)) + 2*MA02*Cos(6*w)*
     (pow(Csc(w), 2)) - 72*Cos(beta)*(pow(MC, 2))*
     (pow(Sin(alp), 3)) - 72*Cos(beta)*(pow(MC, 2))*
     (pow(Cos(alp), 2))*Sin(alp) + 48*MA02*Cos(beta)*(pow(Sin(w), 2))*
     Sin(alp) - 128*MA02*Cos(beta)*(pow(Sin(w), 4))*Sin(alp) + 
    MA02*(pow(Cot(w), 2))*(29 - 2*Cos(beta)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*
       Sin(alp)) + Cos(alp)*(-72*Cos(beta)*Cot(beta)*(pow(MC, 2))*
       (pow(Sin(alp), 2)) + MA02*(29 - 29*Cos(2*w) + 20*Cos(4*w) - 
        2*Cos(6*w))*(pow(Csc(w), 2))*Sin(beta))))/(MA02*Pi2*v2) + 
 (0.010416666666666666*Re(A0i(aa0, MU2))*
   (-29*MA02 - 72*Cos(beta)*Cot(beta)*(pow(MU, 2))*
     (pow(Cos(alp), 3)) - 29*MA02*(pow(Csc(w), 2)) - 
    20*MA02*Cos(4*w)*(pow(Csc(w), 2)) + 2*MA02*Cos(6*w)*
     (pow(Csc(w), 2)) - 72*Cos(beta)*(pow(MU, 2))*
     (pow(Sin(alp), 3)) - 72*Cos(beta)*(pow(MU, 2))*
     (pow(Cos(alp), 2))*Sin(alp) + 48*MA02*Cos(beta)*(pow(Sin(w), 2))*
     Sin(alp) - 128*MA02*Cos(beta)*(pow(Sin(w), 4))*Sin(alp) + 
    MA02*(pow(Cot(w), 2))*(29 - 2*Cos(beta)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*
       Sin(alp)) + Cos(alp)*(-72*Cos(beta)*Cot(beta)*(pow(MU, 2))*
       (pow(Sin(alp), 2)) + MA02*(29 - 29*Cos(2*w) + 20*Cos(4*w) - 
        2*Cos(6*w))*(pow(Csc(w), 2))*Sin(beta))))/(MA02*Pi2*v2) + 
 (0.005208333333333333*Re(A0i(aa0, MB2))*
   (-144*Cos(beta)*Cot(beta)*(pow(MB, 2))*(pow(Cos(alp), 3))*
     (pow(Csc(w), 4)) + MA02*(56 - 11*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
     (pow(Csc(w), 6)) - 64*MA02*Cos(beta)*Sin(alp) + 
    48*MA02*Cos(beta)*(pow(Csc(w), 2))*Sin(alp) - 
    144*Cos(beta)*(pow(MB, 2))*(pow(Cos(alp), 2))*
     (pow(Csc(w), 4))*Sin(alp) - 4*Cos(beta)*(pow(Csc(w), 4))*
     (MA02*(-12 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2)) - 
      18*(MA02 - 2*(pow(MB, 2))*(pow(Sin(alp), 2))))*Sin(alp) - 
    Cos(alp)*(144*Cos(beta)*Cot(beta)*(pow(MB, 2))*(pow(Csc(w), 4))*
       (pow(Sin(alp), 2)) + MA02*(56 - 11*Cos(2*w) - 10*Cos(4*w) + 
        Cos(6*w))*(pow(Csc(w), 6))*Sin(beta))))/
  (MA02*Pi2*v2*(pow(Csc(w), 4))) + 
 (0.005208333333333333*Re(A0i(aa0, MD2))*
   (-144*Cos(beta)*Cot(beta)*(pow(MD, 2))*(pow(Cos(alp), 3))*
     (pow(Csc(w), 4)) + MA02*(56 - 11*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
     (pow(Csc(w), 6)) - 64*MA02*Cos(beta)*Sin(alp) + 
    48*MA02*Cos(beta)*(pow(Csc(w), 2))*Sin(alp) - 
    144*Cos(beta)*(pow(MD, 2))*(pow(Cos(alp), 2))*
     (pow(Csc(w), 4))*Sin(alp) - 4*Cos(beta)*(pow(Csc(w), 4))*
     (MA02*(-12 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2)) - 
      18*(MA02 - 2*(pow(MD, 2))*(pow(Sin(alp), 2))))*Sin(alp) - 
    Cos(alp)*(144*Cos(beta)*Cot(beta)*(pow(MD, 2))*(pow(Csc(w), 4))*
       (pow(Sin(alp), 2)) + MA02*(56 - 11*Cos(2*w) - 10*Cos(4*w) + 
        Cos(6*w))*(pow(Csc(w), 6))*Sin(beta))))/
  (MA02*Pi2*v2*(pow(Csc(w), 4))) + 
 (0.005208333333333333*Re(A0i(aa0, MS2))*
   (-144*Cos(beta)*Cot(beta)*(pow(MS, 2))*(pow(Cos(alp), 3))*
     (pow(Csc(w), 4)) + MA02*(56 - 11*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
     (pow(Csc(w), 6)) - 64*MA02*Cos(beta)*Sin(alp) + 
    48*MA02*Cos(beta)*(pow(Csc(w), 2))*Sin(alp) - 
    144*Cos(beta)*(pow(MS, 2))*(pow(Cos(alp), 2))*
     (pow(Csc(w), 4))*Sin(alp) - 4*Cos(beta)*(pow(Csc(w), 4))*
     (MA02*(-12 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2)) - 
      18*(MA02 - 2*(pow(MS, 2))*(pow(Sin(alp), 2))))*Sin(alp) - 
    Cos(alp)*(144*Cos(beta)*Cot(beta)*(pow(MS, 2))*(pow(Csc(w), 4))*
       (pow(Sin(alp), 2)) + MA02*(56 - 11*Cos(2*w) - 10*Cos(4*w) + 
        Cos(6*w))*(pow(Csc(w), 6))*Sin(beta))))/
  (MA02*Pi2*v2*(pow(Csc(w), 4))) + 
 (0.09375*B0i(bb0, m12, Mh2, Mh2)*Csc(2*beta)*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   (2*(M2 - Mh2)*Cos(alp + beta) + (-M2 + Mh2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (Pi2*v2) - (0.125*C0i(cc00, m12, m32, m22, MHH2, MHH2, MHp2)*Csc(2*beta)*
   (pow(Sin(alp - beta), 3))*((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(Pi2*v2) + (0.03125*B0i(bb0, m12, MHH2, MHH2)*
   Csc(2*beta)*Sin(alp - beta)*((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(Pi2*v2) + 
 (0.0625*MW2*C0i(cc0, m12, m32, m22, MHH2, MHH2, MW2)*Csc(beta)*
   (pow(Cos(alp - beta), 2))*Sec(beta)*Sin(alp - beta)*
   ((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) - 
 (0.0625*C0i(cc00, m12, m32, m22, MHH2, MHH2, MW2)*Csc(beta)*
   (pow(Cos(alp - beta), 2))*Sec(beta)*Sin(alp - beta)*
   ((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) + 
 (0.0234375*(M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh2)*Cos(alp + beta))*(pow(Cos(alp - beta), 2))*
   (pow(Csc(2*beta), 2))*Re(B0i(bb0, Mh2, Mh2, Mh2))*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - 
 (0.0234375*(M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh2)*Cos(alp + beta))*(pow(Cos(alp - beta), 2))*
   (pow(Csc(2*beta), 2))*Re(B0i(bb0, MHH2, Mh2, Mh2))*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.0625*MW2*C0i(cc0, m12, m22, m32, Mh2, MHH2, MW2)*
   Csc(beta)*(pow(Cos(alp - beta), 2))*Sec(beta)*Sin(alp - beta)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) + 
 (0.0625*MW2*C0i(cc0, m12, m32, m22, Mh2, MHH2, MW2)*Csc(beta)*
   (pow(Cos(alp - beta), 2))*Sec(beta)*Sin(alp - beta)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) + 
 (0.0625*C0i(cc00, m12, m22, m32, Mh2, MHH2, MHp2)*Csc(beta)*
   (pow(Cos(alp - beta), 2))*Sec(beta)*Sin(alp - beta)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) - 
 (0.0625*C0i(cc00, m12, m22, m32, Mh2, MHH2, MW2)*Csc(beta)*
   (pow(Cos(alp - beta), 2))*Sec(beta)*Sin(alp - beta)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) + 
 (0.0625*C0i(cc00, m12, m32, m22, Mh2, MHH2, MHp2)*Csc(beta)*
   (pow(Cos(alp - beta), 2))*Sec(beta)*Sin(alp - beta)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) - 
 (0.0625*C0i(cc00, m12, m32, m22, Mh2, MHH2, MW2)*Csc(beta)*
   (pow(Cos(alp - beta), 2))*Sec(beta)*Sin(alp - beta)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) + 
 (0.0078125*Csc(2*beta)*Re(A0i(aa0, MHH2))*Sin(alp - beta)*
   (2*(-5*M2 + Mh2 + 5*MHH2)*Sin(2*alp) - 2*M2*Sin(2*(alp - 2*beta)) + 
    Mh2*Sin(4*alp - 2*beta) - MHH2*Sin(4*alp - 2*beta) - 8*M2*Sin(2*beta) - 
    4*MA02*Sin(2*beta) + Mh2*Sin(2*beta) + 11*MHH2*Sin(2*beta)))/
  (MA02*Pi2*v2) - (0.015625*(pow(Cos(alp - beta), 2))*
   (pow(Csc(2*beta), 2))*Re(B0i(bb0, Mh2, Mh2, MHH2))*Sin(alp - beta)*
   (-9*M2*Mh2 - 9*M2*MHH2 + 5*Mh2*MHH2 + 8*(pow(M2, 2)) + 
    2*(pow(Mh, 4)) + 2*(pow(MHH2, 2)) - 
    Cos(4*alp)*(5*Mh2*MHH2 - 9*M2*(Mh2 + MHH2) + 9*(pow(M2, 2)) + 
      2*(pow(Mh, 4)) + 2*(pow(MHH2, 2))) + 
    (pow(M2, 2))*(pow(Cos(beta), 4)) - 6*(pow(M2, 2))*
     (pow(Cos(beta), 2))*(pow(Sin(beta), 2)) + 
    (pow(M2, 2))*(pow(Sin(beta), 4)) - 2*M2*(Mh2 - MHH2)*Sin(2*alp)*
     Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.015625*(pow(Cos(alp - beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(bb0, MHH2, Mh2, MHH2))*Sin(alp - beta)*
   (-9*M2*Mh2 - 9*M2*MHH2 + 5*Mh2*MHH2 + 8*(pow(M2, 2)) + 
    2*(pow(Mh, 4)) + 2*(pow(MHH2, 2)) - 
    Cos(4*alp)*(5*Mh2*MHH2 - 9*M2*(Mh2 + MHH2) + 9*(pow(M2, 2)) + 
      2*(pow(Mh, 4)) + 2*(pow(MHH2, 2))) + 
    (pow(M2, 2))*(pow(Cos(beta), 4)) - 6*(pow(M2, 2))*
     (pow(Cos(beta), 2))*(pow(Sin(beta), 2)) + 
    (pow(M2, 2))*(pow(Sin(beta), 4)) - 2*M2*(Mh2 - MHH2)*Sin(2*alp)*
     Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0078125*A0i(aa0, MHH2)*Csc(2*beta)*Sin(2*(alp - beta))*
   (2*(M2 - MA02)*Sin(alp - 3*beta) + (-Mh2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 + 2*MA02 - Mh2 - 3*MHH2)*Sin(alp + beta)))/(MA02*Pi2*v2) + 
 (0.01171875*(pow(Csc(2*beta), 2))*Re(B0i(bb0, Mh2, MHH2, MHH2))*
   Sin(2*(alp - beta))*((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta))*
   (-(M2*Sin(alp - 3*beta)) + (M2 - MHH2)*Sin(3*alp - beta) + 
    (-2*M2 + 3*MHH2)*Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.01171875*(pow(Csc(2*beta), 2))*Re(B0i(bb0, MHH2, MHH2, MHH2))*
   Sin(2*(alp - beta))*((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta))*
   (-(M2*Sin(alp - 3*beta)) + (M2 - MHH2)*Sin(3*alp - beta) + 
    (-2*M2 + 3*MHH2)*Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0009765625*Cos(alp - beta)*((2*MA02 - Mh2)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh2, MA02, MA02))*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0009765625*Cos(alp - beta)*((-2*MA02 + Mh2)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MA02, MA02))*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0078125*(MA02 - MHH2)*B0i(bb0, 0, MA02, MHH2)*Csc(2*beta)*
   Sin(2*(alp - beta))*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/(MA02*Pi2*v2) + 
 (0.0078125*(MA02 - MHH2)*B0i(bb0, MA02, MA02, MHH2)*Csc(2*beta)*
   Sin(2*(alp - beta))*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/(MA02*Pi2*v2) - 
 (0.001953125*Cos(alp - beta)*((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh2 + 2*MHp2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh2, MHp2, MHp2))*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.001953125*Cos(alp - beta)*((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh2 + 2*MHp2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MHp2, MHp2))*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.5*Re(B0i(bb00, 0, MW2, MW2))*(5*Cos(w)*Cot(w) + 
    (2 + 3*Cos(2*w))*Csc(w)*Sin(alp - beta) - Sin(w)))/
  (Pi2*v2*(pow(Csc(w), 3))) - 
 (0.25*A0i(aa0, ME2)*Cos(alp - beta)*(pow(ME, 2))*Tan(beta))/
  (MA02*Pi2*v2) - (0.125*B0i(bb1, MA02, ME2, ME2)*Cos(alp - beta)*
   (pow(ME, 2))*Tan(beta))/(Pi2*v2) - 
 (0.25*A0i(aa0, ML2)*Cos(alp - beta)*(pow(ML, 2))*Tan(beta))/
  (MA02*Pi2*v2) - (0.125*B0i(bb1, MA02, ML2, ML2)*Cos(alp - beta)*
   (pow(ML, 2))*Tan(beta))/(Pi2*v2) - 
 (0.25*A0i(aa0, MM2)*Cos(alp - beta)*(pow(MM, 2))*Tan(beta))/
  (MA02*Pi2*v2) - (0.125*B0i(bb1, MA02, MM2, MM2)*Cos(alp - beta)*
   (pow(MM, 2))*Tan(beta))/(Pi2*v2) - 
 (0.25*Cos(alp)*(pow(ME, 4))*Re(B0i(bb0, Mh2, ME2, ME2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.25*Cos(alp)*(pow(ML, 4))*Re(B0i(bb0, Mh2, ML2, ML2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.25*Cos(alp)*(pow(MM, 4))*Re(B0i(bb0, Mh2, MM2, MM2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*Cos(alp)*(pow(ME, 4))*Re(B0i(bb0, MHH2, ME2, ME2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*Cos(alp)*(pow(ML, 4))*Re(B0i(bb0, MHH2, ML2, ML2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*Cos(alp)*(pow(MM, 4))*Re(B0i(bb0, MHH2, MM2, MM2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.125*MHH2*Cos(alp)*(pow(ME, 2))*Re(B0i(bb1, MHH2, ME2, ME2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.125*MHH2*Cos(alp)*(pow(ML, 2))*Re(B0i(bb1, MHH2, ML2, ML2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.125*MHH2*Cos(alp)*(pow(MM, 2))*Re(B0i(bb1, MHH2, MM2, MM2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.015625*Re(A0i(aa0, ME2))*(-21*MA02 + MA02*(-8 - 10*Cos(4*w) + Cos(6*w))*
     (pow(Csc(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
       2)) - 21*MA02*Cos(beta)*Sin(alp) - 8*Cos(beta)*(pow(ME, 2))*
     Sin(alp) + 8*(pow(ME, 2))*Sec(beta)*Sin(alp) + 
    21*MA02*Cos(alp)*Sin(beta) + 16*Cos(alp)*(pow(ME, 2))*Sin(beta) + 
    21*MA02*(pow(Cot(w), 2))*(1 + Cos(beta)*Sin(alp) - 
      Cos(alp)*Sin(beta)) + 8*(pow(ME, 2))*Sin(alp)*Sin(beta)*Tan(beta)))/
  (MA02*Pi2*v2) + (0.015625*Re(A0i(aa0, ML2))*
   (-21*MA02 + MA02*(-8 - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
     (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)) - 
    21*MA02*Cos(beta)*Sin(alp) - 8*Cos(beta)*(pow(ML, 2))*Sin(alp) + 
    8*(pow(ML, 2))*Sec(beta)*Sin(alp) + 21*MA02*Cos(alp)*Sin(beta) + 
    16*Cos(alp)*(pow(ML, 2))*Sin(beta) + 21*MA02*(pow(Cot(w), 2))*
     (1 + Cos(beta)*Sin(alp) - Cos(alp)*Sin(beta)) + 
    8*(pow(ML, 2))*Sin(alp)*Sin(beta)*Tan(beta)))/(MA02*Pi2*v2) + 
 (0.015625*Re(A0i(aa0, MM2))*(-21*MA02 + MA02*(-8 - 10*Cos(4*w) + Cos(6*w))*
     (pow(Csc(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
       2)) - 21*MA02*Cos(beta)*Sin(alp) - 8*Cos(beta)*(pow(MM, 2))*
     Sin(alp) + 8*(pow(MM, 2))*Sec(beta)*Sin(alp) + 
    21*MA02*Cos(alp)*Sin(beta) + 16*Cos(alp)*(pow(MM, 2))*Sin(beta) + 
    21*MA02*(pow(Cot(w), 2))*(1 + Cos(beta)*Sin(alp) - 
      Cos(alp)*Sin(beta)) + 8*(pow(MM, 2))*Sin(alp)*Sin(beta)*Tan(beta)))/
  (MA02*Pi2*v2)
;
}

ComplexType dKappahZZ(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32)
{
return 0. + (0.75*MT2*A0i(aa0, MT2)*Cos(alp - beta)*Cot(beta))/(MA02*Pi2*v2) + 
 (0.375*MT2*B0i(bb1, MA02, MT2, MT2)*Cos(alp - beta)*Cot(beta))/(Pi2*v2) + 
 (0.16666666666666666*MT2*C0i(cc00, m22, m32, m12, MT2, MT2, MT2)*
   (-9 + 4*Cos(2*w) - 4*Cos(4*w))*(-1 + Cos(alp)*Csc(beta)))/(Pi2*v2) + 
 (0.020833333333333332*MT2*C0i(cc1, m22, m32, m12, MT2, MT2, MT2)*
   (9*(m12 + 5*m22 - m32) - 16*m22*Cos(2*w) + 16*m22*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta)))/(Pi2*v2) + 
 (0.041666666666666664*MT2*C0i(cc2, m22, m32, m12, MT2, MT2, MT2)*
   (9*(2*m12 + m22 - m32) - 4*(m12 + m22 - m32)*Cos(2*w) + 
    4*(m12 + m22 - m32)*Cos(4*w))*(-1 + Cos(alp)*Csc(beta)))/(Pi2*v2) + 
 (0.020833333333333332*MT2*C0i(cc0, m22, m32, m12, MT2, MT2, MT2)*
   (9*(m12 + m22 - m32 + 4*MT2) - 4*(m12 + m22 - m32)*Cos(2*w) + 
    4*(m12 + m22 - m32)*Cos(4*w))*(-1 + Cos(alp)*Csc(beta)))/(Pi2*v2) - 
 (0.0234375*A0i(aa0, MA02)*Cos(alp - beta)*(3*(Mh2 - MHH2)*Cos(2*alp) + 
    (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 4*(-2*M2 + Mh2 + MHH2)*Cos(2*beta))*
   Csc(2*beta))/(MA02*Pi2*v2) - (0.015625*A0i(aa0, MHp2)*Cos(alp - beta)*
   (3*(Mh2 - MHH2)*Cos(2*alp) + (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 
    4*(-2*M2 + Mh2 + MHH2)*Cos(2*beta))*Csc(2*beta))/(MA02*Pi2*v2) + 
 (0.015625*B0i(bb0, m12, MA02, MA02)*((2*MA02 - Mh2)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh2)*Cos(alp + beta))*Csc(2*beta))/(Pi2*v2) + 
 (0.75*A0i(aa0, MB2)*Cos(alp - beta)*Cot(beta)*(pow(MB, 2)))/
  (MA02*Pi2*v2) + (0.375*B0i(bb1, MA02, MB2, MB2)*Cos(alp - beta)*Cot(beta)*
   (pow(MB, 2)))/(Pi2*v2) - 
 (0.16666666666666666*C0i(cc00, m22, m32, m12, MB2, MB2, MB2)*
   (6 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2)))/
  (Pi2*v2) + (0.041666666666666664*B0i(bb0, m32, MB2, MB2)*
   (15 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2)))/
  (Pi2*v2) + (0.020833333333333332*C0i(cc1, m22, m32, m12, MB2, MB2, MB2)*
   (9*m12 + 33*m22 - 9*m32 + 8*m22*Cos(2*w) + 4*m22*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(MB, 2)))/(Pi2*v2) + 
 (0.041666666666666664*C0i(cc2, m22, m32, m12, MB2, MB2, MB2)*
   (15*m12 + 6*m22 - 6*m32 + 2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2)))/
  (Pi2*v2) + (0.020833333333333332*C0i(cc0, m22, m32, m12, MB2, MB2, MB2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*(2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w) + 6*(m12 + m22 - m32 + 6*(pow(MB, 2)))))/
  (Pi2*v2) + (0.75*A0i(aa0, MC2)*Cos(alp - beta)*Cot(beta)*(pow(MC, 2)))/
  (MA02*Pi2*v2) + (0.375*B0i(bb1, MA02, MC2, MC2)*Cos(alp - beta)*Cot(beta)*
   (pow(MC, 2)))/(Pi2*v2) + 
 (0.16666666666666666*C0i(cc00, m22, m32, m12, MC2, MC2, MC2)*
   (-9 + 4*Cos(2*w) - 4*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2)))/
  (Pi2*v2) + (0.020833333333333332*C0i(cc1, m22, m32, m12, MC2, MC2, MC2)*
   (9*(m12 + 5*m22 - m32) - 16*m22*Cos(2*w) + 16*m22*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(MC, 2)))/(Pi2*v2) + 
 (0.041666666666666664*C0i(cc2, m22, m32, m12, MC2, MC2, MC2)*
   (9*(2*m12 + m22 - m32) - 4*(m12 + m22 - m32)*Cos(2*w) + 
    4*(m12 + m22 - m32)*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2)))/
  (Pi2*v2) + (0.020833333333333332*C0i(cc0, m22, m32, m12, MC2, MC2, MC2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*
   (-4*(m12 + m22 - m32)*Cos(2*w) + 4*(m12 + m22 - m32)*Cos(4*w) + 
    9*(m12 + m22 - m32 + 4*(pow(MC, 2)))))/(Pi2*v2) + 
 (0.75*A0i(aa0, MD2)*Cos(alp - beta)*Cot(beta)*(pow(MD, 2)))/
  (MA02*Pi2*v2) + (0.375*B0i(bb1, MA02, MD2, MD2)*Cos(alp - beta)*Cot(beta)*
   (pow(MD, 2)))/(Pi2*v2) - 
 (0.16666666666666666*C0i(cc00, m22, m32, m12, MD2, MD2, MD2)*
   (6 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(MD, 2)))/
  (Pi2*v2) + (0.041666666666666664*B0i(bb0, m32, MD2, MD2)*
   (15 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(MD, 2)))/
  (Pi2*v2) + (0.020833333333333332*C0i(cc1, m22, m32, m12, MD2, MD2, MD2)*
   (9*m12 + 33*m22 - 9*m32 + 8*m22*Cos(2*w) + 4*m22*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(MD, 2)))/(Pi2*v2) + 
 (0.041666666666666664*C0i(cc2, m22, m32, m12, MD2, MD2, MD2)*
   (15*m12 + 6*m22 - 6*m32 + 2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(MD, 2)))/
  (Pi2*v2) + (0.020833333333333332*C0i(cc0, m22, m32, m12, MD2, MD2, MD2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MD, 2))*(2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w) + 6*(m12 + m22 - m32 + 6*(pow(MD, 2)))))/
  (Pi2*v2) + (0.75*A0i(aa0, MS2)*Cos(alp - beta)*Cot(beta)*(pow(MS, 2)))/
  (MA02*Pi2*v2) + (0.375*B0i(bb1, MA02, MS2, MS2)*Cos(alp - beta)*Cot(beta)*
   (pow(MS, 2)))/(Pi2*v2) - 
 (0.16666666666666666*C0i(cc00, m22, m32, m12, MS2, MS2, MS2)*
   (6 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(MS, 2)))/
  (Pi2*v2) + (0.041666666666666664*B0i(bb0, m32, MS2, MS2)*
   (15 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(MS, 2)))/
  (Pi2*v2) + (0.020833333333333332*C0i(cc1, m22, m32, m12, MS2, MS2, MS2)*
   (9*m12 + 33*m22 - 9*m32 + 8*m22*Cos(2*w) + 4*m22*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(MS, 2)))/(Pi2*v2) + 
 (0.041666666666666664*C0i(cc2, m22, m32, m12, MS2, MS2, MS2)*
   (15*m12 + 6*m22 - 6*m32 + 2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(MS, 2)))/
  (Pi2*v2) + (0.020833333333333332*C0i(cc0, m22, m32, m12, MS2, MS2, MS2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MS, 2))*(2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w) + 6*(m12 + m22 - m32 + 6*(pow(MS, 2)))))/
  (Pi2*v2) + (0.75*A0i(aa0, MU2)*Cos(alp - beta)*Cot(beta)*(pow(MU, 2)))/
  (MA02*Pi2*v2) + (0.375*B0i(bb1, MA02, MU2, MU2)*Cos(alp - beta)*Cot(beta)*
   (pow(MU, 2)))/(Pi2*v2) + 
 (0.16666666666666666*C0i(cc00, m22, m32, m12, MU2, MU2, MU2)*
   (-9 + 4*Cos(2*w) - 4*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(MU, 2)))/
  (Pi2*v2) + (0.020833333333333332*C0i(cc1, m22, m32, m12, MU2, MU2, MU2)*
   (9*(m12 + 5*m22 - m32) - 16*m22*Cos(2*w) + 16*m22*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(MU, 2)))/(Pi2*v2) + 
 (0.041666666666666664*C0i(cc2, m22, m32, m12, MU2, MU2, MU2)*
   (9*(2*m12 + m22 - m32) - 4*(m12 + m22 - m32)*Cos(2*w) + 
    4*(m12 + m22 - m32)*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(MU, 2)))/
  (Pi2*v2) + (0.020833333333333332*C0i(cc0, m22, m32, m12, MU2, MU2, MU2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MU, 2))*
   (-4*(m12 + m22 - m32)*Cos(2*w) + 4*(m12 + m22 - m32)*Cos(4*w) + 
    9*(m12 + m22 - m32 + 4*(pow(MU, 2)))))/(Pi2*v2) - 
 (0.015625*(MA02 - Mh2)*B0i(bb0, 0, MA02, Mh2)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*Csc(2*beta)*(pow(Cos(alp - beta), 2)))/
  (MA02*Pi2*v2) - (0.015625*(MA02 - Mh2)*B0i(bb0, MA02, MA02, Mh2)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*Csc(2*beta)*(pow(Cos(alp - beta), 2)))/
  (MA02*Pi2*v2) - (0.0625*C0i(cc00, m12, m32, m22, MA02, MA02, Mh2)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*Csc(2*beta)*(pow(Cos(alp - beta), 2)))/(Pi2*v2) + 
 (0.015625*A0i(aa0, Mh2)*(2*(M2 - MA02)*Cos(alp - 3*beta) + 
    (-Mh2 + MHH2)*Cos(3*alp - beta) + (2*M2 + 2*MA02 - 3*Mh2 - MHH2)*
     Cos(alp + beta))*Csc(2*beta)*(pow(Cos(alp - beta), 2)))/
  (MA02*Pi2*v2) - (0.03125*B0i(bb0, m12, MHp2, MHp2)*
   ((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*Csc(2*beta)*(pow(Cos(2*w), 2)))/(Pi2*v2) + 
 (0.125*C0i(cc00, m22, m32, m12, MHp2, MHp2, MHp2)*
   ((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*Csc(2*beta)*(pow(Cos(2*w), 2)))/(Pi2*v2) - 
 1.*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)) - 
 (0.03125*MW2*B0i(bb0, m32, MW2, MW2)*(11 + 20*Cos(2*w) + Cos(4*w))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.03125*MW2*C0i(cc0, m22, m32, m12, MW2, MW2, MW2)*
   (-10*m12 + 19*m22 + 11*m32 - 3*Mh2 + 8*MW2 + 
    4*(-4*m12 + 6*m22 + 4*m32 + Mh2 + 6*MW2)*Cos(2*w) - 
    (6*m12 - 5*m22 - 5*m32 + Mh2)*Cos(4*w))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.125*C0i(cc00, m22, m32, m12, MW2, MW2, MW2)*
   (Mh2 + 14*MW2 + 16*MW2*Cos(2*w) + (Mh2 + 6*MW2)*Cos(4*w))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) + 
 (0.03125*B0i(bb0, m12, MW2, MW2)*(Mh2 + 18*MW2 + 24*MW2*Cos(2*w) + 
    (Mh2 + 6*MW2)*Cos(4*w))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.25*MW2*C0i(cc2, m22, m32, m12, MW2, MW2, MW2)*
   (3*m12 + 2*m22 - 2*m32 + m12*Cos(2*w))*(pow(Cos(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.125*MW2*C0i(cc1, m22, m32, m12, MW2, MW2, MW2)*
   (m12 + 9*m22 - m32 + (m12 + m22 - m32)*Cos(2*w))*(pow(Cos(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) + 
 (0.25*C0i(cc0, m22, m12, m32, Mh2, MZ2, MZ2)*(pow(MW, 4))*
   (pow(Sec(w), 4))*(1 + (pow(Sin(alp - beta), 3))))/(Pi2*v2) + 
 (0.00390625*((Mh2 - MHH2)*Cos(3*alp - 5*beta) + 
    (-8*M2 - 4*MA02 + 5*Mh2 + 3*MHH2)*Cos(alp - 3*beta) + 
    3*Mh2*Cos(3*alp - beta) - 3*MHH2*Cos(3*alp - beta) - 
    8*M2*Cos(alp + beta) + 4*MA02*Cos(alp + beta) + 7*Mh2*Cos(alp + beta) + 
    MHH2*Cos(alp + beta))*Csc(2*beta)*Re(A0i(aa0, MA02)))/(MA02*Pi2*v2) + 
 (0.0078125*((Mh2 - MHH2)*Cos(3*alp - 5*beta) + 
    (-8*M2 - 2*MA02 + 5*Mh2 + 3*MHH2)*Cos(alp - 3*beta) + 
    3*Mh2*Cos(3*alp - beta) - 3*MHH2*Cos(3*alp - beta) - 
    8*M2*Cos(alp + beta) + 2*MA02*Cos(alp + beta) + 7*Mh2*Cos(alp + beta) + 
    MHH2*Cos(alp + beta) - MA02*Cos(alp - 3*beta - 4*w) + 
    MA02*Cos(alp + beta - 4*w) - 8*MA02*Cos(alp - 3*beta - 2*w) + 
    8*MA02*Cos(alp + beta - 2*w) - 8*MA02*Cos(alp - 3*beta + 2*w) + 
    8*MA02*Cos(alp + beta + 2*w) - MA02*Cos(alp - 3*beta + 4*w) + 
    MA02*Cos(alp + beta + 4*w))*Csc(2*beta)*Re(A0i(aa0, MHp2)))/
  (MA02*Pi2*v2) + (0.375*(MT2 - MW2)*Cos(2*w)*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MW2, MB2, MT2)))/(Pi2*v2) + 
 (0.375*Cos(2*w)*(-MW2 + (pow(MU, 2)))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MW2, MD2, MU2)))/(Pi2*v2) - 
 (0.125*MW2*Cos(2*w)*(pow(Csc(w), 2))*(1 + (pow(Sin(alp - beta), 3)))*
   Re(B0i(bb0, MW2, Mh2, MW2)))/(Pi2*v2) + 
 (0.015625*MW2*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MW2, MW2, MZ2)))/(Pi2*v2) - 
 (0.09375*(-1 + 3*Cos(2*w))*(pow(MB, 2))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MZ2, MB2, MB2)))/(Pi2*v2) - 
 (0.09375*(-1 + 3*Cos(2*w))*(pow(MC, 2))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MZ2, MC2, MC2)))/(Pi2*v2) - 
 (0.09375*(-1 + 3*Cos(2*w))*(pow(MD, 2))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MZ2, MD2, MD2)))/(Pi2*v2) - 
 (0.03125*(-1 + 3*Cos(2*w))*(pow(ME, 2))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.0625*MW2*(-1 + 3*Cos(2*w))*(pow(Csc(w), 2))*(pow(Sec(w), 2))*
   (1 + (pow(Sin(alp - beta), 3)))*Re(B0i(bb0, MZ2, Mh2, MZ2)))/
  (Pi2*v2) - (0.03125*(-1 + 3*Cos(2*w))*(pow(ML, 2))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MZ2, ML2, ML2)))/(Pi2*v2) - 
 (0.03125*(-1 + 3*Cos(2*w))*(pow(MM, 2))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MZ2, MM2, MM2)))/(Pi2*v2) - 
 (0.09375*(-1 + 3*Cos(2*w))*(pow(MS, 2))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MZ2, MS2, MS2)))/(Pi2*v2) - 
 (0.09375*MT2*(-1 + 3*Cos(2*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MZ2, MT2, MT2)))/(Pi2*v2) - 
 (0.09375*(-1 + 3*Cos(2*w))*(pow(MU, 2))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MZ2, MU2, MU2)))/(Pi2*v2) - 
 (0.015625*MW2*(27 + 12*Cos(2*w) + 17*Cos(4*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb0, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.5*(2 + 3*Cos(2*w))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.125*Cos(2*w)*(pow(Csc(w), 2))*(1 + (pow(Sin(alp - beta), 3)))*
   Re(B0i(bb00, MW2, Mh2, MW2)))/(Pi2*v2) + 
 (0.125*Cos(2*w)*(pow(Csc(w), 2))*(pow(Sin(alp - beta), 3))*
   Re(B0i(bb00, MW2, MHH2, MHp2)))/(Pi2*v2) + 
 (0.1875*(-1 + 3*Cos(2*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, 0, 0)))/(Pi2*v2) - 
 (0.0625*(-1 + 3*Cos(2*w))*(pow(Csc(w), 2))*(pow(Sin(alp - beta), 3))*
   Re(B0i(bb00, MZ2, MA02, MHH2)))/(Pi2*v2) + 
 (0.010416666666666666*(-6 + 35*Cos(2*w) + 4*Cos(4*w) + 3*Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(bb00, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.020833333333333332*(-15 + 37*Cos(2*w) - 10*Cos(4*w) + 6*Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(bb00, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.010416666666666666*(-6 + 35*Cos(2*w) + 4*Cos(4*w) + 3*Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(bb00, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.03125*(-10 + 19*Cos(2*w) - 8*Cos(4*w) + 3*Cos(6*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.0625*(4 - 2*(pow(Cot(w), 2)) - (-1 + 3*Cos(2*w))*(pow(Csc(w), 2))*
     (pow(Sin(alp - beta), 3)))*Re(B0i(bb00, MZ2, Mh2, MZ2)))/(Pi2*v2) + 
 (0.03125*(-10 + 19*Cos(2*w) - 8*Cos(4*w) + 3*Cos(6*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.03125*(-10 + 19*Cos(2*w) - 8*Cos(4*w) + 3*Cos(6*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.010416666666666666*(-6 + 35*Cos(2*w) + 4*Cos(4*w) + 3*Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(bb00, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.020833333333333332*(-15 + 37*Cos(2*w) - 10*Cos(4*w) + 6*Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(bb00, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.020833333333333332*(-15 + 37*Cos(2*w) - 10*Cos(4*w) + 6*Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(bb00, MZ2, MU2, MU2)))/(Pi2*v2) - 
 (0.015625*(10 + 35*Cos(2*w) + 18*Cos(4*w) + 9*Cos(6*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb00, MZ2, MW2, MW2)))/(Pi2*v2) - 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (pow(Csc(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(bb1, MZ2, MC2, MC2)))/(Pi2*v2) - 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(Pi2*v2) - 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb1, MZ2, ME2, ME2)))/(Pi2*v2) - 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb1, MZ2, ML2, ML2)))/(Pi2*v2) - 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb1, MZ2, MM2, MM2)))/(Pi2*v2) - 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (pow(Csc(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(bb1, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (pow(Csc(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(bb1, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.015625*(pow(Csc(2*beta), 2))*(pow(Sin(alp - beta), 3))*
   (pow((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta), 2))*
   Re(B0i(dbb0, Mh2, MHH2, MHH2)))/(Pi2*v2) + 
 (0.03125*(-2*Mh2*MW2 + (pow(Mh, 4)) + 12*(pow(MW, 4)))*
   (1 + (pow(Sin(alp - beta), 3)))*Re(B0i(dbb0, Mh2, MW2, MW2)))/
  (Pi2*v2) + (0.25*(pow(MW, 4))*(pow(Sec(w), 4))*
   (1 + (pow(Sin(alp - beta), 3)))*Re(B0i(dbb0, MZ2, Mh2, MZ2)))/
  (Pi2*v2) - (0.25*MW2*(pow(Sec(w), 2))*(pow(Sin(alp - beta), 3))*
   Re(B0i(dbb00, MZ2, MA02, MHH2)))/(Pi2*v2) + 
 (0.08333333333333333*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb00, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.08333333333333333*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb00, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.08333333333333333*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb00, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.25*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb00, MZ2, ME2, ME2)))/(Pi2*v2) - 
 (0.25*MW2*(pow(Sec(w), 2))*(1 + (pow(Sin(alp - beta), 3)))*
   Re(B0i(dbb00, MZ2, Mh2, MZ2)))/(Pi2*v2) + 
 (0.25*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb00, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.25*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb00, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.08333333333333333*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb00, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.08333333333333333*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb00, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.08333333333333333*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb00, MZ2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*MW2*(7 + 8*Cos(2*w) + 3*Cos(4*w))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb00, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.125*Mh2*MW2*(1 + (pow(Sin(alp - beta), 3)))*
   Re(B0i(dbb1, Mh2, MW2, MW2)))/(Pi2*v2) + 
 (0.0625*Mh2*MW2*(pow(Sec(w), 2))*(1 + (pow(Sin(alp - beta), 3)))*
   Re(B0i(dbb1, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(MW, 4))*
   (pow(Sec(w), 4))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(dbb1, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.041666666666666664*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*(pow(MW, 4))*
   (pow(Sec(w), 4))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(dbb1, MZ2, MC2, MC2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(MW, 4))*
   (pow(Sec(w), 4))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(dbb1, MZ2, MD2, MD2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(MW, 4))*(pow(Sec(w), 4))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb1, MZ2, ME2, ME2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(MW, 4))*(pow(Sec(w), 4))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb1, MZ2, ML2, ML2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(MW, 4))*(pow(Sec(w), 4))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   Re(B0i(dbb1, MZ2, MM2, MM2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(MW, 4))*
   (pow(Sec(w), 4))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(dbb1, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.041666666666666664*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*(pow(MW, 4))*
   (pow(Sec(w), 4))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(dbb1, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.041666666666666664*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*(pow(MW, 4))*
   (pow(Sec(w), 4))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(dbb1, MZ2, MU2, MU2)))/(Pi2*v2) - 
 (0.09375*C0i(cc00, m22, m12, m32, MA02, Mh2, Mh2)*
   (M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh2)*Cos(alp + beta))*Csc(beta)*(pow(Cos(alp - beta), 2))*
   Sec(beta))/(Pi2*v2) - (0.03125*C0i(cc00, m12, m32, m22, MA02, MA02, MHH2)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*Csc(beta)*(pow(Sin(alp - beta), 2))*Sec(beta))/
  (Pi2*v2) + (0.001953125*(pow(Cos(w), 2))*Re(A0i(aa0, Mh2))*
   (-64*(pow(Csc(2*w), 2)) - (2*M2*Cos(3*alp - 5*beta) + 
      (14*M2 + 4*MA02 - 11*Mh2 - MHH2)*Cos(alp - 3*beta) - 
      Mh2*Cos(5*alp - 3*beta) + MHH2*Cos(5*alp - 3*beta) + 
      10*M2*Cos(3*alp - beta) - 13*Mh2*Cos(3*alp - beta) + 
      MHH2*Cos(3*alp - beta) + 22*M2*Cos(alp + beta) - 
      4*MA02*Cos(alp + beta) - 23*Mh2*Cos(alp + beta) - MHH2*Cos(alp + beta))*
     Csc(beta)*(pow(MA02, -1))*(pow(Csc(w), 2))*(pow(Sec(w), 2))*
     Sec(beta)))/(Pi2*v2*(pow(Csc(w), 2))) + 
 (0.09375*C0i(cc00, m12, m32, m22, Mh2, Mh2, MZ2)*(pow(Cos(w), 2))*
   (-16*Mh2*(pow(Csc(2*w), 2)) - 
    (M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
      (2*M2 - 3*Mh2)*Cos(alp + beta))*Csc(beta)*(pow(Csc(w), 2))*
     (pow(Sec(w), 2))*(pow(Sin(alp - beta), 2))*Sec(beta)))/
  (Pi2*v2*(pow(Csc(w), 2))) + 
 (0.09375*MW2*C0i(cc0, m12, m32, m22, Mh2, Mh2, MZ2)*(pow(Cos(w), 2))*
   (64*Mh2*(pow(Csc(2*w), 4)) + 
    (M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
      (2*M2 - 3*Mh2)*Cos(alp + beta))*Csc(beta)*(pow(Csc(w), 4))*
     (pow(Sec(w), 4))*(pow(Sin(alp - beta), 2))*Sec(beta)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.5*C0i(cc00, m22, m32, m12, ME2, ME2, ME2)*(2 - 2*Cos(2*w) + Cos(4*w))*
   (pow(ME, 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.125*B0i(bb0, m32, ME2, ME2)*(3 - 2*Cos(2*w) + Cos(4*w))*(pow(ME, 2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*C0i(cc1, m22, m32, m12, ME2, ME2, ME2)*
   (m12 + 9*m22 - m32 - 8*m22*Cos(2*w) + 4*m22*Cos(4*w))*(pow(ME, 2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.125*C0i(cc2, m22, m32, m12, ME2, ME2, ME2)*(3*m12 + 2*m22 - 2*m32 - 
    2*(m12 + m22 - m32)*Cos(2*w) + (m12 + m22 - m32)*Cos(4*w))*
   (pow(ME, 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*C0i(cc0, m22, m32, m12, ME2, ME2, ME2)*(pow(ME, 2))*
   (-2*(m12 + m22 - m32)*Cos(2*w) + (m12 + m22 - m32)*Cos(4*w) + 
    2*(m12 + m22 - m32 + 2*(pow(ME, 2))))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.5*C0i(cc00, m22, m32, m12, ML2, ML2, ML2)*(2 - 2*Cos(2*w) + Cos(4*w))*
   (pow(ML, 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.125*B0i(bb0, m32, ML2, ML2)*(3 - 2*Cos(2*w) + Cos(4*w))*(pow(ML, 2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*C0i(cc1, m22, m32, m12, ML2, ML2, ML2)*
   (m12 + 9*m22 - m32 - 8*m22*Cos(2*w) + 4*m22*Cos(4*w))*(pow(ML, 2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.125*C0i(cc2, m22, m32, m12, ML2, ML2, ML2)*(3*m12 + 2*m22 - 2*m32 - 
    2*(m12 + m22 - m32)*Cos(2*w) + (m12 + m22 - m32)*Cos(4*w))*
   (pow(ML, 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*C0i(cc0, m22, m32, m12, ML2, ML2, ML2)*(pow(ML, 2))*
   (-2*(m12 + m22 - m32)*Cos(2*w) + (m12 + m22 - m32)*Cos(4*w) + 
    2*(m12 + m22 - m32 + 2*(pow(ML, 2))))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.5*C0i(cc00, m22, m32, m12, MM2, MM2, MM2)*(2 - 2*Cos(2*w) + Cos(4*w))*
   (pow(MM, 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.125*B0i(bb0, m32, MM2, MM2)*(3 - 2*Cos(2*w) + Cos(4*w))*(pow(MM, 2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*C0i(cc1, m22, m32, m12, MM2, MM2, MM2)*
   (m12 + 9*m22 - m32 - 8*m22*Cos(2*w) + 4*m22*Cos(4*w))*(pow(MM, 2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.125*C0i(cc2, m22, m32, m12, MM2, MM2, MM2)*(3*m12 + 2*m22 - 2*m32 - 
    2*(m12 + m22 - m32)*Cos(2*w) + (m12 + m22 - m32)*Cos(4*w))*
   (pow(MM, 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*C0i(cc0, m22, m32, m12, MM2, MM2, MM2)*(pow(MM, 2))*
   (-2*(m12 + m22 - m32)*Cos(2*w) + (m12 + m22 - m32)*Cos(4*w) + 
    2*(m12 + m22 - m32 + 2*(pow(MM, 2))))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.75*Cos(alp)*Csc(beta)*(pow(MB, 4))*Re(B0i(bb0, Mh2, MB2, MB2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*Csc(beta)*(pow(MC, 4))*Re(B0i(bb0, Mh2, MC2, MC2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*Csc(beta)*(pow(MD, 4))*Re(B0i(bb0, Mh2, MD2, MD2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*Csc(beta)*(pow(MS, 4))*Re(B0i(bb0, Mh2, MS2, MS2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*Csc(beta)*(pow(MT, 4))*Re(B0i(bb0, Mh2, MT2, MT2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*Csc(beta)*(pow(MU, 4))*Re(B0i(bb0, Mh2, MU2, MU2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*Csc(beta)*(pow(MB, 4))*Re(B0i(bb0, MHH2, MB2, MB2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*Csc(beta)*(pow(MC, 4))*Re(B0i(bb0, MHH2, MC2, MC2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*Csc(beta)*(pow(MD, 4))*Re(B0i(bb0, MHH2, MD2, MD2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*Csc(beta)*(pow(MS, 4))*Re(B0i(bb0, MHH2, MS2, MS2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*Csc(beta)*(pow(MT, 4))*Re(B0i(bb0, MHH2, MT2, MT2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*Csc(beta)*(pow(MU, 4))*Re(B0i(bb0, MHH2, MU2, MU2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*Csc(beta)*(pow(MB, 2))*Re(B0i(bb1, MHH2, MB2, MB2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*Csc(beta)*(pow(MC, 2))*Re(B0i(bb1, MHH2, MC2, MC2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*Csc(beta)*(pow(MD, 2))*Re(B0i(bb1, MHH2, MD2, MD2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*Csc(beta)*(pow(MS, 2))*Re(B0i(bb1, MHH2, MS2, MS2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*MT2*Cos(alp)*Csc(beta)*Re(B0i(bb1, MHH2, MT2, MT2))*Sin(alp)*
   (Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*Csc(beta)*(pow(MU, 2))*Re(B0i(bb1, MHH2, MU2, MU2))*
   Sin(alp)*(Cos(alp)*Cot(beta) + Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*(pow(ME, 2))*Re(B0i(bb1, Mh2, ME2, ME2))*
   (2*(Mh2 - MHH2) + MHH2*Cos(2*alp - beta)*(pow(Sec(beta), 2))*
     Sin(alp) + (2*Mh2 - MHH2)*Sec(beta)*Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*(pow(ML, 2))*Re(B0i(bb1, Mh2, ML2, ML2))*
   (2*(Mh2 - MHH2) + MHH2*Cos(2*alp - beta)*(pow(Sec(beta), 2))*
     Sin(alp) + (2*Mh2 - MHH2)*Sec(beta)*Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*(pow(MM, 2))*Re(B0i(bb1, Mh2, MM2, MM2))*
   (2*(Mh2 - MHH2) + MHH2*Cos(2*alp - beta)*(pow(Sec(beta), 2))*
     Sin(alp) + (2*Mh2 - MHH2)*Sec(beta)*Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0625*(Mh2 - MHH2)*A0i(aa0, MW2)*(pow(Cos(alp - beta), 2))*
   Sin(alp - beta))/(MA02*Pi2*v2) + 
 (0.09375*(Mh2 - MHH2)*A0i(aa0, MZ2)*(pow(Cos(alp - beta), 2))*
   Sin(alp - beta))/(MA02*Pi2*v2) + 
 (0.0625*MW2*B0i(bb1, MA02, Mh2, MZ2)*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Sin(alp - beta))/(Pi2*v2) - 
 (0.0625*MW2*B0i(bb1, MA02, MHH2, MZ2)*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.0625*C0i(cc00, m22, m32, m12, MA02, Mh2, MZ2)*
   (MA02 - Mh2 - 2*MW2 + (MA02 - Mh2)*Cos(2*w))*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.0625*C0i(cc00, m32, m22, m12, MA02, Mh2, MZ2)*
   (MA02 - Mh2 - 2*MW2 + (MA02 - Mh2)*Cos(2*w))*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Sin(alp - beta))/(Pi2*v2) - 
 (0.0625*C0i(cc00, m22, m12, m32, MHH2, MZ2, MZ2)*
   (Mh2 + 4*MW2 + Mh2*Cos(2*w))*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.015625*Mh2*B0i(bb0, 0, Mh2, MZ2)*(-MA02 + Mh2 - 2*MW2 + 
    (-MA02 + Mh2)*Cos(2*w))*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Sin(alp - beta))/(MA02*Pi2*v2) + 
 (0.0625*C0i(cc00, m22, m32, m12, MA02, MHH2, MZ2)*
   (-MA02 + Mh2 + 2*MW2 + (-MA02 + Mh2)*Cos(2*w))*
   (pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*Sin(alp - beta))/
  (Pi2*v2) + (0.0625*C0i(cc00, m32, m22, m12, MA02, MHH2, MZ2)*
   (-MA02 + Mh2 + 2*MW2 + (-MA02 + Mh2)*Cos(2*w))*
   (pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*Sin(alp - beta))/
  (Pi2*v2) + (0.015625*MHH2*B0i(bb0, 0, MHH2, MZ2)*
   (MA02 - MHH2 + 2*MW2 + (MA02 - MHH2)*Cos(2*w))*
   (pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*Sin(alp - beta))/
  (MA02*Pi2*v2) + (0.015625*B0i(bb0, MA02, MHH2, MZ2)*
   (-(MHH2*(MHH2 - 2*MW2)) + MA02*(MHH2 + 2*MW2) + 
    (MA02 - MHH2)*MHH2*Cos(2*w))*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Sin(alp - beta))/(MA02*Pi2*v2) - 
 (0.015625*B0i(bb0, MA02, Mh2, MZ2)*(2*Mh2*MW2 + MA02*(Mh2 + 2*MW2) + 
    (MA02 - Mh2)*Mh2*Cos(2*w) - (pow(Mh, 4)))*
   (pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*Sin(alp - beta))/
  (MA02*Pi2*v2) + (0.25*C0i(cc0, m22, m12, m32, MHH2, MZ2, MZ2)*
   (pow(MW, 4))*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 4))*
   Sin(alp - beta))/(Pi2*v2) - 
 (0.015625*(MHH2*(Mh2 - 2*MW2) - MA02*(Mh2 + MHH2 + 2*MW2) + 
    (MA02 - Mh2)*(MA02 - MHH2)*Cos(2*w) + (pow(MA02, 2)))*
   (pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*
   Re(B0i(bb0, Mh2, MA02, MZ2))*Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*(Mh2*(MHH2 - MHp2) + MHp2*(MHp2 - MW2) - MHH2*(MHp2 + MW2))*
   (pow(Cos(alp - beta), 2))*Re(B0i(bb0, Mh2, MHp2, MW2))*
   Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0625*(Mh2*(MHH2 - MHp2) + MHp2*(MHp2 - MW2) - MHH2*(MHp2 + MW2))*
   (pow(Cos(alp - beta), 2))*Re(B0i(bb0, MHH2, MHp2, MW2))*
   Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.03125*(Mh2*MHH2 - 2*MHH2*MW2 + 12*(pow(MW, 4)))*
   (pow(Cos(alp - beta), 2))*Re(B0i(bb0, MHH2, MW2, MW2))*
   Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.001953125*(3*Mh2*MHH2 - 8*MHH2*MW2 + 4*MHH2*(Mh2 - 2*MW2)*Cos(2*w) + 
    Mh2*MHH2*Cos(4*w) + 96*(pow(MW, 4)))*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 4))*Re(B0i(bb0, MHH2, MZ2, MZ2))*Sin(alp - beta))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.125*MW2*Cos(2*w)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb0, MW2, MHH2, MW2))*Sin(alp - beta))/
  (Pi2*v2) + (0.0625*MW2*(-1 + 3*Cos(2*w))*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*(pow(Sec(w), 2))*Re(B0i(bb0, MZ2, MHH2, MZ2))*
   Sin(alp - beta))/(Pi2*v2) + (0.5*Cos(2*w)*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MHp2, MHp2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.125*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, MA02, MHp2))*
   Sin(alp - beta))/(Pi2*v2) + (0.125*Cos(2*w)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MW2, Mh2, MHp2))*Sin(alp - beta))/
  (Pi2*v2) + (0.125*Cos(2*w)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MW2, MHH2, MW2))*Sin(alp - beta))/
  (Pi2*v2) - (0.0625*(-1 + 3*Cos(2*w))*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MZ2, MA02, Mh2))*Sin(alp - beta))/
  (Pi2*v2) - (0.0625*(-1 + 3*Cos(2*w))*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MZ2, MHH2, MZ2))*Sin(alp - beta))/
  (Pi2*v2) - (0.0625*(-1 + 3*Cos(2*w))*(pow(Cos(2*w), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MZ2, MHp2, MHp2))*Sin(alp - beta))/
  (Pi2*v2) - (0.0625*MHH2*MW2*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Re(B0i(bb1, Mh2, MA02, MZ2))*Sin(alp - beta))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.125*MHH2*MW2*(pow(Cos(alp - beta), 2))*
   Re(B0i(bb1, Mh2, MHp2, MW2))*Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0625*MHH2*MW2*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*
   Re(B0i(bb1, MHH2, MA02, MZ2))*Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.125*MHH2*MW2*(pow(Cos(alp - beta), 2))*Re(B0i(bb1, MHH2, MHp2, MW2))*
   Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.125*MHH2*MW2*(pow(Cos(alp - beta), 2))*Re(B0i(bb1, MHH2, MW2, MW2))*
   Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*MHH2*MW2*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*
   Re(B0i(bb1, MHH2, MZ2, MZ2))*Sin(alp - beta))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.00390625*(pow((2*MA02 - Mh2)*Cos(alp - 3*beta) + 
      (4*M2 - 2*MA02 - 3*Mh2)*Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh2, MA02, MA02))*Sin(alp - beta))/(Pi2*v2) + 
 (0.015625*(-2*Mh2*MW2 - 2*MA02*(Mh2 + MW2) + (pow(MA02, 2)) + 
    (pow(Mh, 4)) + Cos(2*w)*(pow(MA02 - Mh2, 2)))*
   (pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, Mh2, MA02, MZ2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.03125*(pow(Cos(alp - beta), 2))*(pow(Csc(2*beta), 2))*
   (pow((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta), 2))*
   Re(B0i(dbb0, Mh2, Mh2, MHH2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.0078125*(pow((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + 
      (-4*M2 + 3*Mh2 + 2*MHp2)*Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(dbb0, Mh2, MHp2, MHp2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.0625*(MHp2*(MHp2 - MW2) - Mh2*(2*MHp2 + MW2) + (pow(Mh, 4)))*
   (pow(Cos(alp - beta), 2))*Re(B0i(dbb0, Mh2, MHp2, MW2))*
   Sin(alp - beta))/(Pi2*v2) + 
 (0.25*(pow(MW, 4))*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 4))*
   Re(B0i(dbb0, MZ2, MHH2, MZ2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MHp2, MHp2))*Sin(alp - beta))/
  (Pi2*v2) - (0.25*MW2*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb00, MZ2, MA02, Mh2))*Sin(alp - beta))/(Pi2*v2) - 
 (0.25*MW2*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb00, MZ2, MHH2, MZ2))*Sin(alp - beta))/(Pi2*v2) - 
 (0.25*MW2*(pow(Cos(2*w), 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb00, MZ2, MHp2, MHp2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.0625*Mh2*MW2*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb1, Mh2, MA02, MZ2))*Sin(alp - beta))/(Pi2*v2) + 
 (0.125*Mh2*MW2*(pow(Cos(alp - beta), 2))*Re(B0i(dbb1, Mh2, MHp2, MW2))*
   Sin(alp - beta))/(Pi2*v2) + (0.03125*Mh2*B0i(bb0, m12, MZ2, MZ2)*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.125*MW2*B0i(bb0, m22, Mh2, MZ2)*(pow(Sec(w), 2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.125*MW2*B0i(bb0, m32, Mh2, MZ2)*(pow(Sec(w), 2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.25*MW2*B0i(bb0, m22, MW2, MW2)*(pow(Sin(w), 4))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (1.125*MW2*(pow(Sin(w), 2))*Re(B0i(bb0, 0, MW2, MW2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.5*MW2*Cos(2*w)*Re(B0i(bb0, MW2, 0, MW2))*(1 + Sin(alp - beta)))/
  (Pi2*v2) + (0.375*Cos(2*w)*(pow(MC, 2))*(pow(Csc(w), 2))*
   Re(B0i(bb0, MW2, MC2, MS2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MB2, MB2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MC2, MC2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MD2, MD2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.5*(-1 + 2*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, ME2, ME2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.5*(-1 + 2*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, ML2, ML2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.5*(-1 + 2*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, MM2, MM2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MS2, MS2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MT2, MT2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MU2, MU2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.25*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, 0, ME2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.25*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, 0, ML2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.25*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, 0, MM2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (1.*Cos(2*w)*Re(B0i(bb00, MW2, 0, MW2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.75*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, MB2, MT2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.75*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, MC2, MS2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.75*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, MD2, MU2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.125*Cos(2*w)*(5 + 4*Cos(2*w))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MW2, MW2, MZ2))*(1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MB2, MB2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MC2, MC2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MD2, MD2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ME2, ME2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ML2, ML2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MM2, MM2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MS2, MS2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MT2, MT2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MU2, MU2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.25*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MW2, MW2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.125*MW2*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb1, MW2, 0, ME2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.125*MW2*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb1, MW2, 0, ML2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.125*MW2*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb1, MW2, 0, MM2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.375*MW2*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb1, MW2, MB2, MT2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.375*MW2*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb1, MW2, MC2, MS2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.375*MW2*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb1, MW2, MD2, MU2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.25*MW2*Cos(2*w)*(pow(Cot(w), 2))*Re(B0i(bb1, MW2, MW2, MZ2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.1875*MW2*(pow(Csc(w), 2))*Re(B0i(bb1, MZ2, 0, 0))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.25*MW2*(pow(Cos(w), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb1, MZ2, MW2, MW2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.375*MW2*(pow(MB, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, MB2, MB2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.375*MW2*(pow(MC, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, MC2, MC2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.375*MW2*(pow(MD, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, MD2, MD2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.125*MW2*(pow(ME, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, ME2, ME2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.125*MW2*(pow(ML, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, ML2, ML2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.125*MW2*(pow(MM, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, MM2, MM2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.375*MW2*(pow(MS, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, MS2, MS2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.375*MT2*MW2*(pow(Sec(w), 2))*Re(B0i(dbb0, MZ2, MT2, MT2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.375*MW2*(pow(MU, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, MU2, MU2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.125*(5 + 9*Cos(2*w))*(pow(MW, 4))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, MW2, MW2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MB2, MB2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MC2, MC2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MD2, MD2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ME2, ME2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ML2, ML2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MM2, MM2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MS2, MS2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MT2, MT2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MU2, MU2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (1.5*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MW2, MW2))*
   (1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.75*MW2*(pow(Sec(w), 2))*Re(B0i(dbb00, MZ2, 0, 0))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.375*(pow(MW, 4))*(pow(Sec(w), 4))*Re(B0i(dbb1, MZ2, 0, 0))*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.5*(pow(MW, 4))*Re(B0i(dbb1, MZ2, MW2, MW2))*(1 + Sin(alp - beta)))/
  (Pi2*v2) + (0.03125*C0i(cc00, m22, m12, m32, Mh2, MZ2, MZ2)*
   (Mh2 + 4*MW2 + Mh2*Cos(2*w))*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   (-3 + Cos(2*(alp - beta)) + 2*Sin(alp - beta)))/(Pi2*v2) - 
 (0.0009765625*(-8*Mh2*MW2 + 4*Mh2*(Mh2 - 2*MW2)*Cos(2*w) + 
    3*(pow(Mh, 4)) + Cos(4*w)*(pow(Mh, 4)) + 96*(pow(MW, 4)))*
   (pow(Sec(w), 4))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2))*Re(B0i(dbb0, Mh2, MZ2, MZ2))*(-3 + Cos(2*(alp - beta)) + 
    2*Sin(alp - beta)))/(Pi2*v2) + 
 (0.0078125*(pow(Csc(w), 2))*Re(A0i(aa0, MW2))*
   (10 + Cos(2*w) + 34*Cos(4*w) + 3*Cos(6*w) + 
    (10*MA02 - 2*Mh2 + 2*MHH2 - 2*(Mh2 - MHH2)*Cos(2*(alp - beta)) + 
      (Mh2 - MHH2)*Cos(2*(alp - beta - w)) + MA02*Cos(2*w) + 2*Mh2*Cos(2*w) - 
      2*MHH2*Cos(2*w) + 34*MA02*Cos(4*w) + 3*MA02*Cos(6*w) + 
      Mh2*Cos(2*(alp - beta + w)) - MHH2*Cos(2*(alp - beta + w)))*
     (pow(MA02, -1))*Sin(alp - beta)))/(Pi2*v2) + 
 (0.00390625*(pow(Csc(w), 2))*Re(A0i(aa0, MZ2))*
   (-4*(7 + 11*Cos(2*w) + 6*Cos(4*w)) - 
    (28*MA02 + 2*Mh2 - 2*MHH2 + 2*(Mh2 - MHH2)*Cos(2*(alp - beta)) + 
      (-Mh2 + MHH2)*Cos(2*(alp - beta - w)) + 44*MA02*Cos(2*w) - 
      2*Mh2*Cos(2*w) + 2*MHH2*Cos(2*w) + 24*MA02*Cos(4*w) - 
      Mh2*Cos(2*(alp - beta + w)) + MHH2*Cos(2*(alp - beta + w)))*
     (pow(MA02, -1))*Sin(alp - beta)))/(Pi2*v2) + 
 (0.015625*Re(B0i(bb0, Mh2, MW2, MW2))*(4*(-Mh2 + MHH2)*MW2 + 
    (Mh2*(MHH2 - 4*MW2) + 2*MW2*(MHH2 + 6*MW2) + Cos(2*(alp - beta))*
       (Mh2*MHH2 - 2*MHH2*MW2 + 12*(pow(MW, 4))))*Sin(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.75*(pow(MB, 4))*Re(B0i(dbb0, Mh2, MB2, MB2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.75*(pow(MC, 4))*Re(B0i(dbb0, Mh2, MC2, MC2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.75*(pow(MD, 4))*Re(B0i(dbb0, Mh2, MD2, MD2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.75*(pow(MS, 4))*Re(B0i(dbb0, Mh2, MS2, MS2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.75*(pow(MT, 4))*Re(B0i(dbb0, Mh2, MT2, MT2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.75*(pow(MU, 4))*Re(B0i(dbb0, Mh2, MU2, MU2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.375*Mh2*(pow(MB, 2))*Re(B0i(dbb1, Mh2, MB2, MB2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.375*Mh2*(pow(MC, 2))*Re(B0i(dbb1, Mh2, MC2, MC2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.375*Mh2*(pow(MD, 2))*Re(B0i(dbb1, Mh2, MD2, MD2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.375*Mh2*(pow(MS, 2))*Re(B0i(dbb1, Mh2, MS2, MS2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.375*Mh2*MT2*Re(B0i(dbb1, Mh2, MT2, MT2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.375*Mh2*(pow(MU, 2))*Re(B0i(dbb1, Mh2, MU2, MU2))*
   (1 + (pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sin(alp - beta)))/
  (Pi2*v2) + (0.00048828125*(pow(Sec(w), 2))*Re(B0i(bb0, Mh2, MZ2, MZ2))*
   (-64*MW2 + (pow(Mh2 - MHH2, -1))*(6*Mh2*MHH2 - 32*Mh2*MW2 + 
      16*MHH2*MW2 + Mh2*MHH2*Cos(2*(alp - beta - 2*w)) + 
      4*Mh2*MHH2*Cos(2*(alp - beta - w)) - 8*MHH2*MW2*
       Cos(2*(alp - beta - w)) + 8*Mh2*MHH2*Cos(2*w) - 32*Mh2*MW2*Cos(2*w) + 
      16*MHH2*MW2*Cos(2*w) + 2*Mh2*MHH2*Cos(4*w) + 
      4*Mh2*MHH2*Cos(2*(alp - beta + w)) - 8*MHH2*MW2*
       Cos(2*(alp - beta + w)) + Mh2*MHH2*Cos(2*(alp - beta + 2*w)) + 
      192*(pow(MW, 4)) + 2*Cos(2*(alp - beta))*(3*Mh2*MHH2 - 8*MHH2*MW2 + 
        96*(pow(MW, 4))))*(pow(Sec(w), 2))*Sin(alp - beta)))/
  (Pi2*v2) + (0.0625*MW2*(pow(Cos(w), 2))*Re(B0i(bb1, Mh2, MW2, MW2))*
   (8*(pow(Csc(2*w), 2)) + (2*Mh2 - MHH2 + MHH2*Cos(2*(alp - beta)))*
     (pow(Mh2 - MHH2, -1))*(pow(Csc(w), 2))*(pow(Sec(w), 2))*
     Sin(alp - beta)))/(Pi2*v2*(pow(Csc(w), 2))) + 
 (0.03515625*(pow(Cos(w), 2))*Re(B0i(dbb0, Mh2, Mh2, Mh2))*
   (16*(pow(Mh, 4))*(pow(Csc(2*w), 2)) + 
    (pow(M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
        (2*M2 - 3*Mh2)*Cos(alp + beta), 2))*(pow(Csc(2*beta), 2))*
     (pow(Csc(w), 2))*(pow(Sec(w), 2))*Sin(alp - beta)))/
  (Pi2*v2*(pow(Csc(w), 2))) + 
 (0.03125*MW2*(pow(Cos(w), 2))*Re(B0i(bb1, Mh2, MZ2, MZ2))*
   (32*(pow(Csc(2*w), 4)) + (2*Mh2 - MHH2 + MHH2*Cos(2*(alp - beta)))*
     (pow(Mh2 - MHH2, -1))*(pow(Csc(w), 4))*(pow(Sec(w), 4))*
     Sin(alp - beta)))/(Pi2*v2*(pow(Csc(w), 4))) - 
 (0.25*(pow(ME, 4))*Re(B0i(dbb0, Mh2, ME2, ME2))*
   (1 + (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.25*(pow(ML, 4))*Re(B0i(dbb0, Mh2, ML2, ML2))*
   (1 + (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.25*(pow(MM, 4))*Re(B0i(dbb0, Mh2, MM2, MM2))*
   (1 + (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.125*Mh2*(pow(ME, 2))*Re(B0i(dbb1, Mh2, ME2, ME2))*
   (1 + (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.125*Mh2*(pow(ML, 2))*Re(B0i(dbb1, Mh2, ML2, ML2))*
   (1 + (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*Sin(alp - beta)))/
  (Pi2*v2) - (0.125*Mh2*(pow(MM, 2))*Re(B0i(dbb1, Mh2, MM2, MM2))*
   (1 + (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*Sin(alp - beta)))/
  (Pi2*v2) + (0.015625*Cos(alp - beta)*((MA02 - Mh2)*(MA02 - MHH2) - 
    (MA02 + MHH2)*MW2*(pow(Sec(w), 2)))*Re(B0i(bb0, MHH2, MA02, MZ2))*
   Sin(2*(alp - beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MB, 2))*Re(B0i(bb1, Mh2, MB2, MB2))*
   (-2*Mh2 + 2*MHH2 + Cos(alp)*Csc(beta)*(2*Mh2 - MHH2 + 
      MHH2*Csc(beta)*Sin(2*alp - beta))))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MC, 2))*Re(B0i(bb1, Mh2, MC2, MC2))*
   (-2*Mh2 + 2*MHH2 + Cos(alp)*Csc(beta)*(2*Mh2 - MHH2 + 
      MHH2*Csc(beta)*Sin(2*alp - beta))))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MD, 2))*Re(B0i(bb1, Mh2, MD2, MD2))*
   (-2*Mh2 + 2*MHH2 + Cos(alp)*Csc(beta)*(2*Mh2 - MHH2 + 
      MHH2*Csc(beta)*Sin(2*alp - beta))))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MS, 2))*Re(B0i(bb1, Mh2, MS2, MS2))*
   (-2*Mh2 + 2*MHH2 + Cos(alp)*Csc(beta)*(2*Mh2 - MHH2 + 
      MHH2*Csc(beta)*Sin(2*alp - beta))))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*MT2*Re(B0i(bb1, Mh2, MT2, MT2))*(-2*Mh2 + 2*MHH2 + 
    Cos(alp)*Csc(beta)*(2*Mh2 - MHH2 + MHH2*Csc(beta)*Sin(2*alp - beta))))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.1875*(pow(MU, 2))*
   Re(B0i(bb1, Mh2, MU2, MU2))*(-2*Mh2 + 2*MHH2 + 
    Cos(alp)*Csc(beta)*(2*Mh2 - MHH2 + MHH2*Csc(beta)*Sin(2*alp - beta))))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.08333333333333333*MT2*B0i(bb0, m32, MT2, MT2)*
   (9 - 2*Cos(2*w) + 2*Cos(4*w))*Csc(beta)*(Cos(alp) - Sin(beta)))/(Pi2*v2) + 
 (0.08333333333333333*B0i(bb0, m32, MC2, MC2)*(9 - 2*Cos(2*w) + 2*Cos(4*w))*
   Csc(beta)*(pow(MC, 2))*(Cos(alp) - Sin(beta)))/(Pi2*v2) + 
 (0.08333333333333333*B0i(bb0, m32, MU2, MU2)*(9 - 2*Cos(2*w) + 2*Cos(4*w))*
   Csc(beta)*(pow(MU, 2))*(Cos(alp) - Sin(beta)))/(Pi2*v2) + 
 (0.09375*B0i(bb0, m12, Mh2, Mh2)*Csc(2*beta)*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   (2*(M2 - Mh2)*Cos(alp + beta) + (-M2 + Mh2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (Pi2*v2) - (0.125*C0i(cc00, m22, m12, m32, MA02, MHH2, MHH2)*Csc(2*beta)*
   (pow(Sin(alp - beta), 3))*((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(Pi2*v2) + (0.03125*B0i(bb0, m12, MHH2, MHH2)*
   Csc(2*beta)*Sin(alp - beta)*((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/(Pi2*v2) - 
 (0.0625*C0i(cc00, m12, m32, m22, MHH2, MHH2, MZ2)*Csc(beta)*
   (pow(Cos(alp - beta), 2))*Sec(beta)*Sin(alp - beta)*
   ((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) + 
 (0.0625*MW2*C0i(cc0, m12, m32, m22, MHH2, MHH2, MZ2)*Csc(beta)*
   (pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*Sec(beta)*
   Sin(alp - beta)*((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (Pi2*v2) + (0.0234375*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh2)*Cos(3*alp - beta) + (2*M2 - 3*Mh2)*Cos(alp + beta))*
   (pow(Cos(alp - beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(bb0, Mh2, Mh2, Mh2))*((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0234375*(M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh2)*Cos(alp + beta))*(pow(Cos(alp - beta), 2))*
   (pow(Csc(2*beta), 2))*Re(B0i(bb0, MHH2, Mh2, Mh2))*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.0625*C0i(cc00, m12, m22, m32, Mh2, MHH2, MZ2)*
   Csc(beta)*(pow(Cos(alp - beta), 2))*Sec(beta)*Sin(alp - beta)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) - 
 (0.0625*C0i(cc00, m12, m32, m22, Mh2, MHH2, MZ2)*Csc(beta)*
   (pow(Cos(alp - beta), 2))*Sec(beta)*Sin(alp - beta)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) + 
 (0.0625*C0i(cc00, m22, m12, m32, MA02, Mh2, MHH2)*Csc(beta)*
   (pow(Cos(alp - beta), 2))*Sec(beta)*Sin(alp - beta)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) + 
 (0.0625*C0i(cc00, m32, m12, m22, MA02, Mh2, MHH2)*Csc(beta)*
   (pow(Cos(alp - beta), 2))*Sec(beta)*Sin(alp - beta)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) + 
 (0.0625*MW2*C0i(cc0, m12, m22, m32, Mh2, MHH2, MZ2)*Csc(beta)*
   (pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*Sec(beta)*
   Sin(alp - beta)*((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (Pi2*v2) + (0.0625*MW2*C0i(cc0, m12, m32, m22, Mh2, MHH2, MZ2)*Csc(beta)*
   (pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*Sec(beta)*
   Sin(alp - beta)*((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (Pi2*v2) + (0.0078125*Csc(2*beta)*Re(A0i(aa0, MHH2))*Sin(alp - beta)*
   (2*(-5*M2 + Mh2 + 5*MHH2)*Sin(2*alp) - 2*M2*Sin(2*(alp - 2*beta)) + 
    Mh2*Sin(4*alp - 2*beta) - MHH2*Sin(4*alp - 2*beta) - 8*M2*Sin(2*beta) - 
    4*MA02*Sin(2*beta) + Mh2*Sin(2*beta) + 11*MHH2*Sin(2*beta)))/
  (MA02*Pi2*v2) - (0.015625*(pow(Cos(alp - beta), 2))*
   (pow(Csc(2*beta), 2))*Re(B0i(bb0, Mh2, Mh2, MHH2))*Sin(alp - beta)*
   (-9*M2*Mh2 - 9*M2*MHH2 + 5*Mh2*MHH2 + 8*(pow(M2, 2)) + 
    2*(pow(Mh, 4)) + 2*(pow(MHH2, 2)) - 
    Cos(4*alp)*(5*Mh2*MHH2 - 9*M2*(Mh2 + MHH2) + 9*(pow(M2, 2)) + 
      2*(pow(Mh, 4)) + 2*(pow(MHH2, 2))) + 
    (pow(M2, 2))*(pow(Cos(beta), 4)) - 6*(pow(M2, 2))*
     (pow(Cos(beta), 2))*(pow(Sin(beta), 2)) + 
    (pow(M2, 2))*(pow(Sin(beta), 4)) - 2*M2*(Mh2 - MHH2)*Sin(2*alp)*
     Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.015625*(pow(Cos(alp - beta), 2))*(pow(Csc(2*beta), 2))*
   Re(B0i(bb0, MHH2, Mh2, MHH2))*Sin(alp - beta)*
   (-9*M2*Mh2 - 9*M2*MHH2 + 5*Mh2*MHH2 + 8*(pow(M2, 2)) + 
    2*(pow(Mh, 4)) + 2*(pow(MHH2, 2)) - 
    Cos(4*alp)*(5*Mh2*MHH2 - 9*M2*(Mh2 + MHH2) + 9*(pow(M2, 2)) + 
      2*(pow(Mh, 4)) + 2*(pow(MHH2, 2))) + 
    (pow(M2, 2))*(pow(Cos(beta), 4)) - 6*(pow(M2, 2))*
     (pow(Cos(beta), 2))*(pow(Sin(beta), 2)) + 
    (pow(M2, 2))*(pow(Sin(beta), 4)) - 2*M2*(Mh2 - MHH2)*Sin(2*alp)*
     Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.01171875*(pow(Csc(2*beta), 2))*Re(B0i(bb0, Mh2, MHH2, MHH2))*
   Sin(2*(alp - beta))*((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta))*
   (-(M2*Sin(alp - 3*beta)) + (M2 - MHH2)*Sin(3*alp - beta) + 
    (-2*M2 + 3*MHH2)*Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.01171875*(pow(Csc(2*beta), 2))*Re(B0i(bb0, MHH2, MHH2, MHH2))*
   Sin(2*(alp - beta))*((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta))*
   (-(M2*Sin(alp - 3*beta)) + (M2 - MHH2)*Sin(3*alp - beta) + 
    (-2*M2 + 3*MHH2)*Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0009765625*Cos(alp - beta)*((2*MA02 - Mh2)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh2, MA02, MA02))*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0009765625*Cos(alp - beta)*((-2*MA02 + Mh2)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MA02, MA02))*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0078125*(MA02 - MHH2)*B0i(bb0, 0, MA02, MHH2)*Csc(2*beta)*
   Sin(2*(alp - beta))*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/(MA02*Pi2*v2) + 
 (0.0078125*(MA02 - MHH2)*B0i(bb0, MA02, MA02, MHH2)*Csc(2*beta)*
   Sin(2*(alp - beta))*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/(MA02*Pi2*v2) - 
 (0.0078125*A0i(aa0, MHH2)*Csc(2*beta)*Sin(2*(alp - beta))*
   (-2*(M2 - MA02)*Sin(alp - 3*beta) + (Mh2 - MHH2)*Sin(3*alp - beta) + 
    (-2*M2 - 2*MA02 + Mh2 + 3*MHH2)*Sin(alp + beta)))/(MA02*Pi2*v2) - 
 (0.001953125*Cos(alp - beta)*((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh2 + 2*MHp2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh2, MHp2, MHp2))*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.001953125*Cos(alp - beta)*((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh2 + 2*MHp2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MHp2, MHp2))*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.125*MW2*Csc(2*w)*Re(B0i(bb1, MW2, 0, MW2))*(1 + Sin(alp - beta))*
   Sin(4*w))/(Pi2*v2) - (0.25*A0i(aa0, ME2)*Cos(alp - beta)*(pow(ME, 2))*
   Tan(beta))/(MA02*Pi2*v2) - (0.125*B0i(bb1, MA02, ME2, ME2)*Cos(alp - beta)*
   (pow(ME, 2))*Tan(beta))/(Pi2*v2) - 
 (0.25*A0i(aa0, ML2)*Cos(alp - beta)*(pow(ML, 2))*Tan(beta))/
  (MA02*Pi2*v2) - (0.125*B0i(bb1, MA02, ML2, ML2)*Cos(alp - beta)*
   (pow(ML, 2))*Tan(beta))/(Pi2*v2) - 
 (0.25*A0i(aa0, MM2)*Cos(alp - beta)*(pow(MM, 2))*Tan(beta))/
  (MA02*Pi2*v2) - (0.125*B0i(bb1, MA02, MM2, MM2)*Cos(alp - beta)*
   (pow(MM, 2))*Tan(beta))/(Pi2*v2) - 
 (0.25*Cos(alp)*(pow(ME, 4))*Re(B0i(bb0, Mh2, ME2, ME2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.25*Cos(alp)*(pow(ML, 4))*Re(B0i(bb0, Mh2, ML2, ML2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.25*Cos(alp)*(pow(MM, 4))*Re(B0i(bb0, Mh2, MM2, MM2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*Cos(alp)*(pow(ME, 4))*Re(B0i(bb0, MHH2, ME2, ME2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*Cos(alp)*(pow(ML, 4))*Re(B0i(bb0, MHH2, ML2, ML2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*Cos(alp)*(pow(MM, 4))*Re(B0i(bb0, MHH2, MM2, MM2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.125*MHH2*Cos(alp)*(pow(ME, 2))*Re(B0i(bb1, MHH2, ME2, ME2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.125*MHH2*Cos(alp)*(pow(ML, 2))*Re(B0i(bb1, MHH2, ML2, ML2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.125*MHH2*Cos(alp)*(pow(MM, 2))*Re(B0i(bb1, MHH2, MM2, MM2))*Sec(beta)*
   Sin(alp)*(Cos(alp) + Sin(alp)*Tan(beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.015625*Cos(w)*(pow(Sin(w), 3))*Re(A0i(aa0, ME2))*
   (-64*MA02*Cos(alp)*Csc(w)*Sec(w)*Sin(beta) + MA02*(pow(Csc(w), 5))*
     Sec(w)*(4 + Cos(6*w) + 4*Cos(beta)*Sin(alp) + Cos(beta)*Cos(6*w)*
       Sin(alp) - 11*Cos(2*w)*(1 + Cos(beta)*Sin(alp)) + 
      2*Cos(4*w)*(1 + Cos(beta)*Sin(alp)) + 4*Cos(alp)*Sin(beta)) + 
    2*(pow(Csc(w), 3))*Sec(w)*(2*Cos(alp)*Cos(beta)*
       (3*MA02 - 4*(pow(ME, 2))) - 8*(pow(ME, 2))*Sin(alp)*Sin(beta))*
     Tan(beta) + 32*MA02*Cos(alp)*Sin(beta)*Tan(w)))/(MA02*Pi2*v2) - 
 (0.015625*Cos(w)*(pow(Sin(w), 3))*Re(A0i(aa0, ML2))*
   (-64*MA02*Cos(alp)*Csc(w)*Sec(w)*Sin(beta) + MA02*(pow(Csc(w), 5))*
     Sec(w)*(4 + Cos(6*w) + 4*Cos(beta)*Sin(alp) + Cos(beta)*Cos(6*w)*
       Sin(alp) - 11*Cos(2*w)*(1 + Cos(beta)*Sin(alp)) + 
      2*Cos(4*w)*(1 + Cos(beta)*Sin(alp)) + 4*Cos(alp)*Sin(beta)) + 
    2*(pow(Csc(w), 3))*Sec(w)*(2*Cos(alp)*Cos(beta)*
       (3*MA02 - 4*(pow(ML, 2))) - 8*(pow(ML, 2))*Sin(alp)*Sin(beta))*
     Tan(beta) + 32*MA02*Cos(alp)*Sin(beta)*Tan(w)))/(MA02*Pi2*v2) - 
 (0.015625*Cos(w)*(pow(Sin(w), 3))*Re(A0i(aa0, MM2))*
   (-64*MA02*Cos(alp)*Csc(w)*Sec(w)*Sin(beta) + MA02*(pow(Csc(w), 5))*
     Sec(w)*(4 + Cos(6*w) + 4*Cos(beta)*Sin(alp) + Cos(beta)*Cos(6*w)*
       Sin(alp) - 11*Cos(2*w)*(1 + Cos(beta)*Sin(alp)) + 
      2*Cos(4*w)*(1 + Cos(beta)*Sin(alp)) + 4*Cos(alp)*Sin(beta)) + 
    2*(pow(Csc(w), 3))*Sec(w)*(2*Cos(alp)*Cos(beta)*
       (3*MA02 - 4*(pow(MM, 2))) - 8*(pow(MM, 2))*Sin(alp)*Sin(beta))*
     Tan(beta) + 32*MA02*Cos(alp)*Sin(beta)*Tan(w)))/(MA02*Pi2*v2) + 
 (0.010416666666666666*Cos(w)*(pow(Sin(w), 3))*Re(A0i(aa0, MT2))*
   (-((7 + Cos(2*w) + 8*Cos(4*w) + 2*Cos(6*w))*(pow(Csc(w), 5))*Sec(w)) + 
    2*(pow(MA02, -1))*(-36*MT2*(pow(Cos(alp), 3))*
       (pow(Cot(beta), 2))*(pow(Csc(w), 4)) + 
      Cos(alp)*(-32*MA02 + 80*MA02*(pow(Csc(w), 2)) + 
        9*MA02*(pow(Csc(w), 6)) - 3*(pow(Csc(w), 4))*
         (17*MA02 + 12*MT2*(pow(Cot(beta), 2))*(pow(Sin(alp), 2)))) - 
      36*MT2*Cot(beta)*(pow(Cos(alp), 2))*(pow(Csc(w), 4))*Sin(alp) + 
      Cot(beta)*(32*MA02 - 80*MA02*(pow(Csc(w), 2)) - 
        9*MA02*(pow(Csc(w), 6)) + (pow(Csc(w), 4))*
         (51*MA02 - 36*MT2*(pow(Sin(alp), 2))))*Sin(alp))*Sin(beta)*
     Tan(w)))/(Pi2*v2) + (0.005208333333333333*Cos(w)*(pow(Sin(w), 3))*
   Re(A0i(aa0, MB2))*(-((-4 - 43*Cos(2*w) + 10*Cos(4*w) + Cos(6*w))*
      (pow(Csc(w), 5))*Sec(w)) - 4*(pow(MA02, -1))*
     (36*(pow(MB, 2))*(pow(Cos(alp), 3))*(pow(Cot(beta), 2))*
       (pow(Csc(w), 4)) + Cos(alp)*(8*MA02 - 
        32*MA02*(pow(Csc(w), 2)) + 9*MA02*(pow(Csc(w), 6)) + 
        3*(pow(Csc(w), 4))*(MA02 + 12*(pow(MB, 2))*
           (pow(Cot(beta), 2))*(pow(Sin(alp), 2)))) + 
      36*Cot(beta)*(pow(MB, 2))*(pow(Cos(alp), 2))*
       (pow(Csc(w), 4))*Sin(alp) + Cot(beta)*
       (-8*MA02 + 32*MA02*(pow(Csc(w), 2)) - 
        9*MA02*(pow(Csc(w), 6)) - 3*(pow(Csc(w), 4))*
         (MA02 - 12*(pow(MB, 2))*(pow(Sin(alp), 2))))*Sin(alp))*
     Sin(beta)*Tan(w)))/(Pi2*v2) + 
 (0.010416666666666666*Cos(w)*(pow(Sin(w), 3))*Re(A0i(aa0, MC2))*
   (-((7 + Cos(2*w) + 8*Cos(4*w) + 2*Cos(6*w))*(pow(Csc(w), 5))*Sec(w)) + 
    2*(pow(MA02, -1))*(-36*(pow(MC, 2))*(pow(Cos(alp), 3))*
       (pow(Cot(beta), 2))*(pow(Csc(w), 4)) + 
      Cos(alp)*(-32*MA02 + 80*MA02*(pow(Csc(w), 2)) + 
        9*MA02*(pow(Csc(w), 6)) - 3*(pow(Csc(w), 4))*
         (17*MA02 + 12*(pow(MC, 2))*(pow(Cot(beta), 2))*
           (pow(Sin(alp), 2)))) - 36*Cot(beta)*(pow(MC, 2))*
       (pow(Cos(alp), 2))*(pow(Csc(w), 4))*Sin(alp) + 
      Cot(beta)*(32*MA02 - 80*MA02*(pow(Csc(w), 2)) - 
        9*MA02*(pow(Csc(w), 6)) + (pow(Csc(w), 4))*
         (51*MA02 - 36*(pow(MC, 2))*(pow(Sin(alp), 2))))*Sin(alp))*
     Sin(beta)*Tan(w)))/(Pi2*v2) + 
 (0.005208333333333333*Cos(w)*(pow(Sin(w), 3))*Re(A0i(aa0, MD2))*
   (-((-4 - 43*Cos(2*w) + 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 5))*
      Sec(w)) - 4*(pow(MA02, -1))*
     (36*(pow(MD, 2))*(pow(Cos(alp), 3))*(pow(Cot(beta), 2))*
       (pow(Csc(w), 4)) + Cos(alp)*(8*MA02 - 
        32*MA02*(pow(Csc(w), 2)) + 9*MA02*(pow(Csc(w), 6)) + 
        3*(pow(Csc(w), 4))*(MA02 + 12*(pow(MD, 2))*
           (pow(Cot(beta), 2))*(pow(Sin(alp), 2)))) + 
      36*Cot(beta)*(pow(MD, 2))*(pow(Cos(alp), 2))*
       (pow(Csc(w), 4))*Sin(alp) + Cot(beta)*
       (-8*MA02 + 32*MA02*(pow(Csc(w), 2)) - 
        9*MA02*(pow(Csc(w), 6)) - 3*(pow(Csc(w), 4))*
         (MA02 - 12*(pow(MD, 2))*(pow(Sin(alp), 2))))*Sin(alp))*
     Sin(beta)*Tan(w)))/(Pi2*v2) + 
 (0.005208333333333333*Cos(w)*(pow(Sin(w), 3))*Re(A0i(aa0, MS2))*
   (-((-4 - 43*Cos(2*w) + 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 5))*
      Sec(w)) - 4*(pow(MA02, -1))*
     (36*(pow(MS, 2))*(pow(Cos(alp), 3))*(pow(Cot(beta), 2))*
       (pow(Csc(w), 4)) + Cos(alp)*(8*MA02 - 
        32*MA02*(pow(Csc(w), 2)) + 9*MA02*(pow(Csc(w), 6)) + 
        3*(pow(Csc(w), 4))*(MA02 + 12*(pow(MS, 2))*
           (pow(Cot(beta), 2))*(pow(Sin(alp), 2)))) + 
      36*Cot(beta)*(pow(MS, 2))*(pow(Cos(alp), 2))*
       (pow(Csc(w), 4))*Sin(alp) + Cot(beta)*
       (-8*MA02 + 32*MA02*(pow(Csc(w), 2)) - 
        9*MA02*(pow(Csc(w), 6)) - 3*(pow(Csc(w), 4))*
         (MA02 - 12*(pow(MS, 2))*(pow(Sin(alp), 2))))*Sin(alp))*
     Sin(beta)*Tan(w)))/(Pi2*v2) + 
 (0.010416666666666666*Cos(w)*(pow(Sin(w), 3))*Re(A0i(aa0, MU2))*
   (-((7 + Cos(2*w) + 8*Cos(4*w) + 2*Cos(6*w))*(pow(Csc(w), 5))*Sec(w)) + 
    2*(pow(MA02, -1))*(-36*(pow(MU, 2))*(pow(Cos(alp), 3))*
       (pow(Cot(beta), 2))*(pow(Csc(w), 4)) + 
      Cos(alp)*(-32*MA02 + 80*MA02*(pow(Csc(w), 2)) + 
        9*MA02*(pow(Csc(w), 6)) - 3*(pow(Csc(w), 4))*
         (17*MA02 + 12*(pow(MU, 2))*(pow(Cot(beta), 2))*
           (pow(Sin(alp), 2)))) - 36*Cot(beta)*(pow(MU, 2))*
       (pow(Cos(alp), 2))*(pow(Csc(w), 4))*Sin(alp) + 
      Cot(beta)*(32*MA02 - 80*MA02*(pow(Csc(w), 2)) - 
        9*MA02*(pow(Csc(w), 6)) + (pow(Csc(w), 4))*
         (51*MA02 - 36*(pow(MU, 2))*(pow(Sin(alp), 2))))*Sin(alp))*
     Sin(beta)*Tan(w)))/(Pi2*v2)
;
}

ComplexType dKappahbb(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32)
{
return 0. + (0.03125*(4*Mh2*MT2 - (m12 + m22 - m32)*MW2 + MHp2*(-4*MT2 + 2*MW2))*
   C0i(cc0, m22, m32, m12, MHp2, MT2, MW2)*Cos(alp - beta)*Cot(beta))/
  (Pi2*v2) + (0.03125*(4*Mh2*MT2 - 4*MHp2*MT2 + 
    MW2*(-3*m12 + m22 - m32 + 2*MW2))*C0i(cc0, m32, m22, m12, MHp2, MT2, MW2)*
   Cos(alp - beta)*Cot(beta))/(Pi2*v2) - 
 (0.03125*(m12 - m22 - m32)*MW2*C0i(cc1, m22, m32, m12, MHp2, MT2, MW2)*
   Cos(alp - beta)*Cot(beta))/(Pi2*v2) - 
 (0.03125*(3*m12 - 3*m22 + m32)*MW2*C0i(cc1, m32, m22, m12, MHp2, MT2, MW2)*
   Cos(alp - beta)*Cot(beta))/(Pi2*v2) - 
 (0.03125*(m12 - m22 + m32)*MW2*C0i(cc2, m22, m32, m12, MHp2, MT2, MW2)*
   Cos(alp - beta)*Cot(beta))/(Pi2*v2) - 
 (0.03125*(5*m12 + m22 - m32)*MW2*C0i(cc2, m32, m22, m12, MHp2, MT2, MW2)*
   Cos(alp - beta)*Cot(beta))/(Pi2*v2) + 1.*(-1 + Cos(alp)*Csc(beta)) + 
 (0.125*m12*MT2*C0i(cc1, m12, m32, m22, MT2, MT2, MW2)*
   (-1 + Cos(alp)*Csc(beta)))/(Pi2*v2) + 
 (0.0625*(m12 + m22 - m32)*MT2*C0i(cc2, m12, m32, m22, MT2, MT2, MW2)*
   (-1 + Cos(alp)*Csc(beta)))/(Pi2*v2) + 
 (0.03125*B0i(bb0, m32, MT2, MW2)*(-4*(MT2 + MW2) + 
    (4*MT2 + 3*MW2)*Cos(alp)*Csc(beta) - MW2*Cos(alp - 2*beta)*Csc(beta)))/
  (Pi2*v2) + (0.25*A0i(aa0, ME2)*Cos(alp)*Csc(beta)*(pow(ME, 2)))/
  (MA02*Pi2*v2) + (0.125*B0i(bb1, MA02, ME2, ME2)*Cos(alp)*Csc(beta)*
   (pow(ME, 2)))/(Pi2*v2) + (0.25*A0i(aa0, ML2)*Cos(alp)*Csc(beta)*
   (pow(ML, 2)))/(MA02*Pi2*v2) + 
 (0.125*B0i(bb1, MA02, ML2, ML2)*Cos(alp)*Csc(beta)*(pow(ML, 2)))/
  (Pi2*v2) + (0.25*A0i(aa0, MM2)*Cos(alp)*Csc(beta)*(pow(MM, 2)))/
  (MA02*Pi2*v2) + (0.125*B0i(bb1, MA02, MM2, MM2)*Cos(alp)*Csc(beta)*
   (pow(MM, 2)))/(Pi2*v2) + (0.25*C0i(cc0, m12, m32, m22, MT2, MT2, MW2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MT, 4)))/(Pi2*v2) - 
 (0.75*MT2*A0i(aa0, MT2)*Cos(alp)*Csc(beta)*(pow(Cot(beta), 2)))/
  (MA02*Pi2*v2) - (0.375*MT2*B0i(bb1, MA02, MT2, MT2)*Cos(alp)*Csc(beta)*
   (pow(Cot(beta), 2)))/(Pi2*v2) + 
 (0.0625*MT2*(-m12 - m22 + m32 + 4*MT2)*C0i(cc0, m22, m12, m32, MHp2, MT2, 
    MT2)*Cos(alp)*Csc(beta)*(pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.0625*(m12 + m22 - m32)*MT2*C0i(cc1, m22, m12, m32, MHp2, MT2, MT2)*
   Cos(alp)*Csc(beta)*(pow(Cot(beta), 2)))/(Pi2*v2) + 
 (0.0625*(m12 - m22 + m32)*MT2*C0i(cc2, m22, m12, m32, MHp2, MT2, MT2)*
   Cos(alp)*Csc(beta)*(pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MB2)*Cos(alp)*Csc(beta)*(pow(MB, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MB2, MB2)*Cos(alp)*Csc(beta)*(pow(MB, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.03125*(m12 + m22 - m32)*C0i(cc1, m22, m12, m32, MA02, MB2, MB2)*Cos(alp)*
   Csc(beta)*(pow(MB, 2))*(pow(Cot(beta), 2)))/(Pi2*v2) + 
 (0.03125*(m12 - m22 + m32)*C0i(cc2, m22, m12, m32, MA02, MB2, MB2)*Cos(alp)*
   Csc(beta)*(pow(MB, 2))*(pow(Cot(beta), 2)))/(Pi2*v2) + 
 (0.03125*C0i(cc0, m22, m12, m32, MA02, MB2, MB2)*Cos(alp)*Csc(beta)*
   (pow(MB, 2))*(-m12 - m22 + m32 + 4*(pow(MB, 2)))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MC2)*Cos(alp)*Csc(beta)*(pow(MC, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MC2, MC2)*Cos(alp)*Csc(beta)*(pow(MC, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MD2)*Cos(alp)*Csc(beta)*(pow(MD, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MD2, MD2)*Cos(alp)*Csc(beta)*(pow(MD, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MS2)*Cos(alp)*Csc(beta)*(pow(MS, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MS2, MS2)*Cos(alp)*Csc(beta)*(pow(MS, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MU2)*Cos(alp)*Csc(beta)*(pow(MU, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MU2, MU2)*Cos(alp)*Csc(beta)*(pow(MU, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) + 
 (0.03125*MT2*C0i(cc0, m12, m32, m22, MHp2, MHp2, MT2)*
   ((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*Cot(beta)*(pow(Csc(beta), 2)))/(Pi2*v2) - 
 (0.015625*C0i(cc0, m12, m32, m22, MA02, MA02, MB2)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*Cot(beta)*(pow(MB, 2))*(pow(Csc(beta), 2)))/
  (Pi2*v2) + (0.01171875*A0i(aa0, MA02)*Cos(alp)*(3*(Mh2 - MHH2)*Cos(2*alp) + 
    (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 4*(-2*M2 + Mh2 + MHH2)*Cos(2*beta))*
   (pow(Csc(beta), 3)))/(MA02*Pi2*v2) + 
 (0.0078125*A0i(aa0, MHp2)*Cos(alp)*(3*(Mh2 - MHH2)*Cos(2*alp) + 
    (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 4*(-2*M2 + Mh2 + MHH2)*Cos(2*beta))*
   (pow(Csc(beta), 3)))/(MA02*Pi2*v2) + 
 (0.0078125*(MA02 - Mh2)*B0i(bb0, 0, MA02, Mh2)*Cos(alp)*Cos(alp - beta)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*(pow(Csc(beta), 3)))/(MA02*Pi2*v2) + 
 (0.0078125*(MA02 - Mh2)*B0i(bb0, MA02, MA02, Mh2)*Cos(alp)*Cos(alp - beta)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*(pow(Csc(beta), 3)))/(MA02*Pi2*v2) + 
 (0.0078125*A0i(aa0, Mh2)*Cos(alp)*Cos(alp - beta)*
   (-2*(M2 - MA02)*Cos(alp - 3*beta) + (Mh2 - MHH2)*Cos(3*alp - beta) + 
    (-2*M2 - 2*MA02 + 3*Mh2 + MHH2)*Cos(alp + beta))*(pow(Csc(beta), 3)))/
  (MA02*Pi2*v2) - (0.0078125*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MB2, 
    MB2, Mh2)*(pow(MB, 2))*(-4 + 3*Cos(alp)*(pow(Csc(beta), 3)) + 
    Cos(3*alp)*(pow(Csc(beta), 3))))/(Pi2*v2) - 
 (0.0625*B0i(bb0, m32, MB2, Mh2)*(pow(MB, 2))*
   (-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3))))/(Pi2*v2) - 
 (0.0625*m12*C0i(cc1, m12, m32, m22, MB2, MB2, Mh2)*(pow(MB, 2))*
   (-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3))))/(Pi2*v2) - 
 (0.125*C0i(cc0, m12, m32, m22, MB2, MB2, Mh2)*(pow(MB, 4))*
   (-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3))))/(Pi2*v2) - 
 (0.015625*(m12 - m22 - m32)*MW2*C0i(cc1, m22, m32, m12, MA02, MB2, MZ2)*
   Cos(alp - beta)*Cot(beta)*(pow(Sec(w), 2)))/(Pi2*v2) - 
 (0.015625*(3*m12 - 3*m22 + m32)*MW2*C0i(cc1, m32, m22, m12, MA02, MB2, MZ2)*
   Cos(alp - beta)*Cot(beta)*(pow(Sec(w), 2)))/(Pi2*v2) - 
 (0.015625*(m12 - m22 + m32)*MW2*C0i(cc2, m22, m32, m12, MA02, MB2, MZ2)*
   Cos(alp - beta)*Cot(beta)*(pow(Sec(w), 2)))/(Pi2*v2) - 
 (0.015625*(5*m12 + m22 - m32)*MW2*C0i(cc2, m32, m22, m12, MA02, MB2, MZ2)*
   Cos(alp - beta)*Cot(beta)*(pow(Sec(w), 2)))/(Pi2*v2) - 
 (0.015625*C0i(cc0, m22, m32, m12, MA02, MB2, MZ2)*Cos(alp - beta)*Cot(beta)*
   (m12*MW2 + m22*MW2 - m32*MW2 - 2*Mh2*(pow(MB, 2)) + 
    2*(MA02 - Mh2)*Cos(2*w)*(pow(MB, 2)) + 
    2*MA02*(-MW2 + (pow(MB, 2))))*(pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.003472222222222222*m12*C0i(cc1, m12, m32, m22, MB2, MB2, MZ2)*
   (-1 + Cos(alp)*Csc(beta))*(-12*MW2 + 4*MW2*Cos(4*w) + 9*(pow(MB, 2)) + 
    Cos(2*w)*(8*MW2 + 9*(pow(MB, 2))))*(pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.001736111111111111*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MB2, MB2, 
    MZ2)*(-1 + Cos(alp)*Csc(beta))*(-12*MW2 + 4*MW2*Cos(4*w) + 
    9*(pow(MB, 2)) + Cos(2*w)*(8*MW2 + 9*(pow(MB, 2))))*
   (pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.006944444444444444*C0i(cc0, m12, m32, m22, MB2, MB2, MZ2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*(-12*MW2 + 4*MW2*Cos(4*w) + 
    9*(pow(MB, 2)) + Cos(2*w)*(8*MW2 + 9*(pow(MB, 2))))*
   (pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.001736111111111111*B0i(bb0, m32, MB2, MZ2)*
   (Csc(beta)*(-9*MW2*Cos(alp - 2*beta) + Cos(alp)*(3*MW2 + 8*MW2*Cos(4*w) + 
        18*(pow(MB, 2)) + 2*Cos(2*w)*(8*MW2 + 9*(pow(MB, 2)))))*
     (pow(Sec(w), 2)) - 4*(9*(pow(MB, 2)) + 
      MW2*(pow(1 + 2*Cos(2*w), 2))*(pow(Sec(w), 2)))))/(Pi2*v2) - 
 (0.0078125*C0i(cc0, m32, m22, m12, MA02, MB2, MZ2)*Cos(alp - beta)*Cot(beta)*
   (3*m12*MW2 - m22*MW2 + m32*MW2 + 3*MA02*(pow(MB, 2)) - 
    3*Mh2*(pow(MB, 2)) + (MA02 - Mh2)*Cos(4*w)*(pow(MB, 2)) + 
    Cos(2*w)*((3*m12 - m22 + m32)*MW2 + 4*MA02*(pow(MB, 2)) - 
      4*Mh2*(pow(MB, 2))) - 4*(pow(MW, 4)))*(pow(Sec(w), 4)))/
  (Pi2*v2) - (0.0625*B0i(bb0, m32, MB2, MHH2)*Cos(alp)*(pow(MB, 2))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2)))/(Pi2*v2) - 
 (0.0625*m12*C0i(cc1, m12, m32, m22, MB2, MB2, MHH2)*Cos(alp)*
   (pow(MB, 2))*(pow(Csc(beta), 3))*(pow(Sin(alp), 2)))/
  (Pi2*v2) - (0.03125*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MB2, MB2, 
    MHH2)*Cos(alp)*(pow(MB, 2))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2)))/(Pi2*v2) - 
 (0.125*C0i(cc0, m12, m32, m22, MB2, MB2, MHH2)*Cos(alp)*(pow(MB, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2)))/(Pi2*v2) + 
 (0.03125*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, MT2, MW2, MW2)*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.125*C0i(cc0, m22, m12, m32, MT2, MW2, MW2)*(Mh2*MT2 - m22*MW2 + 
    (pow(MW, 4)))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2)))/(Pi2*v2) + (0.015625*(m12 + 5*m22 - m32)*MW2*
   C0i(cc1, m22, m12, m32, MB2, MZ2, MZ2)*(pow(Sec(w), 2))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.0008680555555555555*C0i(cc0, m22, m12, m32, MB2, MZ2, MZ2)*
   (-3*(12*m22*MW2 - 9*Mh2*(pow(MB, 2)) + 8*(pow(MW, 4))) + 
    4*Cos(2*w)*(-9*m22*MW2 + 9*Mh2*(pow(MB, 2)) + 16*(pow(MW, 4))) + 
    Cos(4*w)*(9*Mh2*(pow(MB, 2)) + 32*(pow(MW, 4))))*
   (pow(Sec(w), 4))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2)))/(Pi2*v2) + (0.027777777777777776*B0i(bb0, m32, 0, MB2)*
   (-1 + Cos(alp)*Csc(beta))*(48*Alfas*Pi + 4*MW2*(pow(v, -2))*
     (pow(Sin(w), 2))))/Pi2 - (0.013888888888888888*(m12 + m22 - m32)*
   C0i(cc1, m22, m12, m32, 0, MB2, MB2)*(-1 + Cos(alp)*Csc(beta))*
   (48*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2))))/Pi2 + 
 (0.013888888888888888*(m12 - m22 + m32)*C0i(cc2, m22, m12, m32, 0, MB2, MB2)*
   (-1 + Cos(alp)*Csc(beta))*(48*Alfas*Pi + 4*MW2*(pow(v, -2))*
     (pow(Sin(w), 2))))/Pi2 + 
 (0.013888888888888888*C0i(cc0, m22, m12, m32, 0, MB2, MB2)*
   (-1 + Cos(alp)*Csc(beta))*(-m12 - m22 + m32 + 4*(pow(MB, 2)))*
   (48*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2))))/Pi2 - 
 (0.00390625*Cos(alp)*(3*(Mh2 - MHH2)*Cos(2*alp) + 
    (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 
    4*(MA02 + (-2*M2 - MA02 + Mh2 + MHH2)*Cos(2*beta)))*
   (pow(Csc(beta), 3))*Re(A0i(aa0, MA02)))/(MA02*Pi2*v2) + 
 (0.005208333333333333*(144*Csc(beta)*(pow(MB, 2))*(pow(Cos(alp), 3))*
     (pow(Cot(beta), 2)) + 
    MA02*(61 + 4*Cos(alp)*(-12 + 2*Cos(2*w) + Cos(4*w))*Csc(beta))*
     (pow(Cot(w), 2)) + (MA02*(-93 + 61*Cos(2*w) - 20*Cos(4*w) + 
       2*Cos(6*w))*(pow(Csc(w), 2)))/2 + 8*Cos(alp)*Csc(beta)*
     (MA02*(9 - Cos(2*w) + Cos(4*w)) + 18*(pow(MB, 2))*
       (pow(Cot(beta), 2))*(pow(Sin(alp), 2))))*Re(A0i(aa0, MB2)))/
  (MA02*Pi2*v2) + (0.010416666666666666*
   (2*Cos(alp)*Csc(beta)*(64*MA02 - 24*MA02*(pow(Csc(w), 2)) + 
      (36*(pow(MC, 2))*(pow(Cot(beta), 2)) + 
        MA02*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2)))*
       (pow(Csc(w), 4))) + MA02*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 
      2*Cos(6*w))*(pow(Csc(w), 6)))*Re(A0i(aa0, MC2)))/
  (MA02*Pi2*v2*(pow(Csc(w), 4))) + 
 (0.005208333333333333*(144*Csc(beta)*(pow(MD, 2))*(pow(Cos(alp), 3))*
     (pow(Cot(beta), 2)) + 
    MA02*(61 + 4*Cos(alp)*(-12 + 2*Cos(2*w) + Cos(4*w))*Csc(beta))*
     (pow(Cot(w), 2)) + (MA02*(-93 + 61*Cos(2*w) - 20*Cos(4*w) + 
       2*Cos(6*w))*(pow(Csc(w), 2)))/2 + 8*Cos(alp)*Csc(beta)*
     (MA02*(9 - Cos(2*w) + Cos(4*w)) + 18*(pow(MD, 2))*
       (pow(Cot(beta), 2))*(pow(Sin(alp), 2))))*Re(A0i(aa0, MD2)))/
  (MA02*Pi2*v2) - (0.015625*(16*MA02 - MA02*Cos(6*w) - 
    16*MA02*Cos(alp)*Csc(beta) + MA02*Cos(alp)*Cos(6*w)*Csc(beta) - 
    10*MA02*Cos(4*w)*(-1 + Cos(alp)*Csc(beta)) + 
    Cos(2*w)*(-29*MA02 + Cos(alp)*Csc(beta)*(29*MA02 - 8*(pow(ME, 2)))) + 
    8*Cos(alp)*Csc(beta)*(pow(ME, 2)))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ME2)))/(MA02*Pi2*v2) - 
 (0.00390625*Cos(alp)*(4*MA02 + 2*(5*M2 - 6*MHH2)*Cos(2*alp) + 
    2*M2*Cos(2*(alp - 2*beta)) - Mh2*Cos(4*alp - 2*beta) + 
    MHH2*Cos(4*alp - 2*beta) - 12*M2*Cos(2*beta) - 4*MA02*Cos(2*beta) + 
    Mh2*Cos(2*beta) + 11*MHH2*Cos(2*beta))*(pow(Csc(beta), 3))*
   Re(A0i(aa0, MHH2)))/(MA02*Pi2*v2) + 
 (0.0078125*Cos(alp)*(-2*MA02 - 3*(Mh2 - MHH2)*Cos(2*alp) + 
    (-Mh2 + MHH2)*Cos(2*(alp - 2*beta)) + 8*M2*Cos(2*beta) + 
    2*MA02*Cos(2*beta) - 4*Mh2*Cos(2*beta) - 4*MHH2*Cos(2*beta) + 
    MA02*Cos(2*(beta - 2*w)) - 8*MA02*Cos(2*(beta - w)) + 16*MA02*Cos(2*w) - 
    2*MA02*Cos(4*w) - 8*MA02*Cos(2*(beta + w)) + MA02*Cos(2*(beta + 2*w)))*
   (pow(Csc(beta), 3))*Re(A0i(aa0, MHp2)))/(MA02*Pi2*v2) - 
 (0.015625*(16*MA02 - MA02*Cos(6*w) - 16*MA02*Cos(alp)*Csc(beta) + 
    MA02*Cos(alp)*Cos(6*w)*Csc(beta) - 10*MA02*Cos(4*w)*
     (-1 + Cos(alp)*Csc(beta)) + Cos(2*w)*(-29*MA02 + 
      Cos(alp)*Csc(beta)*(29*MA02 - 8*(pow(ML, 2)))) + 
    8*Cos(alp)*Csc(beta)*(pow(ML, 2)))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ML2)))/(MA02*Pi2*v2) - 
 (0.015625*(16*MA02 - MA02*Cos(6*w) - 16*MA02*Cos(alp)*Csc(beta) + 
    MA02*Cos(alp)*Cos(6*w)*Csc(beta) - 10*MA02*Cos(4*w)*
     (-1 + Cos(alp)*Csc(beta)) + Cos(2*w)*(-29*MA02 + 
      Cos(alp)*Csc(beta)*(29*MA02 - 8*(pow(MM, 2)))) + 
    8*Cos(alp)*Csc(beta)*(pow(MM, 2)))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MM2)))/(MA02*Pi2*v2) + 
 (0.005208333333333333*(144*Csc(beta)*(pow(MS, 2))*(pow(Cos(alp), 3))*
     (pow(Cot(beta), 2)) + 
    MA02*(61 + 4*Cos(alp)*(-12 + 2*Cos(2*w) + Cos(4*w))*Csc(beta))*
     (pow(Cot(w), 2)) + (MA02*(-93 + 61*Cos(2*w) - 20*Cos(4*w) + 
       2*Cos(6*w))*(pow(Csc(w), 2)))/2 + 8*Cos(alp)*Csc(beta)*
     (MA02*(9 - Cos(2*w) + Cos(4*w)) + 18*(pow(MS, 2))*
       (pow(Cot(beta), 2))*(pow(Sin(alp), 2))))*Re(A0i(aa0, MS2)))/
  (MA02*Pi2*v2) + (0.010416666666666666*
   (2*Cos(alp)*Csc(beta)*(64*MA02 - 24*MA02*(pow(Csc(w), 2)) + 
      (36*MT2*(pow(Cot(beta), 2)) + MA02*(9 - 4*Cos(2*w) + 4*Cos(4*w))*
         (pow(Cot(w), 2)))*(pow(Csc(w), 4))) + 
    MA02*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*(pow(Csc(w), 6)))*
   Re(A0i(aa0, MT2)))/(MA02*Pi2*v2*(pow(Csc(w), 4))) + 
 (0.010416666666666666*(2*Cos(alp)*Csc(beta)*
     (64*MA02 - 24*MA02*(pow(Csc(w), 2)) + 
      (36*(pow(MU, 2))*(pow(Cot(beta), 2)) + 
        MA02*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2)))*
       (pow(Csc(w), 4))) + MA02*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 
      2*Cos(6*w))*(pow(Csc(w), 6)))*Re(A0i(aa0, MU2)))/
  (MA02*Pi2*v2*(pow(Csc(w), 4))) - 
 (1.125*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb0, 0, MW2, MW2)))/(Pi2*v2) - 
 (0.027777777777777776*(-1 + Cos(alp)*Csc(beta))*
   (48*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*
   Re(B0i(bb0, MB2, 0, MB2)))/Pi2 - 
 (0.0625*Cos(alp)*Csc(beta)*(pow(MB, 2))*(pow(Cot(beta), 2))*
   Re(B0i(bb0, MB2, MA02, MB2)))/(Pi2*v2) + 
 (0.0625*(pow(MB, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(bb0, MB2, MB2, Mh2)))/(Pi2*v2) + 
 (0.0625*Cos(alp)*(pow(MB, 2))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, MB2, MB2, MHH2)))/(Pi2*v2) - 
 (0.003472222222222222*(-1 + Cos(alp)*Csc(beta))*(-12*MW2 + 4*MW2*Cos(4*w) + 
    9*(pow(MB, 2)) + Cos(2*w)*(8*MW2 + 9*(pow(MB, 2))))*
   (pow(Sec(w), 2))*Re(B0i(bb0, MB2, MB2, MZ2)))/(Pi2*v2) - 
 (0.125*MT2*Cos(alp)*Csc(beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb0, MB2, MHp2, MT2)))/(Pi2*v2) - 
 (0.125*MT2*(-1 + Cos(alp)*Csc(beta))*Re(B0i(bb0, MB2, MT2, MW2)))/(Pi2*v2) + 
 (0.75*Cos(alp)*(pow(MB, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, MB2, MB2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.75*Cos(alp)*(pow(MC, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, Mh2, MC2, MC2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*(pow(MD, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, MD2, MD2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.25*Cos(alp)*Csc(beta)*(pow(ME, 4))*
   (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, Mh2, ME2, ME2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.25*Cos(alp)*Csc(beta)*(pow(ML, 4))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, ML2, ML2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.25*Cos(alp)*Csc(beta)*(pow(MM, 4))*
   (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, Mh2, MM2, MM2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*(pow(MS, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, MS2, MS2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.75*Cos(alp)*(pow(MT, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, Mh2, MT2, MT2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*(pow(MU, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, MU2, MU2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.75*Cos(alp)*(pow(MB, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, MB2, MB2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*(pow(MC, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, MHH2, MC2, MC2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.75*Cos(alp)*(pow(MD, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, MD2, MD2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*Cos(alp)*Csc(beta)*(pow(ME, 4))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, MHH2, ME2, ME2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.25*Cos(alp)*Csc(beta)*(pow(ML, 4))*
   (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, ML2, ML2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*Cos(alp)*Csc(beta)*(pow(MM, 4))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, MHH2, MM2, MM2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.75*Cos(alp)*(pow(MS, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, MS2, MS2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*(pow(MT, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, MHH2, MT2, MT2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.75*Cos(alp)*(pow(MU, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, MU2, MU2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.5*MW2*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*Re(B0i(bb0, MW2, 0, MW2)))/
  (Pi2*v2) + (0.375*(-MT2 + MW2)*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*
   (pow(Csc(w), 2))*Re(B0i(bb0, MW2, MB2, MT2)))/(Pi2*v2) - 
 (0.375*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*
   (-1 + (pow(Cot(w), 2)))*Re(B0i(bb0, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*(MW2 - (pow(MU, 2)))*
   (pow(Csc(w), 2))*Re(B0i(bb0, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.125*MW2*Cos(2*w)*(pow(Csc(w), 2))*
   (-1 + Cos(alp)*Csc(beta)*(pow(Sin(alp - beta), 2)))*
   Re(B0i(bb0, MW2, Mh2, MW2)))/(Pi2*v2) + 
 (0.125*MW2*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb0, MW2, MHH2, MW2)))/(Pi2*v2) - 
 (0.015625*MW2*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*(pow(Sec(w), 2))*
   Re(B0i(bb0, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MD, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*(-1 + Cos(alp)*Csc(beta))*(pow(ME, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, ME2, ME2)))/(Pi2*v2) - 
 (0.125*MW2*(pow(Csc(w), 2))*(-1 + Cos(alp)*Csc(beta)*
     (pow(Sin(alp - beta), 2)))*Re(B0i(bb0, MZ2, Mh2, MZ2)))/(Pi2*v2) - 
 (0.125*MW2*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb0, MZ2, MHH2, MZ2)))/(Pi2*v2) + 
 (0.0625*(-1 + Cos(alp)*Csc(beta))*(pow(ML, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*(-1 + Cos(alp)*Csc(beta))*(pow(MM, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MS, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.1875*MT2*(-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MU, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.0625*MW2*(5 + 9*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MB2, MB2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MD2, MD2)))/(Pi2*v2) + 
 (0.5*(-1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, ME2, ME2)))/(Pi2*v2) - 
 (0.5*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MHp2, MHp2)))/(Pi2*v2) + 
 (0.5*(-1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, ML2, ML2)))/(Pi2*v2) + 
 (0.5*(-1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MM2, MM2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MT2, MT2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.5*(2 + 3*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.25*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, 0, ME2)))/(Pi2*v2) + 
 (0.25*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, 0, ML2)))/(Pi2*v2) + 
 (0.25*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, 0, MM2)))/(Pi2*v2) - 
 (1.*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*Re(B0i(bb00, MW2, 0, MW2)))/
  (Pi2*v2) - (0.125*Cos(alp)*Csc(beta)*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, MA02, MHp2)))/(Pi2*v2) + 
 (0.75*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, MB2, MT2)))/(Pi2*v2) + 
 (0.75*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.75*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, MD2, MU2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MW2, Mh2, MHp2)))/(Pi2*v2) - 
 (0.125*Cos(2*w)*(pow(Csc(w), 2))*
   (-1 + Cos(alp)*Csc(beta)*(pow(Sin(alp - beta), 2)))*
   Re(B0i(bb00, MW2, Mh2, MW2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MW2, MHH2, MHp2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MW2, MHH2, MW2)))/(Pi2*v2) - 
 (0.125*(2 + 5*Cos(2*w) + 2*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MW2, MW2, MZ2)))/(Pi2*v2) - 
 (0.375*(-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, 0, 0)))/(Pi2*v2) + 
 (0.125*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MA02, Mh2)))/(Pi2*v2) + 
 (0.125*Cos(alp)*Csc(beta)*(pow(Cot(w), 2))*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb00, MZ2, MA02, MHH2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Cot(w), 2))*Re(B0i(bb00, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.041666666666666664*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MC2, MC2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Cot(w), 2))*Re(B0i(bb00, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*(pow(Cot(w), 2))*(-1 + Cos(alp)*Csc(beta)*
     (pow(Sin(alp - beta), 2)))*Re(B0i(bb00, MZ2, Mh2, MZ2)))/(Pi2*v2) + 
 (0.125*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MHH2, MZ2)))/(Pi2*v2) + 
 (0.125*Cos(alp)*Csc(beta)*(pow(Cos(2*w), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MHp2, MHp2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Cot(w), 2))*Re(B0i(bb00, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.041666666666666664*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.041666666666666664*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.0625*(7 + 8*Cos(2*w) + 3*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Cot(w), 2))*Re(B0i(bb00, MZ2, MW2, MW2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MB2, MB2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MC2, MC2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MD2, MD2)))/(Pi2*v2) - 
 (0.5*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, ME2, ME2)))/(Pi2*v2) - 
 (0.5*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, ML2, ML2)))/(Pi2*v2) - 
 (0.5*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MM2, MM2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MS2, MS2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MT2, MT2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.25*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.03125*MW2*Cos(alp - beta)*(MHH2*Cos(2*alp - beta) + 
    (-2*Mh2 + MHH2)*Cos(beta))*Csc(beta)*(pow(Sec(w), 2))*
   Re(B0i(bb1, Mh2, MA02, MZ2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MB, 2))*(-2*Mh2 + 2*MHH2 - 
    Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*(pow(Csc(beta), 3)))*
   Re(B0i(bb1, Mh2, MB2, MB2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MC, 2))*(-2*Mh2 + 2*MHH2 - 
    Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*(pow(Csc(beta), 3)))*
   Re(B0i(bb1, Mh2, MC2, MC2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MD, 2))*(-2*Mh2 + 2*MHH2 - 
    Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*(pow(Csc(beta), 3)))*
   Re(B0i(bb1, Mh2, MD2, MD2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.125*(pow(ME, 2))*(Mh2 - MHH2 + MHH2*Cos(alp)*Csc(beta)*
     (pow(Sec(beta), 2))*(pow(Sin(alp), 2)))*
   Re(B0i(bb1, Mh2, ME2, ME2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0625*MW2*Cos(alp - beta)*(MHH2*Cos(2*alp - beta) + 
    (-2*Mh2 + MHH2)*Cos(beta))*Csc(beta)*Re(B0i(bb1, Mh2, MHp2, MW2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.125*(pow(ML, 2))*
   (Mh2 - MHH2 + MHH2*Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(bb1, Mh2, ML2, ML2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.125*(pow(MM, 2))*
   (Mh2 - MHH2 + MHH2*Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(bb1, Mh2, MM2, MM2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.1875*(pow(MS, 2))*
   (-2*Mh2 + 2*MHH2 - Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*
     (pow(Csc(beta), 3)))*Re(B0i(bb1, Mh2, MS2, MS2)))/
  ((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*MT2*(-2*Mh2 + 2*MHH2 - Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*
     (pow(Csc(beta), 3)))*Re(B0i(bb1, Mh2, MT2, MT2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.1875*(pow(MU, 2))*
   (-2*Mh2 + 2*MHH2 - Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*
     (pow(Csc(beta), 3)))*Re(B0i(bb1, Mh2, MU2, MU2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.375*MHH2*Cos(alp)*(pow(MB, 2))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, MB2, MB2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*(pow(MC, 2))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb1, MHH2, MC2, MC2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.375*MHH2*Cos(alp)*(pow(MD, 2))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, MD2, MD2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.125*MHH2*Cos(alp)*Csc(beta)*(pow(ME, 2))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb1, MHH2, ME2, ME2)))/
  ((-Mh2 + MHH2)*Pi2*v2) - (0.125*MHH2*Cos(alp)*Csc(beta)*(pow(ML, 2))*
   (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, ML2, ML2)))/((-Mh2 + MHH2)*Pi2*v2) - 
 (0.125*MHH2*Cos(alp)*Csc(beta)*(pow(MM, 2))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb1, MHH2, MM2, MM2)))/
  ((-Mh2 + MHH2)*Pi2*v2) - (0.375*MHH2*Cos(alp)*(pow(MS, 2))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, MS2, MS2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*MT2*Cos(alp)*(pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, MT2, MT2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*(pow(MU, 2))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb1, MHH2, MU2, MU2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.125*MW2*(-1 + Cos(alp)*Csc(beta))*
   (-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ME2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, 0, ML2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, 0, MM2)))/(Pi2*v2) - 
 (0.25*MW2*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*Re(B0i(bb1, MW2, 0, MW2)))/
  (Pi2*v2) + (0.375*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MB2, MT2)))/(Pi2*v2) - 
 (0.375*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.25*MW2*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb1, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.1875*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, 0, 0)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.25*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Cos(w), 2))*
   (pow(Cot(w), 2))*Re(B0i(bb1, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.027777777777777776*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*
   (48*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*
   Re(B0i(dbb0, MB2, 0, MB2)))/Pi2 - 
 (0.125*(pow(MB, 4))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb0, MB2, MB2, Mh2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*(pow(MB, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(dbb0, MB2, MB2, MHH2)))/(Pi2*v2) + 
 (0.006944444444444444*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*
   (-12*MW2 + 4*MW2*Cos(4*w) + 9*(pow(MB, 2)) + 
    Cos(2*w)*(8*MW2 + 9*(pow(MB, 2))))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MB2, MB2, MZ2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Csc(beta)*(pow(MB, 2))*(-MT2 + (pow(MB, 2)))*
   (pow(Cot(beta), 2))*Re(B0i(dbb0, MB2, MHp2, MT2)))/(Pi2*v2) + 
 (0.25*MT2*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*
   Re(B0i(dbb0, MB2, MT2, MW2)))/(Pi2*v2) - 
 (0.0009765625*Cos(alp)*(pow((2*MA02 - Mh2)*Cos(alp - 3*beta) + 
      (4*M2 - 2*MA02 - 3*Mh2)*Cos(alp + beta), 2))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*Re(B0i(dbb0, Mh2, MA02, MA02)))/(Pi2*v2) + 
 (0.75*(pow(MB, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.75*(pow(MC, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.75*(pow(MD, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.25*(pow(ME, 4))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb0, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.0087890625*(16*(pow(Mh, 4)) - 
    Cos(alp)*(pow(M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
        (2*M2 - 3*Mh2)*Cos(alp + beta), 2))*(pow(Csc(beta), 3))*
     (pow(Sec(beta), 2)))*Re(B0i(dbb0, Mh2, Mh2, Mh2)))/(Pi2*v2) - 
 (0.0078125*Cos(alp)*(pow(Cos(alp - beta), 2))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*(pow((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + 
      M2*Sin(2*beta), 2))*Re(B0i(dbb0, Mh2, Mh2, MHH2)))/(Pi2*v2) - 
 (0.00390625*Cos(alp)*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   (pow(Sin(alp - beta), 2))*
   (pow((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta), 2))*
   Re(B0i(dbb0, Mh2, MHH2, MHH2)))/(Pi2*v2) - 
 (0.001953125*Cos(alp)*(pow((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + 
      (-4*M2 + 3*Mh2 + 2*MHp2)*Cos(alp + beta), 2))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*Re(B0i(dbb0, Mh2, MHp2, MHp2)))/(Pi2*v2) - 
 (0.0625*Cos(alp)*Csc(beta)*(MHp2*(MHp2 - MW2) - Mh2*(2*MHp2 + MW2) + 
    (pow(Mh, 4)))*(pow(Cos(alp - beta), 2))*
   Re(B0i(dbb0, Mh2, MHp2, MW2)))/(Pi2*v2) + 
 (0.25*(pow(ML, 4))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb0, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.25*(pow(MM, 4))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb0, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.75*(pow(MS, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.75*(pow(MT, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.75*(pow(MU, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MU2, MU2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MB2, MB2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MD2, MD2)))/(Pi2*v2) + 
 (1.*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, ME2, ME2)))/(Pi2*v2) - 
 (0.5*MW2*Cos(alp)*Csc(beta)*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MHp2, MHp2)))/(Pi2*v2) + 
 (1.*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, ML2, ML2)))/(Pi2*v2) + 
 (1.*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MM2, MM2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MT2, MT2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (1.5*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MW2, MW2)))/(Pi2*v2) - 
 (0.027777777777777776*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*
   (48*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*
   Re(B0i(dbb1, MB2, 0, MB2)))/Pi2 - 
 (0.125*Cos(alp)*Csc(beta)*(pow(MB, 4))*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MB2, MA02, MB2)))/(Pi2*v2) + 
 (0.125*(pow(MB, 4))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, MB2, MB2, Mh2)))/(Pi2*v2) + 
 (0.125*Cos(alp)*(pow(MB, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(dbb1, MB2, MB2, MHH2)))/(Pi2*v2) + 
 (0.006944444444444444*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*
   (12*MW2 + 2*MW2*Cos(4*w) + 9*(pow(MB, 2)) + 
    Cos(2*w)*(4*MW2 + 9*(pow(MB, 2))))*(pow(Sec(w), 2))*
   Re(B0i(dbb1, MB2, MB2, MZ2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Csc(beta)*(pow(MB, 2))*(MT2 + (pow(MB, 2)))*
   (pow(Cot(beta), 2))*Re(B0i(dbb1, MB2, MHp2, MT2)))/(Pi2*v2) + 
 (0.125*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*
   (MT2 + 2*MW2 + (pow(MB, 2)))*Re(B0i(dbb1, MB2, MT2, MW2)))/(Pi2*v2) - 
 (0.0625*Mh2*MW2*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Re(B0i(dbb1, Mh2, MA02, MZ2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MB, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MC, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MD, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ME, 2))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb1, Mh2, ME2, ME2)))/(Pi2*v2) - 
 (0.125*Mh2*MW2*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   Re(B0i(dbb1, Mh2, MHp2, MW2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ML, 2))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(MM, 2))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MS, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*Mh2*MT2*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MU, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*Mh2*MW2*(-1 + Cos(alp)*Csc(beta)*(pow(Sin(alp - beta), 2)))*
   Re(B0i(dbb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*Mh2*MW2*(pow(Sec(w), 2))*
   (-1 + Cos(alp)*Csc(beta)*(pow(Sin(alp - beta), 2)))*
   Re(B0i(dbb1, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (0.046875*C0i(cc0, m22, m12, m32, MB2, Mh2, Mh2)*(pow(MB, 2))*
   (-4*Mh2 - (M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
      (2*M2 - 3*Mh2)*Cos(alp + beta))*(pow(Cos(alp), 2))*
     (pow(Csc(beta), 3))*Sec(beta)))/(Pi2*v2) + 
 (0.0625*B0i(bb0, m32, MHp2, MT2)*
   (Cos(alp)*Cot(beta)*(MW2*Cos(beta) + 2*MT2*Cot(beta)*Csc(beta)) + 
    MW2*Cos(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.03125*B0i(bb0, m32, MA02, MB2)*
   (Cos(alp)*Cot(beta)*(2*Cot(beta)*Csc(beta)*(pow(MB, 2)) + 
      MW2*Cos(beta)*(pow(Sec(w), 2))) + MW2*Cos(beta)*
     (pow(Sec(w), 2))*Sin(alp)))/(Pi2*v2) + 
 (0.00390625*Re(A0i(aa0, Mh2))*(8*MA02 + 5*M2*Cos(3*alp)*
     (pow(Csc(beta), 3)) + Cos(alp)*(5*M2 - 4*MA02 + 
      2*M2*Cos(2*(alp - 2*beta)) + (-Mh2 + MHH2)*Cos(4*alp - 2*beta) + 
      12*M2*Cos(2*beta) + 4*MA02*Cos(2*beta) - 11*Mh2*Cos(2*beta) - 
      MHH2*Cos(2*beta))*(pow(Csc(beta), 3)) - 
    3*Mh2*Csc(alp)*(pow(Csc(beta), 3))*Sin(4*alp)))/(MA02*Pi2*v2) - 
 (0.03125*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, MT2, MW2, MW2)*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.015625*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, MB2, MZ2, MZ2)*
   (pow(Sec(w), 2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.03125*(Mh2 - MHH2)*A0i(aa0, MW2)*Cos(alp)*Cot(beta)*Csc(beta)*
   Sin(2*(alp - beta)))/(MA02*Pi2*v2) - 
 (0.046875*(Mh2 - MHH2)*A0i(aa0, MZ2)*Cos(alp)*Cot(beta)*Csc(beta)*
   Sin(2*(alp - beta)))/(MA02*Pi2*v2) - 
 (0.03125*MW2*B0i(bb1, MA02, Mh2, MZ2)*Cos(alp)*Cot(beta)*Csc(beta)*
   (pow(Sec(w), 2))*Sin(2*(alp - beta)))/(Pi2*v2) + 
 (0.03125*MW2*B0i(bb1, MA02, MHH2, MZ2)*Cos(alp)*Cot(beta)*Csc(beta)*
   (pow(Sec(w), 2))*Sin(2*(alp - beta)))/(Pi2*v2) + 
 (0.0078125*Mh2*B0i(bb0, 0, Mh2, MZ2)*Cos(alp)*(MA02 - Mh2 + 2*MW2 + 
    (MA02 - Mh2)*Cos(2*w))*Cot(beta)*Csc(beta)*(pow(Sec(w), 2))*
   Sin(2*(alp - beta)))/(MA02*Pi2*v2) - 
 (0.0078125*B0i(bb0, MA02, MHH2, MZ2)*Cos(alp)*(-(MHH2*(MHH2 - 2*MW2)) + 
    MA02*(MHH2 + 2*MW2) + (MA02 - MHH2)*MHH2*Cos(2*w))*Cot(beta)*Csc(beta)*
   (pow(Sec(w), 2))*Sin(2*(alp - beta)))/(MA02*Pi2*v2) + 
 (0.0078125*MHH2*B0i(bb0, 0, MHH2, MZ2)*Cos(alp)*
   (-MA02 + MHH2 - 2*MW2 + (-MA02 + MHH2)*Cos(2*w))*Cot(beta)*Csc(beta)*
   (pow(Sec(w), 2))*Sin(2*(alp - beta)))/(MA02*Pi2*v2) + 
 (0.0078125*B0i(bb0, MA02, Mh2, MZ2)*Cos(alp)*Cot(beta)*Csc(beta)*
   (2*Mh2*MW2 + MA02*(Mh2 + 2*MW2) + (MA02 - Mh2)*Mh2*Cos(2*w) - 
    (pow(Mh, 4)))*(pow(Sec(w), 2))*Sin(2*(alp - beta)))/
  (MA02*Pi2*v2) + (0.015625*Csc(beta)*((MA02 - Mh2)*(MA02 - MHH2) - 
    (MA02 + MHH2)*MW2*(pow(Sec(w), 2)))*Re(B0i(bb0, MHH2, MA02, MZ2))*
   Sin(alp)*Sin(2*(alp - beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.03125*(Mh2*(MHH2 - MHp2) + MHp2*(MHp2 - MW2) - MHH2*(MHp2 + MW2))*
   Csc(beta)*Re(B0i(bb0, MHH2, MHp2, MW2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - 
 (0.015625*Csc(beta)*(Mh2*MHH2 - 2*MHH2*MW2 + 12*(pow(MW, 4)))*
   Re(B0i(bb0, MHH2, MW2, MW2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.0009765625*Csc(beta)*
   (3*Mh2*MHH2 - 8*MHH2*MW2 + 4*MHH2*(Mh2 - 2*MW2)*Cos(2*w) + 
    Mh2*MHH2*Cos(4*w) + 96*(pow(MW, 4)))*(pow(Sec(w), 4))*
   Re(B0i(bb0, MHH2, MZ2, MZ2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.03125*MHH2*MW2*Csc(beta)*(pow(Sec(w), 2))*
   Re(B0i(bb1, MHH2, MA02, MZ2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.0625*MHH2*MW2*Csc(beta)*
   Re(B0i(bb1, MHH2, MHp2, MW2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.0625*MHH2*MW2*Csc(beta)*
   Re(B0i(bb1, MHH2, MW2, MW2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.03125*MHH2*MW2*Csc(beta)*(pow(Sec(w), 2))*
   Re(B0i(bb1, MHH2, MZ2, MZ2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.03125*Csc(beta)*Re(B0i(bb0, Mh2, MHp2, MW2))*
   (2*(Mh2 - MHH2)*MW2*Cos(alp)*(pow(Cos(alp - beta), 2)) + 
    (MHp2*(MHH2 - MHp2 + MW2) + Mh2*(-MHH2 + MHp2 + MW2))*Sin(alp)*
     Sin(2*(alp - beta))))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.015625*MW2*Csc(beta)*Re(B0i(bb0, Mh2, MA02, MZ2))*
   (-2*Cos(alp)*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 2)) + 
    (pow(Mh2 - MHH2, -1))*(pow(MW, -2))*((MA02 - Mh2)*(MA02 - MHH2) - 
      (MA02 + Mh2)*MW2*(pow(Sec(w), 2)))*Sin(alp)*Sin(2*(alp - beta))))/
  (Pi2*v2) + (0.015625*MW2*Re(B0i(bb0, Mh2, MZ2, MZ2))*
   (-8*(pow(Csc(2*w), 2)) - Csc(beta)*(pow(Csc(w), 2))*
     (pow(Sec(w), 2))*(-2*Cos(alp)*(pow(Sin(alp - beta), 2)) - 
      ((pow(Mh2 - MHH2, -1))*(pow(MW, -2))*(3*Mh2*MHH2 - 8*Mh2*MW2 + 
         4*Mh2*(MHH2 - 2*MW2)*Cos(2*w) + Mh2*MHH2*Cos(4*w) + 
         96*(pow(MW, 4)))*(pow(Sec(w), 2))*Sin(alp)*
        Sin(2*(alp - beta)))/16)))/(Pi2*v2*(pow(Csc(w), 2))) + 
 (0.0625*MW2*Re(B0i(bb1, Mh2, MW2, MW2))*(2*(Mh2 - MHH2) + 
    Sin(alp - beta)*(2*Mh2 - MHH2 + MHH2*Csc(beta)*Sin(2*alp - beta))))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.03125*MW2*(pow(Sec(w), 2))*
   Re(B0i(bb1, Mh2, MZ2, MZ2))*(2*(Mh2 - MHH2) + 
    Sin(alp - beta)*(2*Mh2 - MHH2 + MHH2*Csc(beta)*Sin(2*alp - beta))))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.015625*Re(B0i(bb0, Mh2, MW2, MW2))*
   (-4*MW2 + (pow(Mh2 - MHH2, -1))*Sin(alp - beta)*
     (Mh2*(MHH2 - 4*MW2) + 2*MW2*(MHH2 + 6*MW2) + 
      Csc(beta)*(Mh2*MHH2 - 2*MHH2*MW2 + 12*(pow(MW, 4)))*
       Sin(2*alp - beta))))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ME2, ME2))*(Cos(alp) - Sin(beta)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ML2, ML2))*(Cos(alp) - Sin(beta)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MM2, MM2))*(Cos(alp) - Sin(beta)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, MC2, MC2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, ME2, ME2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, ML2, ML2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, MM2, MM2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*
   Csc(beta)*(pow(Csc(w), 2))*Re(B0i(bb1, MZ2, MT2, MT2))*
   (Cos(alp) - Sin(beta)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, MU2, MU2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.03125*Cos(alp)*(-(pow(MA02 - Mh2, 2)) + 
    (MA02 + Mh2)*MW2*(pow(Sec(w), 2)))*
   (pow(Cos(alp)*Cot(beta) + Sin(alp), 2))*Re(B0i(dbb0, Mh2, MA02, MZ2))*
   Sin(beta))/(Pi2*v2) + 
 (0.0078125*Csc(beta)*(-2*Mh2*MW2 + (pow(Mh, 4)) + 12*(pow(MW, 4)))*
   Re(B0i(dbb0, Mh2, MW2, MW2))*(-2*Cos(alp) + Cos(alp - 2*beta) + 
    Cos(3*alp - 2*beta) + 4*Sin(beta)))/(Pi2*v2) + 
 (0.00048828125*Csc(beta)*(-8*Mh2*MW2 + 4*Mh2*(Mh2 - 2*MW2)*Cos(2*w) + 
    3*(pow(Mh, 4)) + Cos(4*w)*(pow(Mh, 4)) + 96*(pow(MW, 4)))*
   (pow(Sec(w), 4))*Re(B0i(dbb0, Mh2, MZ2, MZ2))*
   (-2*Cos(alp) + Cos(alp - 2*beta) + Cos(3*alp - 2*beta) + 4*Sin(beta)))/
  (Pi2*v2) + (0.0078125*Re(A0i(aa0, MZ2))*
   (2*MA02*(-13 + 24*Cos(alp)*Csc(beta)*(pow(Cos(w), 2)))*
     (pow(Cot(w), 2)) - 2*MA02*(-13 + (5 + 6*Cos(4*w))*
       (pow(Csc(w), 2))) + Cos(alp)*(pow(Csc(beta), 2))*
     ((Mh2 - MHH2)*Sin(2*alp - 3*beta) + (Mh2 - MHH2)*Sin(2*alp - beta) - 
      4*MA02*(7 + 6*Cos(2*w))*Sin(beta))))/(MA02*Pi2*v2) + 
 (0.03125*C0i(cc0, m22, m12, m32, MB2, MHH2, MHH2)*(pow(MB, 2))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*Sec(beta)*Sin(alp - beta)*
   ((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) + 
 (0.005859375*Cos(alp - beta)*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh2)*Cos(3*alp - beta) + (2*M2 - 3*Mh2)*Cos(alp + beta))*
   (pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, Mh2, Mh2, Mh2))*Sin(alp)*((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.005859375*Cos(alp - beta)*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh2)*Cos(3*alp - beta) + (2*M2 - 3*Mh2)*Cos(alp + beta))*
   (pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, Mh2, Mh2))*Sin(alp)*((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.015625*C0i(cc0, m22, m12, m32, MB2, Mh2, MHH2)*Cos(alp - beta)*
   (pow(MB, 2))*(pow(Csc(beta), 3))*Sec(beta)*Sin(2*alp)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) - 
 (0.015625*C0i(cc0, m32, m12, m22, MB2, Mh2, MHH2)*Cos(alp - beta)*
   (pow(MB, 2))*(pow(Csc(beta), 3))*Sec(beta)*Sin(2*alp)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) - 
 (0.001953125*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, Mh2, Mh2, MHH2))*Sin(alp)*Sin(2*(alp - beta))*
   (-9*M2*Mh2 - 9*M2*MHH2 + 5*Mh2*MHH2 + 8*(pow(M2, 2)) + 
    2*(pow(Mh, 4)) + 2*(pow(MHH2, 2)) - 
    Cos(4*alp)*(5*Mh2*MHH2 - 9*M2*(Mh2 + MHH2) + 9*(pow(M2, 2)) + 
      2*(pow(Mh, 4)) + 2*(pow(MHH2, 2))) + 
    (pow(M2, 2))*(pow(Cos(beta), 4)) - 6*(pow(M2, 2))*
     (pow(Cos(beta), 2))*(pow(Sin(beta), 2)) + 
    (pow(M2, 2))*(pow(Sin(beta), 4)) - 2*M2*(Mh2 - MHH2)*Sin(2*alp)*
     Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.001953125*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, Mh2, MHH2))*Sin(alp)*Sin(2*(alp - beta))*
   (-9*M2*Mh2 - 9*M2*MHH2 + 5*Mh2*MHH2 + 8*(pow(M2, 2)) + 
    2*(pow(Mh, 4)) + 2*(pow(MHH2, 2)) - 
    Cos(4*alp)*(5*Mh2*MHH2 - 9*M2*(Mh2 + MHH2) + 9*(pow(M2, 2)) + 
      2*(pow(Mh, 4)) + 2*(pow(MHH2, 2))) + 
    (pow(M2, 2))*(pow(Cos(beta), 4)) - 6*(pow(M2, 2))*
     (pow(Cos(beta), 2))*(pow(Sin(beta), 2)) + 
    (pow(M2, 2))*(pow(Sin(beta), 4)) - 2*M2*(Mh2 - MHH2)*Sin(2*alp)*
     Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0234375*Csc(beta)*(pow(Csc(2*beta), 2))*Re(B0i(bb0, Mh2, MHH2, MHH2))*
   Sin(alp)*Sin(alp - beta)*((-3*M2 + Mh2 + 2*MHH2)*Sin(2*alp) - 
    M2*Sin(2*beta))*(M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0234375*Csc(beta)*(pow(Csc(2*beta), 2))*
   Re(B0i(bb0, MHH2, MHH2, MHH2))*Sin(alp)*Sin(alp - beta)*
   ((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta))*
   (M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0009765625*((2*MA02 - Mh2)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh2)*Cos(alp + beta))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh2, MA02, MA02))*Sin(alp)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0009765625*((-2*MA02 + Mh2)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh2)*Cos(alp + beta))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MA02, MA02))*Sin(alp)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0078125*(MA02 - MHH2)*B0i(bb0, 0, MA02, MHH2)*Cos(alp)*
   (pow(Csc(beta), 3))*Sin(alp - beta)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(MA02*Pi2*v2) - 
 (0.0078125*(MA02 - MHH2)*B0i(bb0, MA02, MA02, MHH2)*Cos(alp)*
   (pow(Csc(beta), 3))*Sin(alp - beta)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(MA02*Pi2*v2) + 
 (0.0078125*A0i(aa0, MHH2)*Cos(alp)*(pow(Csc(beta), 3))*Sin(alp - beta)*
   (-2*(M2 - MA02)*Sin(alp - 3*beta) + (Mh2 - MHH2)*Sin(3*alp - beta) + 
    (-2*M2 - 2*MA02 + Mh2 + 3*MHH2)*Sin(alp + beta)))/(MA02*Pi2*v2) - 
 (0.001953125*((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, Mh2, MHp2, MHp2))*Sin(alp)*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.001953125*((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, MHp2, MHp2))*Sin(alp)*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.00390625*(pow(Csc(w), 2))*Re(A0i(aa0, MW2))*
   (44 - 2*Cos(2*w) + 60*Cos(4*w) - 6*Cos(6*w) + Cos(alp)*(pow(MA02, -1))*
     (pow(Csc(beta), 2))*(2*(Mh2 - MHH2)*Sin(2*alp - 3*beta) + 
      2*(Mh2 - MHH2)*Sin(2*alp - beta) - 44*MA02*Sin(beta) + 
      3*MA02*Sin(beta - 6*w) - 30*MA02*Sin(beta - 4*w) - 
      Mh2*Sin(2*alp - 3*beta - 2*w) + MHH2*Sin(2*alp - 3*beta - 2*w) - 
      Mh2*Sin(2*alp - beta - 2*w) + MHH2*Sin(2*alp - beta - 2*w) + 
      MA02*Sin(beta - 2*w) - Mh2*Sin(2*alp - 3*beta + 2*w) + 
      MHH2*Sin(2*alp - 3*beta + 2*w) - Mh2*Sin(2*alp - beta + 2*w) + 
      MHH2*Sin(2*alp - beta + 2*w) + MA02*Sin(beta + 2*w) - 
      30*MA02*Sin(beta + 4*w) + 3*MA02*Sin(beta + 6*w))))/(Pi2*v2)
;
}

ComplexType dKappahcc(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32)
{
return 0. - (0.03125*(m12 - m22 - m32)*MW2*C0i(cc1, m22, m32, m12, MHp2, MS2, MW2)*
   Cos(alp - beta)*Cot(beta))/(Pi2*v2) - 
 (0.03125*(3*m12 - 3*m22 + m32)*MW2*C0i(cc1, m32, m22, m12, MHp2, MS2, MW2)*
   Cos(alp - beta)*Cot(beta))/(Pi2*v2) - 
 (0.03125*(m12 - m22 + m32)*MW2*C0i(cc2, m22, m32, m12, MHp2, MS2, MW2)*
   Cos(alp - beta)*Cot(beta))/(Pi2*v2) - 
 (0.03125*(5*m12 + m22 - m32)*MW2*C0i(cc2, m32, m22, m12, MHp2, MS2, MW2)*
   Cos(alp - beta)*Cot(beta))/(Pi2*v2) + 1.*(-1 + Cos(alp)*Csc(beta)) + 
 (0.25*A0i(aa0, ME2)*Cos(alp)*Csc(beta)*(pow(ME, 2)))/(MA02*Pi2*v2) + 
 (0.125*B0i(bb1, MA02, ME2, ME2)*Cos(alp)*Csc(beta)*(pow(ME, 2)))/
  (Pi2*v2) + (0.25*A0i(aa0, ML2)*Cos(alp)*Csc(beta)*(pow(ML, 2)))/
  (MA02*Pi2*v2) + (0.125*B0i(bb1, MA02, ML2, ML2)*Cos(alp)*Csc(beta)*
   (pow(ML, 2)))/(Pi2*v2) + (0.25*A0i(aa0, MM2)*Cos(alp)*Csc(beta)*
   (pow(MM, 2)))/(MA02*Pi2*v2) + 
 (0.125*B0i(bb1, MA02, MM2, MM2)*Cos(alp)*Csc(beta)*(pow(MM, 2)))/
  (Pi2*v2) + (0.125*m12*C0i(cc1, m12, m32, m22, MS2, MS2, MW2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MS, 2)))/(Pi2*v2) + 
 (0.0625*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MS2, MS2, MW2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MS, 2)))/(Pi2*v2) + 
 (0.03125*C0i(cc0, m22, m32, m12, MHp2, MS2, MW2)*Cos(alp - beta)*Cot(beta)*
   (-((m12 + m22 - m32)*MW2) + MHp2*(2*MW2 - 4*(pow(MS, 2))) + 
    4*Mh2*(pow(MS, 2))))/(Pi2*v2) + 
 (0.03125*C0i(cc0, m32, m22, m12, MHp2, MS2, MW2)*Cos(alp - beta)*Cot(beta)*
   (MW2*(-3*m12 + m22 - m32 + 2*MW2) + 4*Mh2*(pow(MS, 2)) - 
    4*MHp2*(pow(MS, 2))))/(Pi2*v2) + 
 (0.03125*B0i(bb0, m32, MS2, MW2)*(-(MW2*Cos(alp - 2*beta)*Csc(beta)) - 
    4*(MW2 + (pow(MS, 2))) + Cos(alp)*Csc(beta)*
     (3*MW2 + 4*(pow(MS, 2)))))/(Pi2*v2) + 
 (0.25*C0i(cc0, m12, m32, m22, MS2, MS2, MW2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MS, 4)))/(Pi2*v2) - (0.75*MT2*A0i(aa0, MT2)*Cos(alp)*Csc(beta)*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*MT2*B0i(bb1, MA02, MT2, MT2)*Cos(alp)*Csc(beta)*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MB2)*Cos(alp)*Csc(beta)*(pow(MB, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MB2, MB2)*Cos(alp)*Csc(beta)*(pow(MB, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MC2)*Cos(alp)*Csc(beta)*(pow(MC, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MC2, MC2)*Cos(alp)*Csc(beta)*(pow(MC, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.03125*(m12 + m22 - m32)*C0i(cc1, m22, m12, m32, MA02, MC2, MC2)*Cos(alp)*
   Csc(beta)*(pow(MC, 2))*(pow(Cot(beta), 2)))/(Pi2*v2) + 
 (0.03125*(m12 - m22 + m32)*C0i(cc2, m22, m12, m32, MA02, MC2, MC2)*Cos(alp)*
   Csc(beta)*(pow(MC, 2))*(pow(Cot(beta), 2)))/(Pi2*v2) + 
 (0.03125*C0i(cc0, m22, m12, m32, MA02, MC2, MC2)*Cos(alp)*Csc(beta)*
   (pow(MC, 2))*(-m12 - m22 + m32 + 4*(pow(MC, 2)))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MD2)*Cos(alp)*Csc(beta)*(pow(MD, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MD2, MD2)*Cos(alp)*Csc(beta)*(pow(MD, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MS2)*Cos(alp)*Csc(beta)*(pow(MS, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MS2, MS2)*Cos(alp)*Csc(beta)*(pow(MS, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.0625*(m12 + m22 - m32)*C0i(cc1, m22, m12, m32, MHp2, MS2, MS2)*Cos(alp)*
   Csc(beta)*(pow(MS, 2))*(pow(Cot(beta), 2)))/(Pi2*v2) + 
 (0.0625*(m12 - m22 + m32)*C0i(cc2, m22, m12, m32, MHp2, MS2, MS2)*Cos(alp)*
   Csc(beta)*(pow(MS, 2))*(pow(Cot(beta), 2)))/(Pi2*v2) + 
 (0.0625*C0i(cc0, m22, m12, m32, MHp2, MS2, MS2)*Cos(alp)*Csc(beta)*
   (pow(MS, 2))*(-m12 - m22 + m32 + 4*(pow(MS, 2)))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MU2)*Cos(alp)*Csc(beta)*(pow(MU, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MU2, MU2)*Cos(alp)*Csc(beta)*(pow(MU, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.015625*C0i(cc0, m12, m32, m22, MA02, MA02, MC2)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*Cot(beta)*(pow(MC, 2))*(pow(Csc(beta), 2)))/
  (Pi2*v2) + (0.03125*C0i(cc0, m12, m32, m22, MHp2, MHp2, MS2)*
   ((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*Cot(beta)*(pow(MS, 2))*(pow(Csc(beta), 2)))/
  (Pi2*v2) + (0.01171875*A0i(aa0, MA02)*Cos(alp)*(3*(Mh2 - MHH2)*Cos(2*alp) + 
    (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 4*(-2*M2 + Mh2 + MHH2)*Cos(2*beta))*
   (pow(Csc(beta), 3)))/(MA02*Pi2*v2) + 
 (0.0078125*A0i(aa0, MHp2)*Cos(alp)*(3*(Mh2 - MHH2)*Cos(2*alp) + 
    (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 4*(-2*M2 + Mh2 + MHH2)*Cos(2*beta))*
   (pow(Csc(beta), 3)))/(MA02*Pi2*v2) + 
 (0.0078125*(MA02 - Mh2)*B0i(bb0, 0, MA02, Mh2)*Cos(alp)*Cos(alp - beta)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*(pow(Csc(beta), 3)))/(MA02*Pi2*v2) + 
 (0.0078125*(MA02 - Mh2)*B0i(bb0, MA02, MA02, Mh2)*Cos(alp)*Cos(alp - beta)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*(pow(Csc(beta), 3)))/(MA02*Pi2*v2) + 
 (0.0078125*A0i(aa0, Mh2)*Cos(alp)*Cos(alp - beta)*
   (-2*(M2 - MA02)*Cos(alp - 3*beta) + (Mh2 - MHH2)*Cos(3*alp - beta) + 
    (-2*M2 - 2*MA02 + 3*Mh2 + MHH2)*Cos(alp + beta))*(pow(Csc(beta), 3)))/
  (MA02*Pi2*v2) - (0.0078125*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MC2, 
    MC2, Mh2)*(pow(MC, 2))*(-4 + 3*Cos(alp)*(pow(Csc(beta), 3)) + 
    Cos(3*alp)*(pow(Csc(beta), 3))))/(Pi2*v2) - 
 (0.0625*B0i(bb0, m32, MC2, Mh2)*(pow(MC, 2))*
   (-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3))))/(Pi2*v2) - 
 (0.0625*m12*C0i(cc1, m12, m32, m22, MC2, MC2, Mh2)*(pow(MC, 2))*
   (-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3))))/(Pi2*v2) - 
 (0.125*C0i(cc0, m12, m32, m22, MC2, MC2, Mh2)*(pow(MC, 4))*
   (-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3))))/(Pi2*v2) - 
 (0.015625*(m12 - m22 - m32)*MW2*C0i(cc1, m22, m32, m12, MA02, MC2, MZ2)*
   Cos(alp - beta)*Cot(beta)*(pow(Sec(w), 2)))/(Pi2*v2) - 
 (0.015625*(3*m12 - 3*m22 + m32)*MW2*C0i(cc1, m32, m22, m12, MA02, MC2, MZ2)*
   Cos(alp - beta)*Cot(beta)*(pow(Sec(w), 2)))/(Pi2*v2) - 
 (0.015625*(m12 - m22 + m32)*MW2*C0i(cc2, m22, m32, m12, MA02, MC2, MZ2)*
   Cos(alp - beta)*Cot(beta)*(pow(Sec(w), 2)))/(Pi2*v2) - 
 (0.015625*(5*m12 + m22 - m32)*MW2*C0i(cc2, m32, m22, m12, MA02, MC2, MZ2)*
   Cos(alp - beta)*Cot(beta)*(pow(Sec(w), 2)))/(Pi2*v2) - 
 (0.015625*C0i(cc0, m22, m32, m12, MA02, MC2, MZ2)*Cos(alp - beta)*Cot(beta)*
   (m12*MW2 + m22*MW2 - m32*MW2 - 2*Mh2*(pow(MC, 2)) + 
    2*(MA02 - Mh2)*Cos(2*w)*(pow(MC, 2)) + 
    2*MA02*(-MW2 + (pow(MC, 2))))*(pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.003472222222222222*m12*C0i(cc1, m12, m32, m22, MC2, MC2, MZ2)*
   (-1 + Cos(alp)*Csc(beta))*(16*MW2*Cos(4*w) + 9*(pow(MC, 2)) + 
    Cos(2*w)*(-16*MW2 + 9*(pow(MC, 2))))*(pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.001736111111111111*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MC2, MC2, 
    MZ2)*(-1 + Cos(alp)*Csc(beta))*(16*MW2*Cos(4*w) + 9*(pow(MC, 2)) + 
    Cos(2*w)*(-16*MW2 + 9*(pow(MC, 2))))*(pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.006944444444444444*C0i(cc0, m12, m32, m22, MC2, MC2, MZ2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*(16*MW2*Cos(4*w) + 
    9*(pow(MC, 2)) + Cos(2*w)*(-16*MW2 + 9*(pow(MC, 2))))*
   (pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.001736111111111111*B0i(bb0, m32, MC2, MZ2)*(-36*(pow(MC, 2)) + 
    Cos(alp)*Csc(beta)*(27*MW2 + 32*MW2*Cos(4*w) + 18*(pow(MC, 2)) + 
      2*Cos(2*w)*(-16*MW2 + 9*(pow(MC, 2))))*(pow(Sec(w), 2)) - 
    MW2*(4 - 32*Cos(2*w) + 9*Cos(alp - 2*beta)*Csc(beta) + 
      64*(pow(Cos(2*w), 2)))*(pow(Sec(w), 2))))/(Pi2*v2) - 
 (0.0078125*C0i(cc0, m32, m22, m12, MA02, MC2, MZ2)*Cos(alp - beta)*Cot(beta)*
   (3*m12*MW2 - m22*MW2 + m32*MW2 + 3*MA02*(pow(MC, 2)) - 
    3*Mh2*(pow(MC, 2)) + (MA02 - Mh2)*Cos(4*w)*(pow(MC, 2)) + 
    Cos(2*w)*((3*m12 - m22 + m32)*MW2 + 4*MA02*(pow(MC, 2)) - 
      4*Mh2*(pow(MC, 2))) - 4*(pow(MW, 4)))*(pow(Sec(w), 4)))/
  (Pi2*v2) - (0.0625*B0i(bb0, m32, MC2, MHH2)*Cos(alp)*(pow(MC, 2))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2)))/(Pi2*v2) - 
 (0.0625*m12*C0i(cc1, m12, m32, m22, MC2, MC2, MHH2)*Cos(alp)*
   (pow(MC, 2))*(pow(Csc(beta), 3))*(pow(Sin(alp), 2)))/
  (Pi2*v2) - (0.03125*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MC2, MC2, 
    MHH2)*Cos(alp)*(pow(MC, 2))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2)))/(Pi2*v2) - 
 (0.125*C0i(cc0, m12, m32, m22, MC2, MC2, MHH2)*Cos(alp)*(pow(MC, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2)))/(Pi2*v2) + 
 (0.03125*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, MS2, MW2, MW2)*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.125*C0i(cc0, m22, m12, m32, MS2, MW2, MW2)*
   (-(m22*MW2) + Mh2*(pow(MS, 2)) + (pow(MW, 4)))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) + 
 (0.015625*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, MC2, MZ2, MZ2)*
   (pow(Sec(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2)))/(Pi2*v2) - (0.0008680555555555555*C0i(cc0, m22, m12, m32, MC2, MZ2, 
    MZ2)*(4*Cos(2*w)*(-9*m22*MW2 + 9*Mh2*(pow(MC, 2)) - 
      32*(pow(MW, 4))) + 9*(-4*m22*MW2 + 3*Mh2*(pow(MC, 2)) + 
      8*(pow(MW, 4))) + Cos(4*w)*(9*Mh2*(pow(MC, 2)) + 
      128*(pow(MW, 4))))*(pow(Sec(w), 4))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.2222222222222222*(m12 + m22 - m32)*C0i(cc1, m22, m12, m32, 0, MC2, MC2)*
   (-1 + Cos(alp)*Csc(beta))*(3*Alfas*Pi*v2 + MW2*(pow(Sin(w), 2))))/
  (Pi2*v2) + (0.2222222222222222*(m12 - m22 + m32)*
   C0i(cc2, m22, m12, m32, 0, MC2, MC2)*(-1 + Cos(alp)*Csc(beta))*
   (3*Alfas*Pi*v2 + MW2*(pow(Sin(w), 2))))/(Pi2*v2) - 
 (0.2222222222222222*C0i(cc0, m22, m12, m32, 0, MC2, MC2)*
   (-1 + Cos(alp)*Csc(beta))*(m12 + m22 - m32 - 4*(pow(MC, 2)))*
   (3*Alfas*Pi*v2 + MW2*(pow(Sin(w), 2))))/(Pi2*v2) + 
 (0.1111111111111111*B0i(bb0, m32, 0, MC2)*(-1 + Cos(alp)*Csc(beta))*
   (12*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2))))/Pi2 - 
 (0.00390625*Cos(alp)*(3*(Mh2 - MHH2)*Cos(2*alp) + 
    (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 
    4*(MA02 + (-2*M2 - MA02 + Mh2 + MHH2)*Cos(2*beta)))*
   (pow(Csc(beta), 3))*Re(A0i(aa0, MA02)))/(MA02*Pi2*v2) + 
 (0.005208333333333333*(144*Csc(beta)*(pow(MB, 2))*(pow(Cos(alp), 3))*
     (pow(Cot(beta), 2)) + 
    MA02*(61 + 4*Cos(alp)*(-12 + 2*Cos(2*w) + Cos(4*w))*Csc(beta))*
     (pow(Cot(w), 2)) + (MA02*(-93 + 61*Cos(2*w) - 20*Cos(4*w) + 
       2*Cos(6*w))*(pow(Csc(w), 2)))/2 + 8*Cos(alp)*Csc(beta)*
     (MA02*(9 - Cos(2*w) + Cos(4*w)) + 18*(pow(MB, 2))*
       (pow(Cot(beta), 2))*(pow(Sin(alp), 2))))*Re(A0i(aa0, MB2)))/
  (MA02*Pi2*v2) + (0.010416666666666666*
   (2*Cos(alp)*Csc(beta)*(64*MA02 - 24*MA02*(pow(Csc(w), 2)) + 
      (36*(pow(MC, 2))*(pow(Cot(beta), 2)) + 
        MA02*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2)))*
       (pow(Csc(w), 4))) + MA02*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 
      2*Cos(6*w))*(pow(Csc(w), 6)))*Re(A0i(aa0, MC2)))/
  (MA02*Pi2*v2*(pow(Csc(w), 4))) + 
 (0.005208333333333333*(144*Csc(beta)*(pow(MD, 2))*(pow(Cos(alp), 3))*
     (pow(Cot(beta), 2)) + 
    MA02*(61 + 4*Cos(alp)*(-12 + 2*Cos(2*w) + Cos(4*w))*Csc(beta))*
     (pow(Cot(w), 2)) + (MA02*(-93 + 61*Cos(2*w) - 20*Cos(4*w) + 
       2*Cos(6*w))*(pow(Csc(w), 2)))/2 + 8*Cos(alp)*Csc(beta)*
     (MA02*(9 - Cos(2*w) + Cos(4*w)) + 18*(pow(MD, 2))*
       (pow(Cot(beta), 2))*(pow(Sin(alp), 2))))*Re(A0i(aa0, MD2)))/
  (MA02*Pi2*v2) - (0.015625*(16*MA02 - MA02*Cos(6*w) - 
    16*MA02*Cos(alp)*Csc(beta) + MA02*Cos(alp)*Cos(6*w)*Csc(beta) - 
    10*MA02*Cos(4*w)*(-1 + Cos(alp)*Csc(beta)) + 
    Cos(2*w)*(-29*MA02 + Cos(alp)*Csc(beta)*(29*MA02 - 8*(pow(ME, 2)))) + 
    8*Cos(alp)*Csc(beta)*(pow(ME, 2)))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ME2)))/(MA02*Pi2*v2) - 
 (0.00390625*Cos(alp)*(4*MA02 + 2*(5*M2 - 6*MHH2)*Cos(2*alp) + 
    2*M2*Cos(2*(alp - 2*beta)) - Mh2*Cos(4*alp - 2*beta) + 
    MHH2*Cos(4*alp - 2*beta) - 12*M2*Cos(2*beta) - 4*MA02*Cos(2*beta) + 
    Mh2*Cos(2*beta) + 11*MHH2*Cos(2*beta))*(pow(Csc(beta), 3))*
   Re(A0i(aa0, MHH2)))/(MA02*Pi2*v2) + 
 (0.0078125*Cos(alp)*(-2*MA02 - 3*(Mh2 - MHH2)*Cos(2*alp) + 
    (-Mh2 + MHH2)*Cos(2*(alp - 2*beta)) + 8*M2*Cos(2*beta) + 
    2*MA02*Cos(2*beta) - 4*Mh2*Cos(2*beta) - 4*MHH2*Cos(2*beta) + 
    MA02*Cos(2*(beta - 2*w)) - 8*MA02*Cos(2*(beta - w)) + 16*MA02*Cos(2*w) - 
    2*MA02*Cos(4*w) - 8*MA02*Cos(2*(beta + w)) + MA02*Cos(2*(beta + 2*w)))*
   (pow(Csc(beta), 3))*Re(A0i(aa0, MHp2)))/(MA02*Pi2*v2) - 
 (0.015625*(16*MA02 - MA02*Cos(6*w) - 16*MA02*Cos(alp)*Csc(beta) + 
    MA02*Cos(alp)*Cos(6*w)*Csc(beta) - 10*MA02*Cos(4*w)*
     (-1 + Cos(alp)*Csc(beta)) + Cos(2*w)*(-29*MA02 + 
      Cos(alp)*Csc(beta)*(29*MA02 - 8*(pow(ML, 2)))) + 
    8*Cos(alp)*Csc(beta)*(pow(ML, 2)))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ML2)))/(MA02*Pi2*v2) - 
 (0.015625*(16*MA02 - MA02*Cos(6*w) - 16*MA02*Cos(alp)*Csc(beta) + 
    MA02*Cos(alp)*Cos(6*w)*Csc(beta) - 10*MA02*Cos(4*w)*
     (-1 + Cos(alp)*Csc(beta)) + Cos(2*w)*(-29*MA02 + 
      Cos(alp)*Csc(beta)*(29*MA02 - 8*(pow(MM, 2)))) + 
    8*Cos(alp)*Csc(beta)*(pow(MM, 2)))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MM2)))/(MA02*Pi2*v2) + 
 (0.005208333333333333*(144*Csc(beta)*(pow(MS, 2))*(pow(Cos(alp), 3))*
     (pow(Cot(beta), 2)) + 
    MA02*(61 + 4*Cos(alp)*(-12 + 2*Cos(2*w) + Cos(4*w))*Csc(beta))*
     (pow(Cot(w), 2)) + (MA02*(-93 + 61*Cos(2*w) - 20*Cos(4*w) + 
       2*Cos(6*w))*(pow(Csc(w), 2)))/2 + 8*Cos(alp)*Csc(beta)*
     (MA02*(9 - Cos(2*w) + Cos(4*w)) + 18*(pow(MS, 2))*
       (pow(Cot(beta), 2))*(pow(Sin(alp), 2))))*Re(A0i(aa0, MS2)))/
  (MA02*Pi2*v2) + (0.010416666666666666*
   (2*Cos(alp)*Csc(beta)*(64*MA02 - 24*MA02*(pow(Csc(w), 2)) + 
      (36*MT2*(pow(Cot(beta), 2)) + MA02*(9 - 4*Cos(2*w) + 4*Cos(4*w))*
         (pow(Cot(w), 2)))*(pow(Csc(w), 4))) + 
    MA02*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*(pow(Csc(w), 6)))*
   Re(A0i(aa0, MT2)))/(MA02*Pi2*v2*(pow(Csc(w), 4))) + 
 (0.010416666666666666*(2*Cos(alp)*Csc(beta)*
     (64*MA02 - 24*MA02*(pow(Csc(w), 2)) + 
      (36*(pow(MU, 2))*(pow(Cot(beta), 2)) + 
        MA02*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2)))*
       (pow(Csc(w), 4))) + MA02*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 
      2*Cos(6*w))*(pow(Csc(w), 6)))*Re(A0i(aa0, MU2)))/
  (MA02*Pi2*v2*(pow(Csc(w), 4))) - 
 (1.125*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb0, 0, MW2, MW2)))/(Pi2*v2) - 
 (0.1111111111111111*(-1 + Cos(alp)*Csc(beta))*
   (12*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*
   Re(B0i(bb0, MC2, 0, MC2)))/Pi2 - 
 (0.0625*Cos(alp)*Csc(beta)*(pow(MC, 2))*(pow(Cot(beta), 2))*
   Re(B0i(bb0, MC2, MA02, MC2)))/(Pi2*v2) + 
 (0.0625*(pow(MC, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(bb0, MC2, MC2, Mh2)))/(Pi2*v2) + 
 (0.0625*Cos(alp)*(pow(MC, 2))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, MC2, MC2, MHH2)))/(Pi2*v2) - 
 (0.003472222222222222*(-1 + Cos(alp)*Csc(beta))*
   (16*MW2*Cos(4*w) + 9*(pow(MC, 2)) + 
    Cos(2*w)*(-16*MW2 + 9*(pow(MC, 2))))*(pow(Sec(w), 2))*
   Re(B0i(bb0, MC2, MC2, MZ2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Csc(beta)*(pow(MS, 2))*(pow(Cot(beta), 2))*
   Re(B0i(bb0, MC2, MHp2, MS2)))/(Pi2*v2) - 
 (0.125*(-1 + Cos(alp)*Csc(beta))*(pow(MS, 2))*
   Re(B0i(bb0, MC2, MS2, MW2)))/(Pi2*v2) + 
 (0.75*Cos(alp)*(pow(MB, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, MB2, MB2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.75*Cos(alp)*(pow(MC, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, Mh2, MC2, MC2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*(pow(MD, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, MD2, MD2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.25*Cos(alp)*Csc(beta)*(pow(ME, 4))*
   (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, Mh2, ME2, ME2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.25*Cos(alp)*Csc(beta)*(pow(ML, 4))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, ML2, ML2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.25*Cos(alp)*Csc(beta)*(pow(MM, 4))*
   (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, Mh2, MM2, MM2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*(pow(MS, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, MS2, MS2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.75*Cos(alp)*(pow(MT, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, Mh2, MT2, MT2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*(pow(MU, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, MU2, MU2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.75*Cos(alp)*(pow(MB, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, MB2, MB2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*(pow(MC, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, MHH2, MC2, MC2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.75*Cos(alp)*(pow(MD, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, MD2, MD2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*Cos(alp)*Csc(beta)*(pow(ME, 4))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, MHH2, ME2, ME2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.25*Cos(alp)*Csc(beta)*(pow(ML, 4))*
   (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, ML2, ML2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*Cos(alp)*Csc(beta)*(pow(MM, 4))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, MHH2, MM2, MM2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.75*Cos(alp)*(pow(MS, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, MS2, MS2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*(pow(MT, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, MHH2, MT2, MT2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.75*Cos(alp)*(pow(MU, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, MU2, MU2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.5*MW2*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*Re(B0i(bb0, MW2, 0, MW2)))/
  (Pi2*v2) + (0.375*(-MT2 + MW2)*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*
   (pow(Csc(w), 2))*Re(B0i(bb0, MW2, MB2, MT2)))/(Pi2*v2) - 
 (0.375*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*
   (-1 + (pow(Cot(w), 2)))*Re(B0i(bb0, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*(MW2 - (pow(MU, 2)))*
   (pow(Csc(w), 2))*Re(B0i(bb0, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.125*MW2*Cos(2*w)*(pow(Csc(w), 2))*
   (-1 + Cos(alp)*Csc(beta)*(pow(Sin(alp - beta), 2)))*
   Re(B0i(bb0, MW2, Mh2, MW2)))/(Pi2*v2) + 
 (0.125*MW2*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb0, MW2, MHH2, MW2)))/(Pi2*v2) - 
 (0.015625*MW2*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*(pow(Sec(w), 2))*
   Re(B0i(bb0, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MD, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*(-1 + Cos(alp)*Csc(beta))*(pow(ME, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, ME2, ME2)))/(Pi2*v2) - 
 (0.125*MW2*(pow(Csc(w), 2))*(-1 + Cos(alp)*Csc(beta)*
     (pow(Sin(alp - beta), 2)))*Re(B0i(bb0, MZ2, Mh2, MZ2)))/(Pi2*v2) - 
 (0.125*MW2*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb0, MZ2, MHH2, MZ2)))/(Pi2*v2) + 
 (0.0625*(-1 + Cos(alp)*Csc(beta))*(pow(ML, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*(-1 + Cos(alp)*Csc(beta))*(pow(MM, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MS, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.1875*MT2*(-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MU, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.0625*MW2*(5 + 9*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MB2, MB2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MD2, MD2)))/(Pi2*v2) + 
 (0.5*(-1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, ME2, ME2)))/(Pi2*v2) - 
 (0.5*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MHp2, MHp2)))/(Pi2*v2) + 
 (0.5*(-1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, ML2, ML2)))/(Pi2*v2) + 
 (0.5*(-1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MM2, MM2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MT2, MT2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.5*(2 + 3*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.25*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, 0, ME2)))/(Pi2*v2) + 
 (0.25*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, 0, ML2)))/(Pi2*v2) + 
 (0.25*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, 0, MM2)))/(Pi2*v2) - 
 (1.*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*Re(B0i(bb00, MW2, 0, MW2)))/
  (Pi2*v2) - (0.125*Cos(alp)*Csc(beta)*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, MA02, MHp2)))/(Pi2*v2) + 
 (0.75*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, MB2, MT2)))/(Pi2*v2) + 
 (0.75*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.75*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, MD2, MU2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MW2, Mh2, MHp2)))/(Pi2*v2) - 
 (0.125*Cos(2*w)*(pow(Csc(w), 2))*
   (-1 + Cos(alp)*Csc(beta)*(pow(Sin(alp - beta), 2)))*
   Re(B0i(bb00, MW2, Mh2, MW2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MW2, MHH2, MHp2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MW2, MHH2, MW2)))/(Pi2*v2) - 
 (0.125*(2 + 5*Cos(2*w) + 2*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MW2, MW2, MZ2)))/(Pi2*v2) - 
 (0.375*(-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, 0, 0)))/(Pi2*v2) + 
 (0.125*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MA02, Mh2)))/(Pi2*v2) + 
 (0.125*Cos(alp)*Csc(beta)*(pow(Cot(w), 2))*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb00, MZ2, MA02, MHH2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Cot(w), 2))*Re(B0i(bb00, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.041666666666666664*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MC2, MC2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Cot(w), 2))*Re(B0i(bb00, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*(pow(Cot(w), 2))*(-1 + Cos(alp)*Csc(beta)*
     (pow(Sin(alp - beta), 2)))*Re(B0i(bb00, MZ2, Mh2, MZ2)))/(Pi2*v2) + 
 (0.125*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MHH2, MZ2)))/(Pi2*v2) + 
 (0.125*Cos(alp)*Csc(beta)*(pow(Cos(2*w), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MHp2, MHp2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Cot(w), 2))*Re(B0i(bb00, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.041666666666666664*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.041666666666666664*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.0625*(7 + 8*Cos(2*w) + 3*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Cot(w), 2))*Re(B0i(bb00, MZ2, MW2, MW2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MB2, MB2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MC2, MC2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MD2, MD2)))/(Pi2*v2) - 
 (0.5*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, ME2, ME2)))/(Pi2*v2) - 
 (0.5*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, ML2, ML2)))/(Pi2*v2) - 
 (0.5*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MM2, MM2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MS2, MS2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MT2, MT2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.25*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.03125*MW2*Cos(alp - beta)*(MHH2*Cos(2*alp - beta) + 
    (-2*Mh2 + MHH2)*Cos(beta))*Csc(beta)*(pow(Sec(w), 2))*
   Re(B0i(bb1, Mh2, MA02, MZ2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MB, 2))*(-2*Mh2 + 2*MHH2 - 
    Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*(pow(Csc(beta), 3)))*
   Re(B0i(bb1, Mh2, MB2, MB2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MC, 2))*(-2*Mh2 + 2*MHH2 - 
    Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*(pow(Csc(beta), 3)))*
   Re(B0i(bb1, Mh2, MC2, MC2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MD, 2))*(-2*Mh2 + 2*MHH2 - 
    Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*(pow(Csc(beta), 3)))*
   Re(B0i(bb1, Mh2, MD2, MD2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.125*(pow(ME, 2))*(Mh2 - MHH2 + MHH2*Cos(alp)*Csc(beta)*
     (pow(Sec(beta), 2))*(pow(Sin(alp), 2)))*
   Re(B0i(bb1, Mh2, ME2, ME2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0625*MW2*Cos(alp - beta)*(MHH2*Cos(2*alp - beta) + 
    (-2*Mh2 + MHH2)*Cos(beta))*Csc(beta)*Re(B0i(bb1, Mh2, MHp2, MW2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.125*(pow(ML, 2))*
   (Mh2 - MHH2 + MHH2*Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(bb1, Mh2, ML2, ML2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.125*(pow(MM, 2))*
   (Mh2 - MHH2 + MHH2*Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(bb1, Mh2, MM2, MM2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.1875*(pow(MS, 2))*
   (-2*Mh2 + 2*MHH2 - Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*
     (pow(Csc(beta), 3)))*Re(B0i(bb1, Mh2, MS2, MS2)))/
  ((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*MT2*(-2*Mh2 + 2*MHH2 - Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*
     (pow(Csc(beta), 3)))*Re(B0i(bb1, Mh2, MT2, MT2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.1875*(pow(MU, 2))*
   (-2*Mh2 + 2*MHH2 - Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*
     (pow(Csc(beta), 3)))*Re(B0i(bb1, Mh2, MU2, MU2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.375*MHH2*Cos(alp)*(pow(MB, 2))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, MB2, MB2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*(pow(MC, 2))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb1, MHH2, MC2, MC2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.375*MHH2*Cos(alp)*(pow(MD, 2))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, MD2, MD2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.125*MHH2*Cos(alp)*Csc(beta)*(pow(ME, 2))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb1, MHH2, ME2, ME2)))/
  ((-Mh2 + MHH2)*Pi2*v2) - (0.125*MHH2*Cos(alp)*Csc(beta)*(pow(ML, 2))*
   (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, ML2, ML2)))/((-Mh2 + MHH2)*Pi2*v2) - 
 (0.125*MHH2*Cos(alp)*Csc(beta)*(pow(MM, 2))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb1, MHH2, MM2, MM2)))/
  ((-Mh2 + MHH2)*Pi2*v2) - (0.375*MHH2*Cos(alp)*(pow(MS, 2))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, MS2, MS2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*MT2*Cos(alp)*(pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, MT2, MT2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*(pow(MU, 2))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb1, MHH2, MU2, MU2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.125*MW2*(-1 + Cos(alp)*Csc(beta))*
   (-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ME2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, 0, ML2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, 0, MM2)))/(Pi2*v2) - 
 (0.25*MW2*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*Re(B0i(bb1, MW2, 0, MW2)))/
  (Pi2*v2) + (0.375*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MB2, MT2)))/(Pi2*v2) - 
 (0.375*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.25*MW2*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb1, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.1875*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, 0, 0)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.25*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Cos(w), 2))*
   (pow(Cot(w), 2))*Re(B0i(bb1, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.1111111111111111*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*
   (12*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*
   Re(B0i(dbb0, MC2, 0, MC2)))/Pi2 - 
 (0.125*(pow(MC, 4))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb0, MC2, MC2, Mh2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*(pow(MC, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(dbb0, MC2, MC2, MHH2)))/(Pi2*v2) + 
 (0.006944444444444444*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*
   (16*MW2*Cos(4*w) + 9*(pow(MC, 2)) + 
    Cos(2*w)*(-16*MW2 + 9*(pow(MC, 2))))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MC2, MC2, MZ2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Csc(beta)*(pow(MC, 2))*((pow(MC, 2)) - 
    (pow(MS, 2)))*(pow(Cot(beta), 2))*Re(B0i(dbb0, MC2, MHp2, MS2)))/
  (Pi2*v2) + (0.25*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*
   (pow(MS, 2))*Re(B0i(dbb0, MC2, MS2, MW2)))/(Pi2*v2) - 
 (0.0009765625*Cos(alp)*(pow((2*MA02 - Mh2)*Cos(alp - 3*beta) + 
      (4*M2 - 2*MA02 - 3*Mh2)*Cos(alp + beta), 2))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*Re(B0i(dbb0, Mh2, MA02, MA02)))/(Pi2*v2) + 
 (0.75*(pow(MB, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.75*(pow(MC, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.75*(pow(MD, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.25*(pow(ME, 4))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb0, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.0087890625*(16*(pow(Mh, 4)) - 
    Cos(alp)*(pow(M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
        (2*M2 - 3*Mh2)*Cos(alp + beta), 2))*(pow(Csc(beta), 3))*
     (pow(Sec(beta), 2)))*Re(B0i(dbb0, Mh2, Mh2, Mh2)))/(Pi2*v2) - 
 (0.0078125*Cos(alp)*(pow(Cos(alp - beta), 2))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*(pow((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + 
      M2*Sin(2*beta), 2))*Re(B0i(dbb0, Mh2, Mh2, MHH2)))/(Pi2*v2) - 
 (0.00390625*Cos(alp)*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   (pow(Sin(alp - beta), 2))*
   (pow((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta), 2))*
   Re(B0i(dbb0, Mh2, MHH2, MHH2)))/(Pi2*v2) - 
 (0.001953125*Cos(alp)*(pow((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + 
      (-4*M2 + 3*Mh2 + 2*MHp2)*Cos(alp + beta), 2))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*Re(B0i(dbb0, Mh2, MHp2, MHp2)))/(Pi2*v2) - 
 (0.0625*Cos(alp)*Csc(beta)*(MHp2*(MHp2 - MW2) - Mh2*(2*MHp2 + MW2) + 
    (pow(Mh, 4)))*(pow(Cos(alp - beta), 2))*
   Re(B0i(dbb0, Mh2, MHp2, MW2)))/(Pi2*v2) + 
 (0.25*(pow(ML, 4))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb0, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.25*(pow(MM, 4))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb0, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.75*(pow(MS, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.75*(pow(MT, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.75*(pow(MU, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MU2, MU2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MB2, MB2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MD2, MD2)))/(Pi2*v2) + 
 (1.*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, ME2, ME2)))/(Pi2*v2) - 
 (0.5*MW2*Cos(alp)*Csc(beta)*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MHp2, MHp2)))/(Pi2*v2) + 
 (1.*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, ML2, ML2)))/(Pi2*v2) + 
 (1.*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MM2, MM2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MT2, MT2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (1.5*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MW2, MW2)))/(Pi2*v2) - 
 (0.1111111111111111*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*
   (12*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*
   Re(B0i(dbb1, MC2, 0, MC2)))/Pi2 - 
 (0.125*Cos(alp)*Csc(beta)*(pow(MC, 4))*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MC2, MA02, MC2)))/(Pi2*v2) + 
 (0.125*(pow(MC, 4))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, MC2, MC2, Mh2)))/(Pi2*v2) + 
 (0.125*Cos(alp)*(pow(MC, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(dbb1, MC2, MC2, MHH2)))/(Pi2*v2) + 
 (0.006944444444444444*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*
   (8*MW2*Cos(4*w) + 9*(2*MW2 + (pow(MC, 2))) + 
    Cos(2*w)*(-8*MW2 + 9*(pow(MC, 2))))*(pow(Sec(w), 2))*
   Re(B0i(dbb1, MC2, MC2, MZ2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Csc(beta)*(pow(MC, 2))*((pow(MC, 2)) + 
    (pow(MS, 2)))*(pow(Cot(beta), 2))*Re(B0i(dbb1, MC2, MHp2, MS2)))/
  (Pi2*v2) + (0.125*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*
   (2*MW2 + (pow(MC, 2)) + (pow(MS, 2)))*
   Re(B0i(dbb1, MC2, MS2, MW2)))/(Pi2*v2) - 
 (0.0625*Mh2*MW2*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Re(B0i(dbb1, Mh2, MA02, MZ2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MB, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MC, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MD, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ME, 2))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb1, Mh2, ME2, ME2)))/(Pi2*v2) - 
 (0.125*Mh2*MW2*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   Re(B0i(dbb1, Mh2, MHp2, MW2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ML, 2))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(MM, 2))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MS, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*Mh2*MT2*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MU, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*Mh2*MW2*(-1 + Cos(alp)*Csc(beta)*(pow(Sin(alp - beta), 2)))*
   Re(B0i(dbb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*Mh2*MW2*(pow(Sec(w), 2))*
   (-1 + Cos(alp)*Csc(beta)*(pow(Sin(alp - beta), 2)))*
   Re(B0i(dbb1, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (0.046875*C0i(cc0, m22, m12, m32, MC2, Mh2, Mh2)*(pow(MC, 2))*
   (-4*Mh2 - (M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
      (2*M2 - 3*Mh2)*Cos(alp + beta))*(pow(Cos(alp), 2))*
     (pow(Csc(beta), 3))*Sec(beta)))/(Pi2*v2) + 
 (0.0625*B0i(bb0, m32, MHp2, MS2)*
   (Cos(alp)*Cot(beta)*(MW2*Cos(beta) + 2*Cot(beta)*Csc(beta)*
       (pow(MS, 2))) + MW2*Cos(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.03125*B0i(bb0, m32, MA02, MC2)*
   (Cos(alp)*Cot(beta)*(2*Cot(beta)*Csc(beta)*(pow(MC, 2)) + 
      MW2*Cos(beta)*(pow(Sec(w), 2))) + MW2*Cos(beta)*
     (pow(Sec(w), 2))*Sin(alp)))/(Pi2*v2) + 
 (0.00390625*Re(A0i(aa0, Mh2))*(8*MA02 + 5*M2*Cos(3*alp)*
     (pow(Csc(beta), 3)) + Cos(alp)*(5*M2 - 4*MA02 + 
      2*M2*Cos(2*(alp - 2*beta)) + (-Mh2 + MHH2)*Cos(4*alp - 2*beta) + 
      12*M2*Cos(2*beta) + 4*MA02*Cos(2*beta) - 11*Mh2*Cos(2*beta) - 
      MHH2*Cos(2*beta))*(pow(Csc(beta), 3)) - 
    3*Mh2*Csc(alp)*(pow(Csc(beta), 3))*Sin(4*alp)))/(MA02*Pi2*v2) - 
 (0.03125*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, MS2, MW2, MW2)*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.015625*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, MC2, MZ2, MZ2)*
   (pow(Sec(w), 2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.03125*(Mh2 - MHH2)*A0i(aa0, MW2)*Cos(alp)*Cot(beta)*Csc(beta)*
   Sin(2*(alp - beta)))/(MA02*Pi2*v2) - 
 (0.046875*(Mh2 - MHH2)*A0i(aa0, MZ2)*Cos(alp)*Cot(beta)*Csc(beta)*
   Sin(2*(alp - beta)))/(MA02*Pi2*v2) - 
 (0.03125*MW2*B0i(bb1, MA02, Mh2, MZ2)*Cos(alp)*Cot(beta)*Csc(beta)*
   (pow(Sec(w), 2))*Sin(2*(alp - beta)))/(Pi2*v2) + 
 (0.03125*MW2*B0i(bb1, MA02, MHH2, MZ2)*Cos(alp)*Cot(beta)*Csc(beta)*
   (pow(Sec(w), 2))*Sin(2*(alp - beta)))/(Pi2*v2) + 
 (0.0078125*Mh2*B0i(bb0, 0, Mh2, MZ2)*Cos(alp)*(MA02 - Mh2 + 2*MW2 + 
    (MA02 - Mh2)*Cos(2*w))*Cot(beta)*Csc(beta)*(pow(Sec(w), 2))*
   Sin(2*(alp - beta)))/(MA02*Pi2*v2) - 
 (0.0078125*B0i(bb0, MA02, MHH2, MZ2)*Cos(alp)*(-(MHH2*(MHH2 - 2*MW2)) + 
    MA02*(MHH2 + 2*MW2) + (MA02 - MHH2)*MHH2*Cos(2*w))*Cot(beta)*Csc(beta)*
   (pow(Sec(w), 2))*Sin(2*(alp - beta)))/(MA02*Pi2*v2) + 
 (0.0078125*MHH2*B0i(bb0, 0, MHH2, MZ2)*Cos(alp)*
   (-MA02 + MHH2 - 2*MW2 + (-MA02 + MHH2)*Cos(2*w))*Cot(beta)*Csc(beta)*
   (pow(Sec(w), 2))*Sin(2*(alp - beta)))/(MA02*Pi2*v2) + 
 (0.0078125*B0i(bb0, MA02, Mh2, MZ2)*Cos(alp)*Cot(beta)*Csc(beta)*
   (2*Mh2*MW2 + MA02*(Mh2 + 2*MW2) + (MA02 - Mh2)*Mh2*Cos(2*w) - 
    (pow(Mh, 4)))*(pow(Sec(w), 2))*Sin(2*(alp - beta)))/
  (MA02*Pi2*v2) + (0.015625*Csc(beta)*((MA02 - Mh2)*(MA02 - MHH2) - 
    (MA02 + MHH2)*MW2*(pow(Sec(w), 2)))*Re(B0i(bb0, MHH2, MA02, MZ2))*
   Sin(alp)*Sin(2*(alp - beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.03125*(Mh2*(MHH2 - MHp2) + MHp2*(MHp2 - MW2) - MHH2*(MHp2 + MW2))*
   Csc(beta)*Re(B0i(bb0, MHH2, MHp2, MW2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - 
 (0.015625*Csc(beta)*(Mh2*MHH2 - 2*MHH2*MW2 + 12*(pow(MW, 4)))*
   Re(B0i(bb0, MHH2, MW2, MW2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.0009765625*Csc(beta)*
   (3*Mh2*MHH2 - 8*MHH2*MW2 + 4*MHH2*(Mh2 - 2*MW2)*Cos(2*w) + 
    Mh2*MHH2*Cos(4*w) + 96*(pow(MW, 4)))*(pow(Sec(w), 4))*
   Re(B0i(bb0, MHH2, MZ2, MZ2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.03125*MHH2*MW2*Csc(beta)*(pow(Sec(w), 2))*
   Re(B0i(bb1, MHH2, MA02, MZ2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.0625*MHH2*MW2*Csc(beta)*
   Re(B0i(bb1, MHH2, MHp2, MW2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.0625*MHH2*MW2*Csc(beta)*
   Re(B0i(bb1, MHH2, MW2, MW2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.03125*MHH2*MW2*Csc(beta)*(pow(Sec(w), 2))*
   Re(B0i(bb1, MHH2, MZ2, MZ2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.03125*Csc(beta)*Re(B0i(bb0, Mh2, MHp2, MW2))*
   (2*(Mh2 - MHH2)*MW2*Cos(alp)*(pow(Cos(alp - beta), 2)) + 
    (MHp2*(MHH2 - MHp2 + MW2) + Mh2*(-MHH2 + MHp2 + MW2))*Sin(alp)*
     Sin(2*(alp - beta))))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.015625*MW2*Csc(beta)*Re(B0i(bb0, Mh2, MA02, MZ2))*
   (-2*Cos(alp)*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 2)) + 
    (pow(Mh2 - MHH2, -1))*(pow(MW, -2))*((MA02 - Mh2)*(MA02 - MHH2) - 
      (MA02 + Mh2)*MW2*(pow(Sec(w), 2)))*Sin(alp)*Sin(2*(alp - beta))))/
  (Pi2*v2) + (0.015625*MW2*Re(B0i(bb0, Mh2, MZ2, MZ2))*
   (-8*(pow(Csc(2*w), 2)) - Csc(beta)*(pow(Csc(w), 2))*
     (pow(Sec(w), 2))*(-2*Cos(alp)*(pow(Sin(alp - beta), 2)) - 
      ((pow(Mh2 - MHH2, -1))*(pow(MW, -2))*(3*Mh2*MHH2 - 8*Mh2*MW2 + 
         4*Mh2*(MHH2 - 2*MW2)*Cos(2*w) + Mh2*MHH2*Cos(4*w) + 
         96*(pow(MW, 4)))*(pow(Sec(w), 2))*Sin(alp)*
        Sin(2*(alp - beta)))/16)))/(Pi2*v2*(pow(Csc(w), 2))) + 
 (0.0625*MW2*Re(B0i(bb1, Mh2, MW2, MW2))*(2*(Mh2 - MHH2) + 
    Sin(alp - beta)*(2*Mh2 - MHH2 + MHH2*Csc(beta)*Sin(2*alp - beta))))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.03125*MW2*(pow(Sec(w), 2))*
   Re(B0i(bb1, Mh2, MZ2, MZ2))*(2*(Mh2 - MHH2) + 
    Sin(alp - beta)*(2*Mh2 - MHH2 + MHH2*Csc(beta)*Sin(2*alp - beta))))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.015625*Re(B0i(bb0, Mh2, MW2, MW2))*
   (-4*MW2 + (pow(Mh2 - MHH2, -1))*Sin(alp - beta)*
     (Mh2*(MHH2 - 4*MW2) + 2*MW2*(MHH2 + 6*MW2) + 
      Csc(beta)*(Mh2*MHH2 - 2*MHH2*MW2 + 12*(pow(MW, 4)))*
       Sin(2*alp - beta))))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ME2, ME2))*(Cos(alp) - Sin(beta)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ML2, ML2))*(Cos(alp) - Sin(beta)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MM2, MM2))*(Cos(alp) - Sin(beta)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, MC2, MC2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, ME2, ME2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, ML2, ML2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, MM2, MM2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*
   Csc(beta)*(pow(Csc(w), 2))*Re(B0i(bb1, MZ2, MT2, MT2))*
   (Cos(alp) - Sin(beta)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, MU2, MU2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.03125*Cos(alp)*(-(pow(MA02 - Mh2, 2)) + 
    (MA02 + Mh2)*MW2*(pow(Sec(w), 2)))*
   (pow(Cos(alp)*Cot(beta) + Sin(alp), 2))*Re(B0i(dbb0, Mh2, MA02, MZ2))*
   Sin(beta))/(Pi2*v2) + 
 (0.0078125*Csc(beta)*(-2*Mh2*MW2 + (pow(Mh, 4)) + 12*(pow(MW, 4)))*
   Re(B0i(dbb0, Mh2, MW2, MW2))*(-2*Cos(alp) + Cos(alp - 2*beta) + 
    Cos(3*alp - 2*beta) + 4*Sin(beta)))/(Pi2*v2) + 
 (0.00048828125*Csc(beta)*(-8*Mh2*MW2 + 4*Mh2*(Mh2 - 2*MW2)*Cos(2*w) + 
    3*(pow(Mh, 4)) + Cos(4*w)*(pow(Mh, 4)) + 96*(pow(MW, 4)))*
   (pow(Sec(w), 4))*Re(B0i(dbb0, Mh2, MZ2, MZ2))*
   (-2*Cos(alp) + Cos(alp - 2*beta) + Cos(3*alp - 2*beta) + 4*Sin(beta)))/
  (Pi2*v2) + (0.0078125*Re(A0i(aa0, MZ2))*
   (2*MA02*(-13 + 24*Cos(alp)*Csc(beta)*(pow(Cos(w), 2)))*
     (pow(Cot(w), 2)) - 2*MA02*(-13 + (5 + 6*Cos(4*w))*
       (pow(Csc(w), 2))) + Cos(alp)*(pow(Csc(beta), 2))*
     ((Mh2 - MHH2)*Sin(2*alp - 3*beta) + (Mh2 - MHH2)*Sin(2*alp - beta) - 
      4*MA02*(7 + 6*Cos(2*w))*Sin(beta))))/(MA02*Pi2*v2) + 
 (0.03125*C0i(cc0, m22, m12, m32, MC2, MHH2, MHH2)*(pow(MC, 2))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*Sec(beta)*Sin(alp - beta)*
   ((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) + 
 (0.005859375*Cos(alp - beta)*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh2)*Cos(3*alp - beta) + (2*M2 - 3*Mh2)*Cos(alp + beta))*
   (pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, Mh2, Mh2, Mh2))*Sin(alp)*((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.005859375*Cos(alp - beta)*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh2)*Cos(3*alp - beta) + (2*M2 - 3*Mh2)*Cos(alp + beta))*
   (pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, Mh2, Mh2))*Sin(alp)*((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.015625*C0i(cc0, m22, m12, m32, MC2, Mh2, MHH2)*Cos(alp - beta)*
   (pow(MC, 2))*(pow(Csc(beta), 3))*Sec(beta)*Sin(2*alp)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) - 
 (0.015625*C0i(cc0, m32, m12, m22, MC2, Mh2, MHH2)*Cos(alp - beta)*
   (pow(MC, 2))*(pow(Csc(beta), 3))*Sec(beta)*Sin(2*alp)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) - 
 (0.001953125*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, Mh2, Mh2, MHH2))*Sin(alp)*Sin(2*(alp - beta))*
   (-9*M2*Mh2 - 9*M2*MHH2 + 5*Mh2*MHH2 + 8*(pow(M2, 2)) + 
    2*(pow(Mh, 4)) + 2*(pow(MHH2, 2)) - 
    Cos(4*alp)*(5*Mh2*MHH2 - 9*M2*(Mh2 + MHH2) + 9*(pow(M2, 2)) + 
      2*(pow(Mh, 4)) + 2*(pow(MHH2, 2))) + 
    (pow(M2, 2))*(pow(Cos(beta), 4)) - 6*(pow(M2, 2))*
     (pow(Cos(beta), 2))*(pow(Sin(beta), 2)) + 
    (pow(M2, 2))*(pow(Sin(beta), 4)) - 2*M2*(Mh2 - MHH2)*Sin(2*alp)*
     Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.001953125*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, Mh2, MHH2))*Sin(alp)*Sin(2*(alp - beta))*
   (-9*M2*Mh2 - 9*M2*MHH2 + 5*Mh2*MHH2 + 8*(pow(M2, 2)) + 
    2*(pow(Mh, 4)) + 2*(pow(MHH2, 2)) - 
    Cos(4*alp)*(5*Mh2*MHH2 - 9*M2*(Mh2 + MHH2) + 9*(pow(M2, 2)) + 
      2*(pow(Mh, 4)) + 2*(pow(MHH2, 2))) + 
    (pow(M2, 2))*(pow(Cos(beta), 4)) - 6*(pow(M2, 2))*
     (pow(Cos(beta), 2))*(pow(Sin(beta), 2)) + 
    (pow(M2, 2))*(pow(Sin(beta), 4)) - 2*M2*(Mh2 - MHH2)*Sin(2*alp)*
     Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0234375*Csc(beta)*(pow(Csc(2*beta), 2))*Re(B0i(bb0, Mh2, MHH2, MHH2))*
   Sin(alp)*Sin(alp - beta)*((-3*M2 + Mh2 + 2*MHH2)*Sin(2*alp) - 
    M2*Sin(2*beta))*(M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0234375*Csc(beta)*(pow(Csc(2*beta), 2))*
   Re(B0i(bb0, MHH2, MHH2, MHH2))*Sin(alp)*Sin(alp - beta)*
   ((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta))*
   (M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0009765625*((2*MA02 - Mh2)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh2)*Cos(alp + beta))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh2, MA02, MA02))*Sin(alp)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0009765625*((-2*MA02 + Mh2)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh2)*Cos(alp + beta))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MA02, MA02))*Sin(alp)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0078125*(MA02 - MHH2)*B0i(bb0, 0, MA02, MHH2)*Cos(alp)*
   (pow(Csc(beta), 3))*Sin(alp - beta)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(MA02*Pi2*v2) - 
 (0.0078125*(MA02 - MHH2)*B0i(bb0, MA02, MA02, MHH2)*Cos(alp)*
   (pow(Csc(beta), 3))*Sin(alp - beta)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(MA02*Pi2*v2) + 
 (0.0078125*A0i(aa0, MHH2)*Cos(alp)*(pow(Csc(beta), 3))*Sin(alp - beta)*
   (-2*(M2 - MA02)*Sin(alp - 3*beta) + (Mh2 - MHH2)*Sin(3*alp - beta) + 
    (-2*M2 - 2*MA02 + Mh2 + 3*MHH2)*Sin(alp + beta)))/(MA02*Pi2*v2) - 
 (0.001953125*((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, Mh2, MHp2, MHp2))*Sin(alp)*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.001953125*((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, MHp2, MHp2))*Sin(alp)*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.00390625*(pow(Csc(w), 2))*Re(A0i(aa0, MW2))*
   (44 - 2*Cos(2*w) + 60*Cos(4*w) - 6*Cos(6*w) + Cos(alp)*(pow(MA02, -1))*
     (pow(Csc(beta), 2))*(2*(Mh2 - MHH2)*Sin(2*alp - 3*beta) + 
      2*(Mh2 - MHH2)*Sin(2*alp - beta) - 44*MA02*Sin(beta) + 
      3*MA02*Sin(beta - 6*w) - 30*MA02*Sin(beta - 4*w) - 
      Mh2*Sin(2*alp - 3*beta - 2*w) + MHH2*Sin(2*alp - 3*beta - 2*w) - 
      Mh2*Sin(2*alp - beta - 2*w) + MHH2*Sin(2*alp - beta - 2*w) + 
      MA02*Sin(beta - 2*w) - Mh2*Sin(2*alp - 3*beta + 2*w) + 
      MHH2*Sin(2*alp - 3*beta + 2*w) - Mh2*Sin(2*alp - beta + 2*w) + 
      MHH2*Sin(2*alp - beta + 2*w) + MA02*Sin(beta + 2*w) - 
      30*MA02*Sin(beta + 4*w) + 3*MA02*Sin(beta + 6*w))))/(Pi2*v2)
;
}

ComplexType dKappahtautau(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32)
{
return 0. + (0.046875*C0i(cc0, m12, m32, m22, Mh2, Mh2, ML2)*(pow(ML, 2))*
   (4*Mh2 + (M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
      (2*M2 - 3*Mh2)*Cos(alp + beta))*Csc(beta)*(pow(Sec(beta), 3))*
     (pow(Sin(alp), 2))))/(Pi2*v2) - 
 (0.0078125*(m12 + m22 - m32)*C0i(cc1, m22, m12, m32, Mh2, ML2, ML2)*
   (pow(ML, 2))*(pow(Sec(beta), 3))*(3*Cos(beta) + Cos(3*beta) + 
    4*(pow(Sin(alp), 3))))/(Pi2*v2) + 
 (0.0078125*(m12 - m22 + m32)*C0i(cc2, m22, m12, m32, Mh2, ML2, ML2)*
   (pow(ML, 2))*(pow(Sec(beta), 3))*(3*Cos(beta) + Cos(3*beta) + 
    4*(pow(Sin(alp), 3))))/(Pi2*v2) + 
 (0.0078125*C0i(cc0, m22, m12, m32, Mh2, ML2, ML2)*(pow(ML, 2))*
   (-m12 - m22 + m32 + 4*(pow(ML, 2)))*(pow(Sec(beta), 3))*
   (3*Cos(beta) + Cos(3*beta) + 4*(pow(Sin(alp), 3))))/(Pi2*v2) + 
 (0.0625*B0i(bb0, m32, Mh2, ML2)*(pow(ML, 2))*
   (1 + (pow(Sec(beta), 3))*(pow(Sin(alp), 3))))/(Pi2*v2) + 
 (0.03125*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, 0, MW2, MW2)*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) + 
 (0.015625*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, ML2, MZ2, MZ2)*
   (pow(Sec(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2)))/(Pi2*v2) - (0.0078125*C0i(cc0, m22, m12, m32, ML2, MZ2, MZ2)*
   (-4*m22*MW2 + 3*Mh2*(pow(ML, 2)) + 
    4*Cos(2*w)*(-(MW2*(m22 + 16*MW2)) + Mh2*(pow(ML, 2))) + 
    40*(pow(MW, 4)) + Cos(4*w)*(Mh2*(pow(ML, 2)) + 
      32*(pow(MW, 4))))*(pow(Sec(w), 4))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.0625*(pow(ML, 2))*(1 + (pow(Sec(beta), 3))*
     (pow(Sin(alp), 3)))*Re(B0i(bb0, ML2, Mh2, ML2)))/(Pi2*v2) - 
 (0.25*(pow(ME, 4))*(1 + (pow(Sec(beta), 3))*(pow(Sin(alp), 3)))*
   Re(B0i(dbb0, Mh2, ME2, ME2)))/(Pi2*v2) - 
 (0.25*(pow(ML, 4))*(1 + (pow(Sec(beta), 3))*(pow(Sin(alp), 3)))*
   Re(B0i(dbb0, Mh2, ML2, ML2)))/(Pi2*v2) - 
 (0.25*(pow(MM, 4))*(1 + (pow(Sec(beta), 3))*(pow(Sin(alp), 3)))*
   Re(B0i(dbb0, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.25*(pow(ML, 4))*(1 + (pow(Sec(beta), 3))*(pow(Sin(alp), 3)))*
   Re(B0i(dbb0, ML2, Mh2, ML2)))/(Pi2*v2) - 
 (0.125*Mh2*(pow(ME, 2))*(1 + (pow(Sec(beta), 3))*
     (pow(Sin(alp), 3)))*Re(B0i(dbb1, Mh2, ME2, ME2)))/(Pi2*v2) - 
 (0.125*Mh2*(pow(ML, 2))*(1 + (pow(Sec(beta), 3))*
     (pow(Sin(alp), 3)))*Re(B0i(dbb1, Mh2, ML2, ML2)))/(Pi2*v2) - 
 (0.125*Mh2*(pow(MM, 2))*(1 + (pow(Sec(beta), 3))*
     (pow(Sin(alp), 3)))*Re(B0i(dbb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.125*(pow(ML, 4))*(1 + (pow(Sec(beta), 3))*(pow(Sin(alp), 3)))*
   Re(B0i(dbb1, ML2, Mh2, ML2)))/(Pi2*v2) + 
 (0.01171875*A0i(aa0, MA02)*(3*(Mh2 - MHH2)*Cos(2*alp) + 
    (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 4*(-2*M2 + Mh2 + MHH2)*Cos(2*beta))*
   (pow(Sec(beta), 3))*Sin(alp))/(MA02*Pi2*v2) + 
 (0.0078125*A0i(aa0, MHp2)*(3*(Mh2 - MHH2)*Cos(2*alp) + 
    (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 4*(-2*M2 + Mh2 + MHH2)*Cos(2*beta))*
   (pow(Sec(beta), 3))*Sin(alp))/(MA02*Pi2*v2) + 
 (0.0078125*(MA02 - Mh2)*B0i(bb0, 0, MA02, Mh2)*Cos(alp - beta)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*(pow(Sec(beta), 3))*Sin(alp))/(MA02*Pi2*v2) + 
 (0.0078125*(MA02 - Mh2)*B0i(bb0, MA02, MA02, Mh2)*Cos(alp - beta)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*(pow(Sec(beta), 3))*Sin(alp))/(MA02*Pi2*v2) + 
 (0.0078125*A0i(aa0, Mh2)*Cos(alp - beta)*(-2*(M2 - MA02)*Cos(alp - 3*beta) + 
    (Mh2 - MHH2)*Cos(3*alp - beta) + (-2*M2 - 2*MA02 + 3*Mh2 + MHH2)*
     Cos(alp + beta))*(pow(Sec(beta), 3))*Sin(alp))/(MA02*Pi2*v2) + 
 (0.0625*B0i(bb0, m32, MHH2, ML2)*(pow(ML, 2))*(pow(Cos(alp), 2))*
   (pow(Sec(beta), 3))*Sin(alp))/(Pi2*v2) - 
 (0.03125*(m12 + m22 - m32)*C0i(cc1, m22, m12, m32, MHH2, ML2, ML2)*
   (pow(ML, 2))*(pow(Cos(alp), 2))*(pow(Sec(beta), 3))*Sin(alp))/
  (Pi2*v2) + (0.03125*(m12 - m22 + m32)*C0i(cc2, m22, m12, m32, MHH2, ML2, 
    ML2)*(pow(ML, 2))*(pow(Cos(alp), 2))*(pow(Sec(beta), 3))*
   Sin(alp))/(Pi2*v2) + (0.03125*C0i(cc0, m22, m12, m32, MHH2, ML2, ML2)*
   (pow(ML, 2))*(-m12 - m22 + m32 + 4*(pow(ML, 2)))*
   (pow(Cos(alp), 2))*(pow(Sec(beta), 3))*Sin(alp))/(Pi2*v2) + 
 (0.00390625*(-3*(Mh2 - MHH2)*Cos(2*alp) + (-Mh2 + MHH2)*
     Cos(2*(alp - 2*beta)) + 4*(MA02 + (2*M2 + MA02 - Mh2 - MHH2)*
       Cos(2*beta)))*(pow(Sec(beta), 3))*Re(A0i(aa0, MA02))*Sin(alp))/
  (MA02*Pi2*v2) + (0.00390625*(4*MA02 - 2*(5*M2 - 6*MHH2)*Cos(2*alp) - 
    2*M2*Cos(2*(alp - 2*beta)) + Mh2*Cos(4*alp - 2*beta) - 
    MHH2*Cos(4*alp - 2*beta) + 12*M2*Cos(2*beta) + 4*MA02*Cos(2*beta) - 
    Mh2*Cos(2*beta) - 11*MHH2*Cos(2*beta))*(pow(Sec(beta), 3))*
   Re(A0i(aa0, MHH2))*Sin(alp))/(MA02*Pi2*v2) + 
 (0.0078125*(2*MA02 - 3*(Mh2 - MHH2)*Cos(2*alp) + 
    (-Mh2 + MHH2)*Cos(2*(alp - 2*beta)) + 8*M2*Cos(2*beta) + 
    2*MA02*Cos(2*beta) - 4*Mh2*Cos(2*beta) - 4*MHH2*Cos(2*beta) + 
    MA02*Cos(2*(beta - 2*w)) - 8*MA02*Cos(2*(beta - w)) - 16*MA02*Cos(2*w) + 
    2*MA02*Cos(4*w) - 8*MA02*Cos(2*(beta + w)) + MA02*Cos(2*(beta + 2*w)))*
   (pow(Sec(beta), 3))*Re(A0i(aa0, MHp2))*Sin(alp))/(MA02*Pi2*v2) - 
 (0.25*(pow(ME, 4))*(pow(Cos(alp), 2))*(pow(Sec(beta), 3))*
   Re(B0i(bb0, Mh2, ME2, ME2))*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.25*(pow(ML, 4))*(pow(Cos(alp), 2))*(pow(Sec(beta), 3))*
   Re(B0i(bb0, Mh2, ML2, ML2))*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.25*(pow(MM, 4))*(pow(Cos(alp), 2))*(pow(Sec(beta), 3))*
   Re(B0i(bb0, Mh2, MM2, MM2))*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*(pow(ME, 4))*(pow(Cos(alp), 2))*(pow(Sec(beta), 3))*
   Re(B0i(bb0, MHH2, ME2, ME2))*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*(pow(ML, 4))*(pow(Cos(alp), 2))*(pow(Sec(beta), 3))*
   Re(B0i(bb0, MHH2, ML2, ML2))*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*(pow(MM, 4))*(pow(Cos(alp), 2))*(pow(Sec(beta), 3))*
   Re(B0i(bb0, MHH2, MM2, MM2))*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*(pow(ML, 2))*(pow(Cos(alp), 2))*(pow(Sec(beta), 3))*
   Re(B0i(bb0, ML2, MHH2, ML2))*Sin(alp))/(Pi2*v2) - 
 (0.125*MHH2*(pow(ME, 2))*(pow(Cos(alp), 2))*(pow(Sec(beta), 3))*
   Re(B0i(bb1, MHH2, ME2, ME2))*Sin(alp))/((-Mh2 + MHH2)*Pi2*v2) - 
 (0.125*MHH2*(pow(ML, 2))*(pow(Cos(alp), 2))*(pow(Sec(beta), 3))*
   Re(B0i(bb1, MHH2, ML2, ML2))*Sin(alp))/((-Mh2 + MHH2)*Pi2*v2) - 
 (0.125*MHH2*(pow(MM, 2))*(pow(Cos(alp), 2))*(pow(Sec(beta), 3))*
   Re(B0i(bb1, MHH2, MM2, MM2))*Sin(alp))/((-Mh2 + MHH2)*Pi2*v2) + 
 (0.0009765625*(pow((2*MA02 - Mh2)*Cos(alp - 3*beta) + 
      (4*M2 - 2*MA02 - 3*Mh2)*Cos(alp + beta), 2))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 3))*Re(B0i(dbb0, Mh2, MA02, MA02))*Sin(alp))/
  (Pi2*v2) + (0.03125*Cos(beta)*((pow(MA02 - Mh2, 2)) - 
    (MA02 + Mh2)*MW2*(pow(Sec(w), 2)))*
   (pow(Cos(alp) + Sin(alp)*Tan(beta), 2))*Re(B0i(dbb0, Mh2, MA02, MZ2))*
   Sin(alp))/(Pi2*v2) + (0.0078125*(pow(Cos(alp - beta), 2))*
   (pow(Csc(beta), 2))*(pow(Sec(beta), 3))*
   (pow((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta), 2))*
   Re(B0i(dbb0, Mh2, Mh2, MHH2))*Sin(alp))/(Pi2*v2) + 
 (0.001953125*(pow((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + 
      (-4*M2 + 3*Mh2 + 2*MHp2)*Cos(alp + beta), 2))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 3))*Re(B0i(dbb0, Mh2, MHp2, MHp2))*Sin(alp))/
  (Pi2*v2) + (0.25*(pow(ML, 4))*(pow(Cos(alp), 2))*
   (pow(Sec(beta), 3))*Re(B0i(dbb0, ML2, MHH2, ML2))*Sin(alp))/(Pi2*v2) + 
 (0.125*(pow(ML, 4))*(pow(Cos(alp), 2))*(pow(Sec(beta), 3))*
   Re(B0i(dbb1, ML2, MHH2, ML2))*Sin(alp))/(Pi2*v2) - 
 (0.75*MT2*A0i(aa0, MT2)*Sec(beta)*Sin(alp))/(MA02*Pi2*v2) - 
 (0.375*MT2*B0i(bb1, MA02, MT2, MT2)*Sec(beta)*Sin(alp))/(Pi2*v2) - 
 (0.75*A0i(aa0, MB2)*(pow(MB, 2))*Sec(beta)*Sin(alp))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MB2, MB2)*(pow(MB, 2))*Sec(beta)*Sin(alp))/
  (Pi2*v2) - (0.75*A0i(aa0, MC2)*(pow(MC, 2))*Sec(beta)*Sin(alp))/
  (MA02*Pi2*v2) - (0.375*B0i(bb1, MA02, MC2, MC2)*(pow(MC, 2))*Sec(beta)*
   Sin(alp))/(Pi2*v2) - (0.75*A0i(aa0, MD2)*(pow(MD, 2))*Sec(beta)*
   Sin(alp))/(MA02*Pi2*v2) - (0.375*B0i(bb1, MA02, MD2, MD2)*(pow(MD, 2))*
   Sec(beta)*Sin(alp))/(Pi2*v2) - 
 (0.75*A0i(aa0, MS2)*(pow(MS, 2))*Sec(beta)*Sin(alp))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MS2, MS2)*(pow(MS, 2))*Sec(beta)*Sin(alp))/
  (Pi2*v2) - (0.75*A0i(aa0, MU2)*(pow(MU, 2))*Sec(beta)*Sin(alp))/
  (MA02*Pi2*v2) - (0.375*B0i(bb1, MA02, MU2, MU2)*(pow(MU, 2))*Sec(beta)*
   Sin(alp))/(Pi2*v2) + (0.25*A0i(aa0, ME2)*(pow(ME, 2))*
   (pow(Tan(beta), 2))*Sec(beta)*Sin(alp))/(MA02*Pi2*v2) + 
 (0.125*B0i(bb1, MA02, ME2, ME2)*(pow(ME, 2))*(pow(Tan(beta), 2))*
   Sec(beta)*Sin(alp))/(Pi2*v2) + 
 (0.25*A0i(aa0, ML2)*(pow(ML, 2))*(pow(Tan(beta), 2))*Sec(beta)*
   Sin(alp))/(MA02*Pi2*v2) + (0.125*B0i(bb1, MA02, ML2, ML2)*(pow(ML, 2))*
   (pow(Tan(beta), 2))*Sec(beta)*Sin(alp))/(Pi2*v2) + 
 (0.03125*(m12 + m22 - m32)*C0i(cc1, m22, m12, m32, MA02, ML2, ML2)*
   (pow(ML, 2))*(pow(Tan(beta), 2))*Sec(beta)*Sin(alp))/(Pi2*v2) - 
 (0.03125*(m12 - m22 + m32)*C0i(cc2, m22, m12, m32, MA02, ML2, ML2)*
   (pow(ML, 2))*(pow(Tan(beta), 2))*Sec(beta)*Sin(alp))/(Pi2*v2) + 
 (0.03125*C0i(cc0, m22, m12, m32, MA02, ML2, ML2)*
   (m12 + m22 - m32 - 4*(pow(ML, 2)))*(pow(ML, 2))*
   (pow(Tan(beta), 2))*Sec(beta)*Sin(alp))/(Pi2*v2) + 
 (0.25*A0i(aa0, MM2)*(pow(MM, 2))*(pow(Tan(beta), 2))*Sec(beta)*
   Sin(alp))/(MA02*Pi2*v2) + (0.125*B0i(bb1, MA02, MM2, MM2)*(pow(MM, 2))*
   (pow(Tan(beta), 2))*Sec(beta)*Sin(alp))/(Pi2*v2) + 
 (0.75*(pow(MB, 4))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, Mh2, MB2, MB2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*(pow(MC, 4))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, Mh2, MC2, MC2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*(pow(MD, 4))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, Mh2, MD2, MD2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*(pow(MS, 4))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, Mh2, MS2, MS2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*(pow(MT, 4))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, Mh2, MT2, MT2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*(pow(MU, 4))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, Mh2, MU2, MU2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*(pow(MB, 4))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MB2, MB2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*(pow(MC, 4))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MC2, MC2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*(pow(MD, 4))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MD2, MD2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*(pow(MS, 4))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MS2, MS2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*(pow(MT, 4))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MT2, MT2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*(pow(MU, 4))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb0, MHH2, MU2, MU2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0625*(pow(ML, 2))*(pow(Tan(beta), 2))*
   Re(B0i(bb0, ML2, MA02, ML2))*Sec(beta)*Sin(alp))/(Pi2*v2) - 
 (0.125*MW2*Cos(2*w)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Re(B0i(bb0, MW2, MHH2, MW2))*Sec(beta)*Sin(alp))/(Pi2*v2) + 
 (0.125*MW2*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Re(B0i(bb0, MZ2, MHH2, MZ2))*Sec(beta)*Sin(alp))/(Pi2*v2) + 
 (0.5*Cos(2*w)*(pow(Sin(w), 2))*Re(B0i(bb00, 0, MHp2, MHp2))*Sec(beta)*
   Sin(alp))/(Pi2*v2) + (0.125*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, MA02, MHp2))*Sec(beta)*Sin(alp))/(Pi2*v2) + 
 (0.125*Cos(2*w)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MW2, Mh2, MHp2))*Sec(beta)*Sin(alp))/(Pi2*v2) + 
 (0.125*Cos(2*w)*(pow(Csc(w), 2))*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb00, MW2, MHH2, MHp2))*Sec(beta)*Sin(alp))/(Pi2*v2) + 
 (0.125*Cos(2*w)*(pow(Cos(alp - beta), 2))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MW2, MHH2, MW2))*Sec(beta)*Sin(alp))/(Pi2*v2) - 
 (0.125*(pow(Cos(alp - beta), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MA02, Mh2))*Sec(beta)*Sin(alp))/(Pi2*v2) - 
 (0.125*(pow(Cot(w), 2))*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb00, MZ2, MA02, MHH2))*Sec(beta)*Sin(alp))/(Pi2*v2) - 
 (0.125*(pow(Cos(alp - beta), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MHH2, MZ2))*Sec(beta)*Sin(alp))/(Pi2*v2) - 
 (0.125*(pow(Cos(2*w), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MHp2, MHp2))*Sec(beta)*Sin(alp))/(Pi2*v2) - 
 (0.375*MHH2*(pow(MB, 2))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MB2, MB2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*(pow(MC, 2))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MC2, MC2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*(pow(MD, 2))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MD2, MD2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*(pow(MS, 2))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MS2, MS2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*MT2*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MT2, MT2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*(pow(MU, 2))*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
   Re(B0i(bb1, MHH2, MU2, MU2))*Sec(beta)*Sin(alp))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0625*(MHp2*(MHp2 - MW2) - Mh2*(2*MHp2 + MW2) + (pow(Mh, 4)))*
   (pow(Cos(alp - beta), 2))*Re(B0i(dbb0, Mh2, MHp2, MW2))*Sec(beta)*
   Sin(alp))/(Pi2*v2) + (0.5*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MHp2, MHp2))*Sec(beta)*Sin(alp))/(Pi2*v2) + 
 (0.0625*Mh2*MW2*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb1, Mh2, MA02, MZ2))*Sec(beta)*Sin(alp))/(Pi2*v2) + 
 (0.125*Mh2*MW2*(pow(Cos(alp - beta), 2))*Re(B0i(dbb1, Mh2, MHp2, MW2))*
   Sec(beta)*Sin(alp))/(Pi2*v2) - 
 (0.125*(pow(ML, 4))*(pow(Tan(beta), 2))*Re(B0i(dbb1, ML2, 0, MHp2))*
   Sec(beta)*Sin(alp))/(Pi2*v2) + 
 (0.125*(pow(ML, 4))*(pow(Tan(beta), 2))*
   Re(B0i(dbb1, ML2, MA02, ML2))*Sec(beta)*Sin(alp))/(Pi2*v2) - 
 1.*Sec(beta)*(Cos(beta) + Sin(alp)) - 
 (0.03125*m12*C0i(cc1, m12, m32, m22, ML2, ML2, MZ2)*
   (4*MW2 + 4*MW2*Cos(4*w) + (pow(ML, 2)) + 
    Cos(2*w)*(-8*MW2 + (pow(ML, 2))))*(pow(Sec(w), 2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.015625*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, ML2, ML2, MZ2)*
   (4*MW2 + 4*MW2*Cos(4*w) + (pow(ML, 2)) + 
    Cos(2*w)*(-8*MW2 + (pow(ML, 2))))*(pow(Sec(w), 2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*C0i(cc0, m12, m32, m22, ML2, ML2, MZ2)*(pow(ML, 2))*
   (4*MW2 + 4*MW2*Cos(4*w) + (pow(ML, 2)) + 
    Cos(2*w)*(-8*MW2 + (pow(ML, 2))))*(pow(Sec(w), 2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (1.*MW2*B0i(bb0, m32, 0, ML2)*(pow(Sin(w), 2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.5*(m12 + m22 - m32)*MW2*C0i(cc1, m22, m12, m32, 0, ML2, ML2)*
   (pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.5*(m12 - m22 + m32)*MW2*C0i(cc2, m22, m12, m32, 0, ML2, ML2)*
   (pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.5*MW2*C0i(cc0, m22, m12, m32, 0, ML2, ML2)*
   (m12 + m22 - m32 - 4*(pow(ML, 2)))*(pow(Sin(w), 2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (1.125*MW2*(pow(Sin(w), 2))*Re(B0i(bb0, 0, MW2, MW2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(bb0, ML2, 0, ML2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.03125*(4*MW2 + 4*MW2*Cos(4*w) + (pow(ML, 2)) + 
    Cos(2*w)*(-8*MW2 + (pow(ML, 2))))*(pow(Sec(w), 2))*
   Re(B0i(bb0, ML2, ML2, MZ2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.5*MW2*Cos(2*w)*Re(B0i(bb0, MW2, 0, MW2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.015625*MW2*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Sec(w), 2))*Re(B0i(bb0, MW2, MW2, MZ2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.1875*(pow(MB, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MB2, MB2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.1875*(pow(MC, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MC2, MC2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.1875*(pow(MD, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MD2, MD2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*(pow(ME, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, ME2, ME2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*(pow(ML, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, ML2, ML2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*(pow(MM, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MM2, MM2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.1875*(pow(MS, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MS2, MS2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.1875*MT2*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MT2, MT2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.1875*(pow(MU, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MU2, MU2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*MW2*(5 + 9*Cos(2*w))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MW2, MW2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MB2, MB2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MC2, MC2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MD2, MD2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.5*(-1 + 2*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, ME2, ME2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.5*(-1 + 2*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, ML2, ML2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.5*(-1 + 2*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, MM2, MM2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MS2, MS2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MT2, MT2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MU2, MU2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.5*(2 + 3*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, MW2, MW2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (1.*Cos(2*w)*Re(B0i(bb00, MW2, 0, MW2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*v2) + (0.125*(2 + 5*Cos(2*w) + 2*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MW2, MW2, MZ2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.375*(pow(Cot(w), 2))*Re(B0i(bb00, MZ2, 0, 0))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MB2, MB2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MC2, MC2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MD2, MD2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ME2, ME2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ML2, ML2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MM2, MM2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MS2, MS2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MT2, MT2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MU2, MU2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*(7 + 8*Cos(2*w) + 3*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MW2, MW2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MB2, MB2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MC2, MC2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MD2, MD2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ME2, ME2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ML2, ML2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MM2, MM2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MS2, MS2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MT2, MT2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MU2, MU2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.25*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MW2, MW2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (0.25*MW2*Cos(2*w)*Re(B0i(bb1, MW2, 0, MW2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.25*MW2*Cos(2*w)*(pow(Cot(w), 2))*Re(B0i(bb1, MW2, MW2, MZ2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.1875*MW2*(pow(Csc(w), 2))*Re(B0i(bb1, MZ2, 0, 0))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MB2, MB2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MC2, MC2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MD2, MD2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, ME2, ME2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, ML2, ML2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MM2, MM2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MS2, MS2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MT2, MT2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MU2, MU2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.25*MW2*(pow(Cos(w), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb1, MZ2, MW2, MW2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (1.*MW2*(pow(ML, 2))*(pow(Sin(w), 2))*Re(B0i(dbb0, ML2, 0, ML2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*(pow(ML, 2))*(4*MW2 + 4*MW2*Cos(4*w) + (pow(ML, 2)) + 
    Cos(2*w)*(-8*MW2 + (pow(ML, 2))))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, ML2, ML2, MZ2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MB2, MB2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MC2, MC2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MD2, MD2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ME2, ME2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ML2, ML2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MM2, MM2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MS2, MS2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MT2, MT2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MU2, MU2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (1.5*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MW2, MW2))*Sec(beta)*
   (Cos(beta) + Sin(alp)))/(Pi2*v2) + 
 (1.*MW2*(pow(ML, 2))*(pow(Sin(w), 2))*Re(B0i(dbb1, ML2, 0, ML2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.125*(pow(ML, 2))*(2*MW2 + (pow(ML, 2)))*
   Re(B0i(dbb1, ML2, 0, MW2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*(pow(ML, 2))*(4*MW2 + 2*MW2*Cos(4*w) + (pow(ML, 2)) + 
    Cos(2*w)*(-4*MW2 + (pow(ML, 2))))*(pow(Sec(w), 2))*
   Re(B0i(dbb1, ML2, ML2, MZ2))*Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*v2) - 
 (0.0625*(pow(ME, 2))*Re(B0i(bb1, Mh2, ME2, ME2))*
   (2*(Mh2 - MHH2) + (2*Mh2 - MHH2 + MHH2*Cos(2*alp))*(pow(Sec(beta), 3))*
     Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*(pow(ML, 2))*Re(B0i(bb1, Mh2, ML2, ML2))*
   (2*(Mh2 - MHH2) + (2*Mh2 - MHH2 + MHH2*Cos(2*alp))*(pow(Sec(beta), 3))*
     Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*(pow(MM, 2))*Re(B0i(bb1, Mh2, MM2, MM2))*
   (2*(Mh2 - MHH2) + (2*Mh2 - MHH2 + MHH2*Cos(2*alp))*(pow(Sec(beta), 3))*
     Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0087890625*Re(B0i(dbb0, Mh2, Mh2, Mh2))*(16*(pow(Mh, 4)) + 
    (pow(M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
        (2*M2 - 3*Mh2)*Cos(alp + beta), 2))*(pow(Csc(beta), 2))*
     (pow(Sec(beta), 3))*Sin(alp)))/(Pi2*v2) + 
 (0.375*(MT2 - MW2)*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb0, MW2, MB2, MT2))*
   (1 + Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.375*(pow(MC, 2))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MW2, MC2, MS2))*(1 + Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.375*(-MW2 + (pow(MU, 2)))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MW2, MD2, MU2))*(1 + Sec(beta)*Sin(alp)))/(Pi2*v2) - 
 (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, ME2))*
   (1 + Sec(beta)*Sin(alp)))/(Pi2*v2) - 
 (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, ML2))*
   (1 + Sec(beta)*Sin(alp)))/(Pi2*v2) - 
 (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, MM2))*
   (1 + Sec(beta)*Sin(alp)))/(Pi2*v2) - 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MB2, MT2))*
   (1 + Sec(beta)*Sin(alp)))/(Pi2*v2) - 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MC2, MS2))*
   (1 + Sec(beta)*Sin(alp)))/(Pi2*v2) - 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MD2, MU2))*
   (1 + Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ME2))*
   (1 + Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ML2))*
   (1 + Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, MM2))*
   (1 + Sec(beta)*Sin(alp)))/(Pi2*v2) - 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MB2, MT2))*
   (1 + Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MC2, MS2))*
   (1 + Sec(beta)*Sin(alp)))/(Pi2*v2) - 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MD2, MU2))*
   (1 + Sec(beta)*Sin(alp)))/(Pi2*v2) - 
 (0.375*(pow(MB, 2))*Re(B0i(bb1, Mh2, MB2, MB2))*
   (Mh2 - MHH2 - MHH2*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
     Sec(beta)*Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*(pow(MC, 2))*Re(B0i(bb1, Mh2, MC2, MC2))*
   (Mh2 - MHH2 - MHH2*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
     Sec(beta)*Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*(pow(MD, 2))*Re(B0i(bb1, Mh2, MD2, MD2))*
   (Mh2 - MHH2 - MHH2*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
     Sec(beta)*Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*(pow(MS, 2))*Re(B0i(bb1, Mh2, MS2, MS2))*
   (Mh2 - MHH2 - MHH2*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
     Sec(beta)*Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MT2*Re(B0i(bb1, Mh2, MT2, MT2))*(Mh2 - MHH2 - 
    MHH2*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*Sec(beta)*Sin(alp)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.375*(pow(MU, 2))*Re(B0i(bb1, Mh2, MU2, MU2))*
   (Mh2 - MHH2 - MHH2*(pow(Cos(alp), 2))*(pow(Csc(beta), 2))*
     Sec(beta)*Sin(alp)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.125*MW2*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb0, MW2, Mh2, MW2))*
   (1 + (pow(Sin(alp - beta), 2))*Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.125*MW2*(pow(Csc(w), 2))*Re(B0i(bb0, MZ2, Mh2, MZ2))*
   (1 + (pow(Sin(alp - beta), 2))*Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.125*Cos(2*w)*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, Mh2, MW2))*
   (1 + (pow(Sin(alp - beta), 2))*Sec(beta)*Sin(alp)))/(Pi2*v2) - 
 (0.125*(pow(Cot(w), 2))*Re(B0i(bb00, MZ2, Mh2, MZ2))*
   (1 + (pow(Sin(alp - beta), 2))*Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.03125*(-2*Mh2*MW2 + (pow(Mh, 4)) + 12*(pow(MW, 4)))*
   Re(B0i(dbb0, Mh2, MW2, MW2))*(1 + (pow(Sin(alp - beta), 2))*Sec(beta)*
     Sin(alp)))/(Pi2*v2) + 
 (0.001953125*(-8*Mh2*MW2 + 4*Mh2*(Mh2 - 2*MW2)*Cos(2*w) + 
    3*(pow(Mh, 4)) + Cos(4*w)*(pow(Mh, 4)) + 96*(pow(MW, 4)))*
   (pow(Sec(w), 4))*Re(B0i(dbb0, Mh2, MZ2, MZ2))*
   (1 + (pow(Sin(alp - beta), 2))*Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.125*Mh2*MW2*Re(B0i(dbb1, Mh2, MW2, MW2))*
   (1 + (pow(Sin(alp - beta), 2))*Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.0625*Mh2*MW2*(pow(Sec(w), 2))*Re(B0i(dbb1, Mh2, MZ2, MZ2))*
   (1 + (pow(Sin(alp - beta), 2))*Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.00390625*Re(A0i(aa0, Mh2))*(8*MA02 + 
    (4*MA02 + 2*(5*M2 - 6*Mh2)*Cos(2*alp) + 2*M2*Cos(2*(alp - 2*beta)) - 
      Mh2*Cos(4*alp - 2*beta) + MHH2*Cos(4*alp - 2*beta) + 
      12*M2*Cos(2*beta) + 4*MA02*Cos(2*beta) - 11*Mh2*Cos(2*beta))*
     (pow(Sec(beta), 3))*Sin(alp) + MHH2*(-1 + (pow(Tan(beta), 2)))*
     Sec(beta)*Sin(alp)))/(MA02*Pi2*v2) + 
 (0.015625*Re(A0i(aa0, ME2))*((-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
     (pow(Csc(w), 2)) - 4*(pow(MA02, -1))*
     (2*MA02 - 2*MA02*Cos(2*w)*(pow(Cot(w), 2)) + 
      MA02*Cos(4*w)*(pow(Cot(w), 2)) - 4*MA02*(pow(Sin(w), 2)) + 
      16*MA02*(pow(Sin(w), 4)) + 4*(pow(ME, 2))*
       (pow(Cos(alp), 2))*(pow(Tan(beta), 2)) + 
      4*(pow(ME, 2))*(pow(Sin(alp), 2))*(pow(Tan(beta), 2)))*
     Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.015625*Re(A0i(aa0, ML2))*((-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
     (pow(Csc(w), 2)) - 4*(pow(MA02, -1))*
     (2*MA02 - 2*MA02*Cos(2*w)*(pow(Cot(w), 2)) + 
      MA02*Cos(4*w)*(pow(Cot(w), 2)) - 4*MA02*(pow(Sin(w), 2)) + 
      16*MA02*(pow(Sin(w), 4)) + 4*(pow(ML, 2))*
       (pow(Cos(alp), 2))*(pow(Tan(beta), 2)) + 
      4*(pow(ML, 2))*(pow(Sin(alp), 2))*(pow(Tan(beta), 2)))*
     Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.015625*Re(A0i(aa0, MM2))*((-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
     (pow(Csc(w), 2)) - 4*(pow(MA02, -1))*
     (2*MA02 - 2*MA02*Cos(2*w)*(pow(Cot(w), 2)) + 
      MA02*Cos(4*w)*(pow(Cot(w), 2)) - 4*MA02*(pow(Sin(w), 2)) + 
      16*MA02*(pow(Sin(w), 4)) + 4*(pow(MM, 2))*
       (pow(Cos(alp), 2))*(pow(Tan(beta), 2)) + 
      4*(pow(MM, 2))*(pow(Sin(alp), 2))*(pow(Tan(beta), 2)))*
     Sec(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.010416666666666666*(pow(Csc(w), 2))*Re(A0i(aa0, MT2))*
   (-29*MA02 + 2*MA02*Cos(6*w) - 29*MA02*Sec(beta)*Sin(alp) + 
    36*MT2*Sec(beta)*Sin(alp) + 2*MA02*Cos(6*w)*Sec(beta)*Sin(alp) - 
    20*MA02*Cos(4*w)*Sec(beta)*(Cos(beta) + Sin(alp)) + 
    Cos(2*w)*Sec(beta)*(29*MA02*Cos(beta) + (29*MA02 - 36*MT2)*Sin(alp))))/
  (MA02*Pi2*v2) + (0.005208333333333333*(pow(Csc(w), 2))*
   Re(A0i(aa0, MB2))*(-16*MA02 + MA02*Cos(6*w) - 16*MA02*Sec(beta)*Sin(alp) + 
    MA02*Cos(6*w)*Sec(beta)*Sin(alp) + 72*(pow(MB, 2))*Sec(beta)*
     Sin(alp) - 10*MA02*Cos(4*w)*Sec(beta)*(Cos(beta) + Sin(alp)) + 
    Cos(2*w)*Sec(beta)*(61*MA02*Cos(beta) + (61*MA02 - 72*(pow(MB, 2)))*
       Sin(alp))))/(MA02*Pi2*v2) + 
 (0.010416666666666666*(pow(Csc(w), 2))*Re(A0i(aa0, MC2))*
   (-29*MA02 + 2*MA02*Cos(6*w) - 29*MA02*Sec(beta)*Sin(alp) + 
    2*MA02*Cos(6*w)*Sec(beta)*Sin(alp) + 36*(pow(MC, 2))*Sec(beta)*
     Sin(alp) - 20*MA02*Cos(4*w)*Sec(beta)*(Cos(beta) + Sin(alp)) + 
    Cos(2*w)*Sec(beta)*(29*MA02*Cos(beta) + (29*MA02 - 36*(pow(MC, 2)))*
       Sin(alp))))/(MA02*Pi2*v2) + 
 (0.005208333333333333*(pow(Csc(w), 2))*Re(A0i(aa0, MD2))*
   (-16*MA02 + MA02*Cos(6*w) - 16*MA02*Sec(beta)*Sin(alp) + 
    MA02*Cos(6*w)*Sec(beta)*Sin(alp) + 72*(pow(MD, 2))*Sec(beta)*
     Sin(alp) - 10*MA02*Cos(4*w)*Sec(beta)*(Cos(beta) + Sin(alp)) + 
    Cos(2*w)*Sec(beta)*(61*MA02*Cos(beta) + (61*MA02 - 72*(pow(MD, 2)))*
       Sin(alp))))/(MA02*Pi2*v2) + 
 (0.005208333333333333*(pow(Csc(w), 2))*Re(A0i(aa0, MS2))*
   (-16*MA02 + MA02*Cos(6*w) - 16*MA02*Sec(beta)*Sin(alp) + 
    MA02*Cos(6*w)*Sec(beta)*Sin(alp) + 72*(pow(MS, 2))*Sec(beta)*
     Sin(alp) - 10*MA02*Cos(4*w)*Sec(beta)*(Cos(beta) + Sin(alp)) + 
    Cos(2*w)*Sec(beta)*(61*MA02*Cos(beta) + (61*MA02 - 72*(pow(MS, 2)))*
       Sin(alp))))/(MA02*Pi2*v2) + 
 (0.010416666666666666*(pow(Csc(w), 2))*Re(A0i(aa0, MU2))*
   (-29*MA02 + 2*MA02*Cos(6*w) - 29*MA02*Sec(beta)*Sin(alp) + 
    2*MA02*Cos(6*w)*Sec(beta)*Sin(alp) + 36*(pow(MU, 2))*Sec(beta)*
     Sin(alp) - 20*MA02*Cos(4*w)*Sec(beta)*(Cos(beta) + Sin(alp)) + 
    Cos(2*w)*Sec(beta)*(29*MA02*Cos(beta) + (29*MA02 - 36*(pow(MU, 2)))*
       Sin(alp))))/(MA02*Pi2*v2) - 
 (0.1875*(pow(MB, 4))*(pow(Csc(beta), 2))*
   Re(B0i(dbb0, Mh2, MB2, MB2))*Sec(beta)*(Cos(beta) - Cos(3*beta) + 
    Sin(alp) + Sin(3*alp)))/(Pi2*v2) - 
 (0.1875*(pow(MC, 4))*(pow(Csc(beta), 2))*
   Re(B0i(dbb0, Mh2, MC2, MC2))*Sec(beta)*(Cos(beta) - Cos(3*beta) + 
    Sin(alp) + Sin(3*alp)))/(Pi2*v2) - 
 (0.1875*(pow(MD, 4))*(pow(Csc(beta), 2))*
   Re(B0i(dbb0, Mh2, MD2, MD2))*Sec(beta)*(Cos(beta) - Cos(3*beta) + 
    Sin(alp) + Sin(3*alp)))/(Pi2*v2) - 
 (0.1875*(pow(MS, 4))*(pow(Csc(beta), 2))*
   Re(B0i(dbb0, Mh2, MS2, MS2))*Sec(beta)*(Cos(beta) - Cos(3*beta) + 
    Sin(alp) + Sin(3*alp)))/(Pi2*v2) - 
 (0.1875*(pow(MT, 4))*(pow(Csc(beta), 2))*
   Re(B0i(dbb0, Mh2, MT2, MT2))*Sec(beta)*(Cos(beta) - Cos(3*beta) + 
    Sin(alp) + Sin(3*alp)))/(Pi2*v2) - 
 (0.1875*(pow(MU, 4))*(pow(Csc(beta), 2))*
   Re(B0i(dbb0, Mh2, MU2, MU2))*Sec(beta)*(Cos(beta) - Cos(3*beta) + 
    Sin(alp) + Sin(3*alp)))/(Pi2*v2) - 
 (0.09375*Mh2*(pow(MB, 2))*(pow(Csc(beta), 2))*
   Re(B0i(dbb1, Mh2, MB2, MB2))*Sec(beta)*(Cos(beta) - Cos(3*beta) + 
    Sin(alp) + Sin(3*alp)))/(Pi2*v2) - 
 (0.09375*Mh2*(pow(MC, 2))*(pow(Csc(beta), 2))*
   Re(B0i(dbb1, Mh2, MC2, MC2))*Sec(beta)*(Cos(beta) - Cos(3*beta) + 
    Sin(alp) + Sin(3*alp)))/(Pi2*v2) - 
 (0.09375*Mh2*(pow(MD, 2))*(pow(Csc(beta), 2))*
   Re(B0i(dbb1, Mh2, MD2, MD2))*Sec(beta)*(Cos(beta) - Cos(3*beta) + 
    Sin(alp) + Sin(3*alp)))/(Pi2*v2) - 
 (0.09375*Mh2*(pow(MS, 2))*(pow(Csc(beta), 2))*
   Re(B0i(dbb1, Mh2, MS2, MS2))*Sec(beta)*(Cos(beta) - Cos(3*beta) + 
    Sin(alp) + Sin(3*alp)))/(Pi2*v2) - 
 (0.09375*Mh2*MT2*(pow(Csc(beta), 2))*Re(B0i(dbb1, Mh2, MT2, MT2))*
   Sec(beta)*(Cos(beta) - Cos(3*beta) + Sin(alp) + Sin(3*alp)))/(Pi2*v2) - 
 (0.09375*Mh2*(pow(MU, 2))*(pow(Csc(beta), 2))*
   Re(B0i(dbb1, Mh2, MU2, MU2))*Sec(beta)*(Cos(beta) - Cos(3*beta) + 
    Sin(alp) + Sin(3*alp)))/(Pi2*v2) - 
 (0.03125*MW2*B0i(bb0, m32, 0, MW2)*
   (4 + Sec(beta)*(3*Sin(alp) + Sin(alp - 2*beta))))/(Pi2*v2) + 
 (0.015625*B0i(bb0, m32, ML2, MZ2)*
   (-4*((pow(ML, 2)) + MW2*(pow(1 - 2*Cos(2*w), 2))*
       (pow(Sec(w), 2))) - (pow(Sec(w), 2))*Sec(beta)*
     ((11*MW2 + 8*MW2*Cos(4*w) + 2*(pow(ML, 2)) + 
        2*Cos(2*w)*(-8*MW2 + (pow(ML, 2))))*Sin(alp) + 
      MW2*Sin(alp - 2*beta))))/(Pi2*v2) + 
 (0.125*(m22 - MW2)*MW2*C0i(cc0, m22, m12, m32, 0, MW2, MW2)*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.03125*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, 0, MW2, MW2)*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.015625*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, ML2, MZ2, MZ2)*
   (pow(Sec(w), 2))*(1 + Sin(alp - beta)))/(Pi2*v2) + 
 (0.015625*Re(B0i(bb0, Mh2, MW2, MW2))*(4*(-Mh2 + MHH2)*MW2 - 
    ((Mh2*MHH2 - 2*MHH2*MW2 + 12*(pow(MW, 4)))*Sec(beta)*
      (Sin(alp) - Sin(3*alp - 2*beta)))/2 + 
    (Mh2*(MHH2 - 4*MW2) + 2*MW2*(MHH2 + 6*MW2))*Sin(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.0625*MW2*Re(B0i(bb1, Mh2, MW2, MW2))*
   (2*(Mh2 - MHH2) + (2*Mh2 - MHH2 + MHH2*Cos(2*alp - beta)*Sec(beta))*
     Sin(alp - beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.03125*MW2*Re(B0i(bb1, Mh2, MZ2, MZ2))*
   (MHH2*Cos(2*alp - beta)*(pow(Csc(w), 2))*(pow(Sec(w), 2))*
     Sec(beta)*Sin(alp - beta) + (pow(Csc(2*w), 2))*
     (8*(Mh2 - MHH2) + (8*Mh2 - 4*MHH2)*Sin(alp - beta))))/
  ((Mh2 - MHH2)*Pi2*v2*(pow(Csc(w), 2))) + 
 (0.015625*Cos(alp)*((MA02 - Mh2)*(MA02 - MHH2) - 
    (MA02 + MHH2)*MW2*(pow(Sec(w), 2)))*Re(B0i(bb0, MHH2, MA02, MZ2))*
   Sec(beta)*Sin(2*(alp - beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.03125*(Mh2*(MHH2 - MHp2) + MHp2*(MHp2 - MW2) - MHH2*(MHp2 + MW2))*
   Cos(alp)*Re(B0i(bb0, MHH2, MHp2, MW2))*Sec(beta)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - 
 (0.015625*Cos(alp)*(Mh2*MHH2 - 2*MHH2*MW2 + 12*(pow(MW, 4)))*
   Re(B0i(bb0, MHH2, MW2, MW2))*Sec(beta)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - 
 (0.0009765625*Cos(alp)*(3*Mh2*MHH2 - 8*MHH2*MW2 + 
    4*MHH2*(Mh2 - 2*MW2)*Cos(2*w) + Mh2*MHH2*Cos(4*w) + 96*(pow(MW, 4)))*
   (pow(Sec(w), 4))*Re(B0i(bb0, MHH2, MZ2, MZ2))*Sec(beta)*
   Sin(2*(alp - beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.03125*MHH2*MW2*Cos(alp)*(pow(Sec(w), 2))*
   Re(B0i(bb1, MHH2, MA02, MZ2))*Sec(beta)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.0625*MHH2*MW2*Cos(alp)*
   Re(B0i(bb1, MHH2, MHp2, MW2))*Sec(beta)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.0625*MHH2*MW2*Cos(alp)*
   Re(B0i(bb1, MHH2, MW2, MW2))*Sec(beta)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.03125*MHH2*MW2*Cos(alp)*(pow(Sec(w), 2))*
   Re(B0i(bb1, MHH2, MZ2, MZ2))*Sec(beta)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.03125*MW2*Re(B0i(bb0, Mh2, MHp2, MW2))*Sec(beta)*
   (-2*(pow(Cos(alp - beta), 2))*Sin(alp) - 
    (Mh2*(MHH2 - MHp2 - MW2) + MHp2*(-MHH2 + MHp2 - MW2))*Cos(alp)*
     (pow(Mh2 - MHH2, -1))*(pow(MW, -2))*Sin(2*(alp - beta))))/
  (Pi2*v2) + (0.015625*MW2*Re(B0i(bb0, Mh2, MA02, MZ2))*Sec(beta)*
   (-2*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 2))*Sin(alp) - 
    Cos(alp)*(pow(Mh2 - MHH2, -1))*(pow(MW, -2))*
     ((MA02 - Mh2)*(MA02 - MHH2) - (MA02 + Mh2)*MW2*(pow(Sec(w), 2)))*
     Sin(2*(alp - beta))))/(Pi2*v2) + 
 (0.015625*MW2*Re(B0i(bb0, Mh2, MZ2, MZ2))*(-8*(pow(Csc(2*w), 2)) + 
    (pow(Csc(w), 2))*(pow(Sec(w), 2))*Sec(beta)*
     (-2*(pow(Sin(alp - beta), 2))*Sin(alp) + 
      (Cos(alp)*(pow(Mh2 - MHH2, -1))*(pow(MW, -2))*
        (3*Mh2*MHH2 - 8*Mh2*MW2 + 4*Mh2*(MHH2 - 2*MW2)*Cos(2*w) + 
         Mh2*MHH2*Cos(4*w) + 96*(pow(MW, 4)))*(pow(Sec(w), 2))*
        Sin(2*(alp - beta)))/16)))/(Pi2*v2*(pow(Csc(w), 2))) + 
 (0.03125*(pow(Csc(2*beta), 3))*(pow(Sin(alp - beta), 2))*
   (pow((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta), 2))*
   Re(B0i(dbb0, Mh2, MHH2, MHH2))*Sin(alp)*Sin(beta))/(Pi2*v2) - 
 (0.03125*MW2*Cos(alp - beta)*(pow(Sec(w), 2))*
   Re(B0i(bb1, Mh2, MA02, MZ2))*Sec(beta)*(MHH2*Sin(2*alp - beta) + 
    (-2*Mh2 + MHH2)*Sin(beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*MW2*Cos(alp - beta)*Re(B0i(bb1, Mh2, MHp2, MW2))*Sec(beta)*
   (MHH2*Sin(2*alp - beta) + (-2*Mh2 + MHH2)*Sin(beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.03125*C0i(cc0, m12, m32, m22, MHH2, MHH2, ML2)*
   Csc(beta)*(pow(ML, 2))*(pow(Cos(alp), 2))*(pow(Sec(beta), 3))*
   Sin(alp - beta)*((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (Pi2*v2) + (0.005859375*Cos(alp)*Cos(alp - beta)*
   (M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 3))*Re(B0i(bb0, Mh2, Mh2, Mh2))*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.005859375*Cos(alp)*Cos(alp - beta)*
   (M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
    (2*M2 - 3*Mh2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 3))*Re(B0i(bb0, MHH2, Mh2, Mh2))*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.015625*C0i(cc0, m12, m22, m32, Mh2, MHH2, ML2)*
   Cos(alp - beta)*Csc(beta)*(pow(ML, 2))*(pow(Sec(beta), 3))*
   Sin(2*alp)*((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (Pi2*v2) + (0.015625*C0i(cc0, m12, m32, m22, Mh2, MHH2, ML2)*
   Cos(alp - beta)*Csc(beta)*(pow(ML, 2))*(pow(Sec(beta), 3))*
   Sin(2*alp)*((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/
  (Pi2*v2) - (0.001953125*Cos(alp)*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 3))*Re(B0i(bb0, Mh2, Mh2, MHH2))*Sin(2*(alp - beta))*
   (-9*M2*Mh2 - 9*M2*MHH2 + 5*Mh2*MHH2 + 8*(pow(M2, 2)) + 
    2*(pow(Mh, 4)) + 2*(pow(MHH2, 2)) - 
    Cos(4*alp)*(5*Mh2*MHH2 - 9*M2*(Mh2 + MHH2) + 9*(pow(M2, 2)) + 
      2*(pow(Mh, 4)) + 2*(pow(MHH2, 2))) + 
    (pow(M2, 2))*(pow(Cos(beta), 4)) - 6*(pow(M2, 2))*
     (pow(Cos(beta), 2))*(pow(Sin(beta), 2)) + 
    (pow(M2, 2))*(pow(Sin(beta), 4)) - 2*M2*(Mh2 - MHH2)*Sin(2*alp)*
     Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.001953125*Cos(alp)*(pow(Csc(beta), 2))*(pow(Sec(beta), 3))*
   Re(B0i(bb0, MHH2, Mh2, MHH2))*Sin(2*(alp - beta))*
   (-9*M2*Mh2 - 9*M2*MHH2 + 5*Mh2*MHH2 + 8*(pow(M2, 2)) + 
    2*(pow(Mh, 4)) + 2*(pow(MHH2, 2)) - 
    Cos(4*alp)*(5*Mh2*MHH2 - 9*M2*(Mh2 + MHH2) + 9*(pow(M2, 2)) + 
      2*(pow(Mh, 4)) + 2*(pow(MHH2, 2))) + 
    (pow(M2, 2))*(pow(Cos(beta), 4)) - 6*(pow(M2, 2))*
     (pow(Cos(beta), 2))*(pow(Sin(beta), 2)) + 
    (pow(M2, 2))*(pow(Sin(beta), 4)) - 2*M2*(Mh2 - MHH2)*Sin(2*alp)*
     Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.005859375*Cos(alp)*(pow(Csc(beta), 2))*(pow(Sec(beta), 3))*
   Re(B0i(bb0, Mh2, MHH2, MHH2))*Sin(alp - beta)*
   ((-3*M2 + Mh2 + 2*MHH2)*Sin(2*alp) - M2*Sin(2*beta))*
   (M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.005859375*Cos(alp)*(pow(Csc(beta), 2))*(pow(Sec(beta), 3))*
   Re(B0i(bb0, MHH2, MHH2, MHH2))*Sin(alp - beta)*
   ((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta))*
   (M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0009765625*Cos(alp)*((-2*MA02 + Mh2)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 3))*Re(B0i(bb0, Mh2, MA02, MA02))*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0009765625*Cos(alp)*((-2*MA02 + Mh2)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 3))*Re(B0i(bb0, MHH2, MA02, MA02))*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0078125*(MA02 - MHH2)*B0i(bb0, 0, MA02, MHH2)*(pow(Sec(beta), 3))*
   Sin(alp)*Sin(alp - beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/(MA02*Pi2*v2) - 
 (0.0078125*(MA02 - MHH2)*B0i(bb0, MA02, MA02, MHH2)*(pow(Sec(beta), 3))*
   Sin(alp)*Sin(alp - beta)*((-2*MA02 + MHH2)*Sin(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*MHH2)*Sin(alp + beta)))/(MA02*Pi2*v2) + 
 (0.0078125*A0i(aa0, MHH2)*(pow(Sec(beta), 3))*Sin(alp)*Sin(alp - beta)*
   (-2*(M2 - MA02)*Sin(alp - 3*beta) + (Mh2 - MHH2)*Sin(3*alp - beta) + 
    (-2*M2 - 2*MA02 + Mh2 + 3*MHH2)*Sin(alp + beta)))/(MA02*Pi2*v2) - 
 (0.001953125*Cos(alp)*((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh2 + 2*MHp2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 3))*Re(B0i(bb0, Mh2, MHp2, MHp2))*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.001953125*Cos(alp)*((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + 
    (-4*M2 + 3*Mh2 + 2*MHp2)*Cos(alp + beta))*(pow(Csc(beta), 2))*
   (pow(Sec(beta), 3))*Re(B0i(bb0, MHH2, MHp2, MHp2))*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0625*MW2*B0i(bb0, m32, 0, MHp2)*Cos(alp - beta)*Tan(beta))/(Pi2*v2) + 
 (0.0625*(m22 - MHp2)*MW2*C0i(cc0, m22, m12, m32, 0, MHp2, MW2)*
   Cos(alp - beta)*Tan(beta))/(Pi2*v2) + 
 (0.0625*(m22 - MW2)*MW2*C0i(cc0, m32, m12, m22, 0, MHp2, MW2)*
   Cos(alp - beta)*Tan(beta))/(Pi2*v2) - 
 (0.03125*(m12 - m22 - m32)*MW2*C0i(cc1, m22, m12, m32, 0, MHp2, MW2)*
   Cos(alp - beta)*Tan(beta))/(Pi2*v2) - 
 (0.03125*(3*m12 - 3*m22 + m32)*MW2*C0i(cc1, m32, m12, m22, 0, MHp2, MW2)*
   Cos(alp - beta)*Tan(beta))/(Pi2*v2) + 
 (0.0625*m32*MW2*C0i(cc2, m22, m12, m32, 0, MHp2, MW2)*Cos(alp - beta)*
   Tan(beta))/(Pi2*v2) + (0.0625*(m12 + 2*m22 - m32)*MW2*
   C0i(cc2, m32, m12, m22, 0, MHp2, MW2)*Cos(alp - beta)*Tan(beta))/
  (Pi2*v2) + (0.015625*C0i(cc0, m12, m32, m22, MA02, MA02, ML2)*
   ((-2*MA02 + Mh2)*Cos(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*Mh2)*
     Cos(alp + beta))*(pow(ML, 2))*(pow(Sec(beta), 2))*Tan(beta))/
  (Pi2*v2) + (0.015625*(m12 - m22 - m32)*MW2*C0i(cc1, m22, m32, m12, MA02, 
    ML2, MZ2)*Cos(alp - beta)*(pow(Sec(w), 2))*Tan(beta))/(Pi2*v2) + 
 (0.015625*(3*m12 - 3*m22 + m32)*MW2*C0i(cc1, m32, m22, m12, MA02, ML2, MZ2)*
   Cos(alp - beta)*(pow(Sec(w), 2))*Tan(beta))/(Pi2*v2) + 
 (0.015625*(m12 - m22 + m32)*MW2*C0i(cc2, m22, m32, m12, MA02, ML2, MZ2)*
   Cos(alp - beta)*(pow(Sec(w), 2))*Tan(beta))/(Pi2*v2) + 
 (0.015625*(5*m12 + m22 - m32)*MW2*C0i(cc2, m32, m22, m12, MA02, ML2, MZ2)*
   Cos(alp - beta)*(pow(Sec(w), 2))*Tan(beta))/(Pi2*v2) + 
 (0.015625*C0i(cc0, m22, m32, m12, MA02, ML2, MZ2)*Cos(alp - beta)*
   (m12*MW2 + m22*MW2 - m32*MW2 - 2*Mh2*(pow(ML, 2)) + 
    2*(MA02 - Mh2)*Cos(2*w)*(pow(ML, 2)) + 
    2*MA02*(-MW2 + (pow(ML, 2))))*(pow(Sec(w), 2))*Tan(beta))/
  (Pi2*v2) + (0.0078125*C0i(cc0, m32, m22, m12, MA02, ML2, MZ2)*
   Cos(alp - beta)*(3*m12*MW2 - m22*MW2 + m32*MW2 + 3*MA02*(pow(ML, 2)) - 
    3*Mh2*(pow(ML, 2)) + (MA02 - Mh2)*Cos(4*w)*(pow(ML, 2)) + 
    Cos(2*w)*((3*m12 - m22 + m32)*MW2 + 4*MA02*(pow(ML, 2)) - 
      4*Mh2*(pow(ML, 2))) - 4*(pow(MW, 4)))*(pow(Sec(w), 4))*
   Tan(beta))/(Pi2*v2) - (0.03125*(Mh2 - MHH2)*A0i(aa0, MW2)*Sec(beta)*
   Sin(alp)*Sin(2*(alp - beta))*Tan(beta))/(MA02*Pi2*v2) - 
 (0.046875*(Mh2 - MHH2)*A0i(aa0, MZ2)*Sec(beta)*Sin(alp)*Sin(2*(alp - beta))*
   Tan(beta))/(MA02*Pi2*v2) - (0.03125*MW2*B0i(bb1, MA02, Mh2, MZ2)*
   (pow(Sec(w), 2))*Sec(beta)*Sin(alp)*Sin(2*(alp - beta))*Tan(beta))/
  (Pi2*v2) + (0.03125*MW2*B0i(bb1, MA02, MHH2, MZ2)*(pow(Sec(w), 2))*
   Sec(beta)*Sin(alp)*Sin(2*(alp - beta))*Tan(beta))/(Pi2*v2) + 
 (0.0078125*Mh2*B0i(bb0, 0, Mh2, MZ2)*(MA02 - Mh2 + 2*MW2 + 
    (MA02 - Mh2)*Cos(2*w))*(pow(Sec(w), 2))*Sec(beta)*Sin(alp)*
   Sin(2*(alp - beta))*Tan(beta))/(MA02*Pi2*v2) - 
 (0.0078125*B0i(bb0, MA02, MHH2, MZ2)*(-(MHH2*(MHH2 - 2*MW2)) + 
    MA02*(MHH2 + 2*MW2) + (MA02 - MHH2)*MHH2*Cos(2*w))*(pow(Sec(w), 2))*
   Sec(beta)*Sin(alp)*Sin(2*(alp - beta))*Tan(beta))/(MA02*Pi2*v2) + 
 (0.0078125*MHH2*B0i(bb0, 0, MHH2, MZ2)*(-MA02 + MHH2 - 2*MW2 + 
    (-MA02 + MHH2)*Cos(2*w))*(pow(Sec(w), 2))*Sec(beta)*Sin(alp)*
   Sin(2*(alp - beta))*Tan(beta))/(MA02*Pi2*v2) + 
 (0.0078125*B0i(bb0, MA02, Mh2, MZ2)*(2*Mh2*MW2 + MA02*(Mh2 + 2*MW2) + 
    (MA02 - Mh2)*Mh2*Cos(2*w) - (pow(Mh, 4)))*(pow(Sec(w), 2))*
   Sec(beta)*Sin(alp)*Sin(2*(alp - beta))*Tan(beta))/(MA02*Pi2*v2) + 
 (0.015625*Re(A0i(aa0, MZ2))*Sec(beta)*
   (-((5 + 13*Cos(2*w) + 6*Cos(4*w))*(pow(Csc(w), 2))*
      (Cos(beta) + Sin(alp))) + (Mh2 - MHH2)*(pow(MA02, -1))*Sin(alp)*
     Sin(2*(alp - beta))*Tan(beta)))/(Pi2*v2) + 
 (0.0078125*Re(A0i(aa0, MW2))*Sec(beta)*
   (-((-22 + Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*(pow(Csc(w), 2))*
      (Cos(beta) + Sin(alp))) + 4*(Mh2 - MHH2)*(pow(MA02, -1))*Sin(alp)*
     Sin(2*(alp - beta))*Tan(beta)))/(Pi2*v2) + 
 (0.03125*B0i(bb0, m32, MA02, ML2)*
   (-(MW2*Cos(alp)*(pow(Sec(w), 2))*Sin(beta)) - 
    Sin(alp)*Tan(beta)*(MW2*(pow(Sec(w), 2))*Sin(beta) + 
      2*(pow(ML, 2))*Sec(beta)*Tan(beta))))/(Pi2*v2)
;
}

ComplexType dKappahtt(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32)
{
return 0. + (0.03125*(m12 - m22 - m32)*MW2*C0i(cc1, m22, m12, m32, MB2, MHp2, MW2)*
   Cos(alp - beta)*Cot(beta))/(Pi2*v2) + 
 (0.03125*(3*m12 - 3*m22 + m32)*MW2*C0i(cc1, m32, m12, m22, MB2, MHp2, MW2)*
   Cos(alp - beta)*Cot(beta))/(Pi2*v2) - 
 (0.0625*m32*MW2*C0i(cc2, m22, m12, m32, MB2, MHp2, MW2)*Cos(alp - beta)*
   Cot(beta))/(Pi2*v2) - (0.0625*(m12 + 2*m22 - m32)*MW2*
   C0i(cc2, m32, m12, m22, MB2, MHp2, MW2)*Cos(alp - beta)*Cot(beta))/
  (Pi2*v2) + 1.*(-1 + Cos(alp)*Csc(beta)) + 
 (0.125*m12*C0i(cc1, m12, m32, m22, MB2, MB2, MW2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MB, 2)))/(Pi2*v2) + 
 (0.0625*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MB2, MB2, MW2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MB, 2)))/(Pi2*v2) + 
 (0.0625*C0i(cc0, m22, m12, m32, MB2, MHp2, MW2)*Cos(alp - beta)*Cot(beta)*
   ((-m22 + MHp2)*MW2 + 2*(Mh2 - MHp2)*(pow(MB, 2))))/(Pi2*v2) + 
 (0.03125*B0i(bb0, m32, MB2, MW2)*(-(MW2*Cos(alp - 2*beta)*Csc(beta)) - 
    4*(MW2 + (pow(MB, 2))) + Cos(alp)*Csc(beta)*
     (3*MW2 + 4*(pow(MB, 2)))))/(Pi2*v2) + 
 (0.25*C0i(cc0, m12, m32, m22, MB2, MB2, MW2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MB, 4)))/(Pi2*v2) + (0.25*A0i(aa0, ME2)*Cos(alp)*Csc(beta)*
   (pow(ME, 2)))/(MA02*Pi2*v2) + 
 (0.125*B0i(bb1, MA02, ME2, ME2)*Cos(alp)*Csc(beta)*(pow(ME, 2)))/
  (Pi2*v2) + (0.25*A0i(aa0, ML2)*Cos(alp)*Csc(beta)*(pow(ML, 2)))/
  (MA02*Pi2*v2) + (0.125*B0i(bb1, MA02, ML2, ML2)*Cos(alp)*Csc(beta)*
   (pow(ML, 2)))/(Pi2*v2) + (0.25*A0i(aa0, MM2)*Cos(alp)*Csc(beta)*
   (pow(MM, 2)))/(MA02*Pi2*v2) + 
 (0.125*B0i(bb1, MA02, MM2, MM2)*Cos(alp)*Csc(beta)*(pow(MM, 2)))/
  (Pi2*v2) + (0.0625*C0i(cc0, m32, m12, m22, MB2, MHp2, MW2)*Cos(alp - beta)*
   Cot(beta)*(-(m22*MW2) + 2*(Mh2 - MHp2)*(pow(MB, 2)) + 
    (pow(MW, 4))))/(Pi2*v2) - (0.75*MT2*A0i(aa0, MT2)*Cos(alp)*Csc(beta)*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*MT2*B0i(bb1, MA02, MT2, MT2)*Cos(alp)*Csc(beta)*
   (pow(Cot(beta), 2)))/(Pi2*v2) + 
 (0.03125*MT2*(-m12 - m22 + m32 + 4*MT2)*C0i(cc0, m22, m12, m32, MA02, MT2, 
    MT2)*Cos(alp)*Csc(beta)*(pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.03125*(m12 + m22 - m32)*MT2*C0i(cc1, m22, m12, m32, MA02, MT2, MT2)*
   Cos(alp)*Csc(beta)*(pow(Cot(beta), 2)))/(Pi2*v2) + 
 (0.03125*(m12 - m22 + m32)*MT2*C0i(cc2, m22, m12, m32, MA02, MT2, MT2)*
   Cos(alp)*Csc(beta)*(pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MB2)*Cos(alp)*Csc(beta)*(pow(MB, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MB2, MB2)*Cos(alp)*Csc(beta)*(pow(MB, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) + 
 (0.125*m12*C0i(cc1, m12, m32, m22, MB2, MB2, MHp2)*Cos(alp)*Csc(beta)*
   (pow(MB, 2))*(pow(Cot(beta), 2)))/(Pi2*v2) + 
 (0.0625*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MB2, MB2, MHp2)*Cos(alp)*
   Csc(beta)*(pow(MB, 2))*(pow(Cot(beta), 2)))/(Pi2*v2) + 
 (0.25*C0i(cc0, m12, m32, m22, MB2, MB2, MHp2)*Cos(alp)*Csc(beta)*
   (pow(MB, 4))*(pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MC2)*Cos(alp)*Csc(beta)*(pow(MC, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MC2, MC2)*Cos(alp)*Csc(beta)*(pow(MC, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MD2)*Cos(alp)*Csc(beta)*(pow(MD, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MD2, MD2)*Cos(alp)*Csc(beta)*(pow(MD, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MS2)*Cos(alp)*Csc(beta)*(pow(MS, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MS2, MS2)*Cos(alp)*Csc(beta)*(pow(MS, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.75*A0i(aa0, MU2)*Cos(alp)*Csc(beta)*(pow(MU, 2))*
   (pow(Cot(beta), 2)))/(MA02*Pi2*v2) - 
 (0.375*B0i(bb1, MA02, MU2, MU2)*Cos(alp)*Csc(beta)*(pow(MU, 2))*
   (pow(Cot(beta), 2)))/(Pi2*v2) - 
 (0.015625*MT2*C0i(cc0, m12, m32, m22, MA02, MA02, MT2)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*Cot(beta)*(pow(Csc(beta), 2)))/(Pi2*v2) + 
 (0.03125*C0i(cc0, m22, m12, m32, MB2, MHp2, MHp2)*
   ((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*Cot(beta)*(pow(MB, 2))*(pow(Csc(beta), 2)))/
  (Pi2*v2) + (0.01171875*A0i(aa0, MA02)*Cos(alp)*(3*(Mh2 - MHH2)*Cos(2*alp) + 
    (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 4*(-2*M2 + Mh2 + MHH2)*Cos(2*beta))*
   (pow(Csc(beta), 3)))/(MA02*Pi2*v2) + 
 (0.0078125*A0i(aa0, MHp2)*Cos(alp)*(3*(Mh2 - MHH2)*Cos(2*alp) + 
    (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 4*(-2*M2 + Mh2 + MHH2)*Cos(2*beta))*
   (pow(Csc(beta), 3)))/(MA02*Pi2*v2) + 
 (0.0078125*(MA02 - Mh2)*B0i(bb0, 0, MA02, Mh2)*Cos(alp)*Cos(alp - beta)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*(pow(Csc(beta), 3)))/(MA02*Pi2*v2) + 
 (0.0078125*(MA02 - Mh2)*B0i(bb0, MA02, MA02, Mh2)*Cos(alp)*Cos(alp - beta)*
   ((2*MA02 - Mh2)*Cos(alp - 3*beta) + (4*M2 - 2*MA02 - 3*Mh2)*
     Cos(alp + beta))*(pow(Csc(beta), 3)))/(MA02*Pi2*v2) + 
 (0.0078125*A0i(aa0, Mh2)*Cos(alp)*Cos(alp - beta)*
   (-2*(M2 - MA02)*Cos(alp - 3*beta) + (Mh2 - MHH2)*Cos(3*alp - beta) + 
    (-2*M2 - 2*MA02 + 3*Mh2 + MHH2)*Cos(alp + beta))*(pow(Csc(beta), 3)))/
  (MA02*Pi2*v2) - (0.0078125*MT2*(-m12 - m22 + m32 + 4*MT2)*
   C0i(cc0, m22, m12, m32, Mh2, MT2, MT2)*
   (-4 + 3*Cos(alp)*(pow(Csc(beta), 3)) + 
    Cos(3*alp)*(pow(Csc(beta), 3))))/(Pi2*v2) + 
 (0.0078125*(m12 + m22 - m32)*MT2*C0i(cc1, m22, m12, m32, Mh2, MT2, MT2)*
   (-4 + 3*Cos(alp)*(pow(Csc(beta), 3)) + 
    Cos(3*alp)*(pow(Csc(beta), 3))))/(Pi2*v2) - 
 (0.0078125*(m12 - m22 + m32)*MT2*C0i(cc2, m22, m12, m32, Mh2, MT2, MT2)*
   (-4 + 3*Cos(alp)*(pow(Csc(beta), 3)) + 
    Cos(3*alp)*(pow(Csc(beta), 3))))/(Pi2*v2) - 
 (0.0625*MT2*B0i(bb0, m32, Mh2, MT2)*
   (-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3))))/(Pi2*v2) - 
 (0.015625*(m12 - m22 - m32)*MW2*C0i(cc1, m22, m32, m12, MA02, MT2, MZ2)*
   Cos(alp - beta)*Cot(beta)*(pow(Sec(w), 2)))/(Pi2*v2) - 
 (0.015625*(3*m12 - 3*m22 + m32)*MW2*C0i(cc1, m32, m22, m12, MA02, MT2, MZ2)*
   Cos(alp - beta)*Cot(beta)*(pow(Sec(w), 2)))/(Pi2*v2) - 
 (0.015625*(m12 - m22 + m32)*MW2*C0i(cc2, m22, m32, m12, MA02, MT2, MZ2)*
   Cos(alp - beta)*Cot(beta)*(pow(Sec(w), 2)))/(Pi2*v2) - 
 (0.015625*(5*m12 + m22 - m32)*MW2*C0i(cc2, m32, m22, m12, MA02, MT2, MZ2)*
   Cos(alp - beta)*Cot(beta)*(pow(Sec(w), 2)))/(Pi2*v2) - 
 (0.015625*C0i(cc0, m22, m32, m12, MA02, MT2, MZ2)*Cos(alp - beta)*
   (-2*Mh2*MT2 + 2*MA02*(MT2 - MW2) + m12*MW2 + m22*MW2 - m32*MW2 + 
    2*(MA02 - Mh2)*MT2*Cos(2*w))*Cot(beta)*(pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.006944444444444444*MT2*C0i(cc0, m12, m32, m22, MT2, MT2, MZ2)*
   (9*MT2 + (9*MT2 - 16*MW2)*Cos(2*w) + 16*MW2*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.003472222222222222*m12*C0i(cc1, m12, m32, m22, MT2, MT2, MZ2)*
   (9*MT2 + (9*MT2 - 16*MW2)*Cos(2*w) + 16*MW2*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.001736111111111111*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MT2, MT2, 
    MZ2)*(9*MT2 + (9*MT2 - 16*MW2)*Cos(2*w) + 16*MW2*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.001736111111111111*B0i(bb0, m32, MT2, MZ2)*
   (-36*MT2 + Cos(alp)*(18*MT2 + 27*MW2 + 2*(9*MT2 - 16*MW2)*Cos(2*w) + 
      32*MW2*Cos(4*w))*Csc(beta)*(pow(Sec(w), 2)) - 
    MW2*(4 - 32*Cos(2*w) + 9*Cos(alp - 2*beta)*Csc(beta) + 
      64*(pow(Cos(2*w), 2)))*(pow(Sec(w), 2))))/(Pi2*v2) - 
 (0.0078125*C0i(cc0, m32, m22, m12, MA02, MT2, MZ2)*Cos(alp - beta)*Cot(beta)*
   (3*MA02*MT2 - 3*Mh2*MT2 + 3*m12*MW2 - m22*MW2 + m32*MW2 + 
    (4*MA02*MT2 - 4*Mh2*MT2 + (3*m12 - m22 + m32)*MW2)*Cos(2*w) + 
    (MA02 - Mh2)*MT2*Cos(4*w) - 4*(pow(MW, 4)))*(pow(Sec(w), 4)))/
  (Pi2*v2) - (0.0625*MT2*B0i(bb0, m32, MHH2, MT2)*Cos(alp)*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2)))/(Pi2*v2) + 
 (0.03125*(m12 + m22 - m32 - 4*MT2)*MT2*C0i(cc0, m22, m12, m32, MHH2, MT2, 
    MT2)*Cos(alp)*(pow(Csc(beta), 3))*(pow(Sin(alp), 2)))/(Pi2*v2) + 
 (0.03125*(m12 + m22 - m32)*MT2*C0i(cc1, m22, m12, m32, MHH2, MT2, MT2)*
   Cos(alp)*(pow(Csc(beta), 3))*(pow(Sin(alp), 2)))/(Pi2*v2) - 
 (0.03125*(m12 - m22 + m32)*MT2*C0i(cc2, m22, m12, m32, MHH2, MT2, MT2)*
   Cos(alp)*(pow(Csc(beta), 3))*(pow(Sin(alp), 2)))/(Pi2*v2) + 
 (0.03125*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, MB2, MW2, MW2)*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.125*C0i(cc0, m22, m12, m32, MB2, MW2, MW2)*
   (-(m22*MW2) + Mh2*(pow(MB, 2)) + (pow(MW, 4)))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) + 
 (0.015625*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, MT2, MZ2, MZ2)*
   (pow(Sec(w), 2))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 
     2)))/(Pi2*v2) - (0.0008680555555555555*C0i(cc0, m22, m12, m32, MT2, MZ2, 
    MZ2)*(4*Cos(2*w)*(9*Mh2*MT2 - 9*m22*MW2 - 32*(pow(MW, 4))) + 
    9*(3*Mh2*MT2 - 4*m22*MW2 + 8*(pow(MW, 4))) + 
    Cos(4*w)*(9*Mh2*MT2 + 128*(pow(MW, 4))))*(pow(Sec(w), 4))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2)))/(Pi2*v2) - 
 (0.2222222222222222*(m12 + m22 - m32 - 4*MT2)*C0i(cc0, m22, m12, m32, 0, 
    MT2, MT2)*(-1 + Cos(alp)*Csc(beta))*(3*Alfas*Pi*v2 + 
    MW2*(pow(Sin(w), 2))))/(Pi2*v2) - 
 (0.2222222222222222*(m12 + m22 - m32)*C0i(cc1, m22, m12, m32, 0, MT2, MT2)*
   (-1 + Cos(alp)*Csc(beta))*(3*Alfas*Pi*v2 + MW2*(pow(Sin(w), 2))))/
  (Pi2*v2) + (0.2222222222222222*(m12 - m22 + m32)*
   C0i(cc2, m22, m12, m32, 0, MT2, MT2)*(-1 + Cos(alp)*Csc(beta))*
   (3*Alfas*Pi*v2 + MW2*(pow(Sin(w), 2))))/(Pi2*v2) + 
 (0.1111111111111111*B0i(bb0, m32, 0, MT2)*(-1 + Cos(alp)*Csc(beta))*
   (12*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2))))/Pi2 - 
 (0.00390625*Cos(alp)*(3*(Mh2 - MHH2)*Cos(2*alp) + 
    (Mh2 - MHH2)*Cos(2*(alp - 2*beta)) + 
    4*(MA02 + (-2*M2 - MA02 + Mh2 + MHH2)*Cos(2*beta)))*
   (pow(Csc(beta), 3))*Re(A0i(aa0, MA02)))/(MA02*Pi2*v2) + 
 (0.005208333333333333*(144*Csc(beta)*(pow(MB, 2))*(pow(Cos(alp), 3))*
     (pow(Cot(beta), 2)) + 
    MA02*(61 + 4*Cos(alp)*(-12 + 2*Cos(2*w) + Cos(4*w))*Csc(beta))*
     (pow(Cot(w), 2)) + (MA02*(-93 + 61*Cos(2*w) - 20*Cos(4*w) + 
       2*Cos(6*w))*(pow(Csc(w), 2)))/2 + 8*Cos(alp)*Csc(beta)*
     (MA02*(9 - Cos(2*w) + Cos(4*w)) + 18*(pow(MB, 2))*
       (pow(Cot(beta), 2))*(pow(Sin(alp), 2))))*Re(A0i(aa0, MB2)))/
  (MA02*Pi2*v2) + (0.010416666666666666*
   (2*Cos(alp)*Csc(beta)*(64*MA02 - 24*MA02*(pow(Csc(w), 2)) + 
      (36*(pow(MC, 2))*(pow(Cot(beta), 2)) + 
        MA02*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2)))*
       (pow(Csc(w), 4))) + MA02*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 
      2*Cos(6*w))*(pow(Csc(w), 6)))*Re(A0i(aa0, MC2)))/
  (MA02*Pi2*v2*(pow(Csc(w), 4))) + 
 (0.005208333333333333*(144*Csc(beta)*(pow(MD, 2))*(pow(Cos(alp), 3))*
     (pow(Cot(beta), 2)) + 
    MA02*(61 + 4*Cos(alp)*(-12 + 2*Cos(2*w) + Cos(4*w))*Csc(beta))*
     (pow(Cot(w), 2)) + (MA02*(-93 + 61*Cos(2*w) - 20*Cos(4*w) + 
       2*Cos(6*w))*(pow(Csc(w), 2)))/2 + 8*Cos(alp)*Csc(beta)*
     (MA02*(9 - Cos(2*w) + Cos(4*w)) + 18*(pow(MD, 2))*
       (pow(Cot(beta), 2))*(pow(Sin(alp), 2))))*Re(A0i(aa0, MD2)))/
  (MA02*Pi2*v2) - (0.015625*(16*MA02 - MA02*Cos(6*w) - 
    16*MA02*Cos(alp)*Csc(beta) + MA02*Cos(alp)*Cos(6*w)*Csc(beta) - 
    10*MA02*Cos(4*w)*(-1 + Cos(alp)*Csc(beta)) + 
    Cos(2*w)*(-29*MA02 + Cos(alp)*Csc(beta)*(29*MA02 - 8*(pow(ME, 2)))) + 
    8*Cos(alp)*Csc(beta)*(pow(ME, 2)))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ME2)))/(MA02*Pi2*v2) - 
 (0.00390625*Cos(alp)*(4*MA02 + 2*(5*M2 - 6*MHH2)*Cos(2*alp) + 
    2*M2*Cos(2*(alp - 2*beta)) - Mh2*Cos(4*alp - 2*beta) + 
    MHH2*Cos(4*alp - 2*beta) - 12*M2*Cos(2*beta) - 4*MA02*Cos(2*beta) + 
    Mh2*Cos(2*beta) + 11*MHH2*Cos(2*beta))*(pow(Csc(beta), 3))*
   Re(A0i(aa0, MHH2)))/(MA02*Pi2*v2) + 
 (0.0078125*Cos(alp)*(-2*MA02 - 3*(Mh2 - MHH2)*Cos(2*alp) + 
    (-Mh2 + MHH2)*Cos(2*(alp - 2*beta)) + 8*M2*Cos(2*beta) + 
    2*MA02*Cos(2*beta) - 4*Mh2*Cos(2*beta) - 4*MHH2*Cos(2*beta) + 
    MA02*Cos(2*(beta - 2*w)) - 8*MA02*Cos(2*(beta - w)) + 16*MA02*Cos(2*w) - 
    2*MA02*Cos(4*w) - 8*MA02*Cos(2*(beta + w)) + MA02*Cos(2*(beta + 2*w)))*
   (pow(Csc(beta), 3))*Re(A0i(aa0, MHp2)))/(MA02*Pi2*v2) - 
 (0.015625*(16*MA02 - MA02*Cos(6*w) - 16*MA02*Cos(alp)*Csc(beta) + 
    MA02*Cos(alp)*Cos(6*w)*Csc(beta) - 10*MA02*Cos(4*w)*
     (-1 + Cos(alp)*Csc(beta)) + Cos(2*w)*(-29*MA02 + 
      Cos(alp)*Csc(beta)*(29*MA02 - 8*(pow(ML, 2)))) + 
    8*Cos(alp)*Csc(beta)*(pow(ML, 2)))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ML2)))/(MA02*Pi2*v2) - 
 (0.015625*(16*MA02 - MA02*Cos(6*w) - 16*MA02*Cos(alp)*Csc(beta) + 
    MA02*Cos(alp)*Cos(6*w)*Csc(beta) - 10*MA02*Cos(4*w)*
     (-1 + Cos(alp)*Csc(beta)) + Cos(2*w)*(-29*MA02 + 
      Cos(alp)*Csc(beta)*(29*MA02 - 8*(pow(MM, 2)))) + 
    8*Cos(alp)*Csc(beta)*(pow(MM, 2)))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MM2)))/(MA02*Pi2*v2) + 
 (0.005208333333333333*(144*Csc(beta)*(pow(MS, 2))*(pow(Cos(alp), 3))*
     (pow(Cot(beta), 2)) + 
    MA02*(61 + 4*Cos(alp)*(-12 + 2*Cos(2*w) + Cos(4*w))*Csc(beta))*
     (pow(Cot(w), 2)) + (MA02*(-93 + 61*Cos(2*w) - 20*Cos(4*w) + 
       2*Cos(6*w))*(pow(Csc(w), 2)))/2 + 8*Cos(alp)*Csc(beta)*
     (MA02*(9 - Cos(2*w) + Cos(4*w)) + 18*(pow(MS, 2))*
       (pow(Cot(beta), 2))*(pow(Sin(alp), 2))))*Re(A0i(aa0, MS2)))/
  (MA02*Pi2*v2) + (0.010416666666666666*
   (2*Cos(alp)*Csc(beta)*(64*MA02 - 24*MA02*(pow(Csc(w), 2)) + 
      (36*MT2*(pow(Cot(beta), 2)) + MA02*(9 - 4*Cos(2*w) + 4*Cos(4*w))*
         (pow(Cot(w), 2)))*(pow(Csc(w), 4))) + 
    MA02*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*(pow(Csc(w), 6)))*
   Re(A0i(aa0, MT2)))/(MA02*Pi2*v2*(pow(Csc(w), 4))) + 
 (0.010416666666666666*(2*Cos(alp)*Csc(beta)*
     (64*MA02 - 24*MA02*(pow(Csc(w), 2)) + 
      (36*(pow(MU, 2))*(pow(Cot(beta), 2)) + 
        MA02*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2)))*
       (pow(Csc(w), 4))) + MA02*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 
      2*Cos(6*w))*(pow(Csc(w), 6)))*Re(A0i(aa0, MU2)))/
  (MA02*Pi2*v2*(pow(Csc(w), 4))) - 
 (1.125*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb0, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.75*Cos(alp)*(pow(MB, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, MB2, MB2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.75*Cos(alp)*(pow(MC, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, Mh2, MC2, MC2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*(pow(MD, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, MD2, MD2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.25*Cos(alp)*Csc(beta)*(pow(ME, 4))*
   (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, Mh2, ME2, ME2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.25*Cos(alp)*Csc(beta)*(pow(ML, 4))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, ML2, ML2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.25*Cos(alp)*Csc(beta)*(pow(MM, 4))*
   (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, Mh2, MM2, MM2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*(pow(MS, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, MS2, MS2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.75*Cos(alp)*(pow(MT, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, Mh2, MT2, MT2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.75*Cos(alp)*(pow(MU, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, Mh2, MU2, MU2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.75*Cos(alp)*(pow(MB, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, MB2, MB2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*(pow(MC, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, MHH2, MC2, MC2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.75*Cos(alp)*(pow(MD, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, MD2, MD2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*Cos(alp)*Csc(beta)*(pow(ME, 4))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, MHH2, ME2, ME2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.25*Cos(alp)*Csc(beta)*(pow(ML, 4))*
   (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, ML2, ML2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.25*Cos(alp)*Csc(beta)*(pow(MM, 4))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, MHH2, MM2, MM2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.75*Cos(alp)*(pow(MS, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, MS2, MS2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.75*Cos(alp)*(pow(MT, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb0, MHH2, MT2, MT2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.75*Cos(alp)*(pow(MU, 4))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MHH2, MU2, MU2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.1111111111111111*(-1 + Cos(alp)*Csc(beta))*
   (12*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*
   Re(B0i(bb0, MT2, 0, MT2)))/Pi2 - 
 (0.0625*MT2*Cos(alp)*Csc(beta)*(pow(Cot(beta), 2))*
   Re(B0i(bb0, MT2, MA02, MT2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Csc(beta)*(pow(MB, 2))*(pow(Cot(beta), 2))*
   Re(B0i(bb0, MT2, MB2, MHp2)))/(Pi2*v2) - 
 (0.125*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*
   Re(B0i(bb0, MT2, MB2, MW2)))/(Pi2*v2) + 
 (0.0625*MT2*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(bb0, MT2, Mh2, MT2)))/(Pi2*v2) + 
 (0.0625*MT2*Cos(alp)*(pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb0, MT2, MHH2, MT2)))/(Pi2*v2) - 
 (0.003472222222222222*(9*MT2 + (9*MT2 - 16*MW2)*Cos(2*w) + 16*MW2*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Sec(w), 2))*
   Re(B0i(bb0, MT2, MT2, MZ2)))/(Pi2*v2) - 
 (0.5*MW2*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*Re(B0i(bb0, MW2, 0, MW2)))/
  (Pi2*v2) - (0.375*(MT2 - MW2)*(-1 + Cos(alp)*Csc(beta))*
   (-1 + (pow(Cot(w), 2)))*Re(B0i(bb0, MW2, MB2, MT2)))/(Pi2*v2) - 
 (0.375*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*
   (-1 + (pow(Cot(w), 2)))*Re(B0i(bb0, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*(MW2 - (pow(MU, 2)))*
   (pow(Csc(w), 2))*Re(B0i(bb0, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.125*MW2*Cos(2*w)*(pow(Csc(w), 2))*
   (-1 + Cos(alp)*Csc(beta)*(pow(Sin(alp - beta), 2)))*
   Re(B0i(bb0, MW2, Mh2, MW2)))/(Pi2*v2) + 
 (0.125*MW2*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb0, MW2, MHH2, MW2)))/(Pi2*v2) - 
 (0.015625*MW2*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*(pow(Sec(w), 2))*
   Re(B0i(bb0, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MD, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*(-1 + Cos(alp)*Csc(beta))*(pow(ME, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, ME2, ME2)))/(Pi2*v2) - 
 (0.125*MW2*(pow(Csc(w), 2))*(-1 + Cos(alp)*Csc(beta)*
     (pow(Sin(alp - beta), 2)))*Re(B0i(bb0, MZ2, Mh2, MZ2)))/(Pi2*v2) - 
 (0.125*MW2*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb0, MZ2, MHH2, MZ2)))/(Pi2*v2) + 
 (0.0625*(-1 + Cos(alp)*Csc(beta))*(pow(ML, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*(-1 + Cos(alp)*Csc(beta))*(pow(MM, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MS, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.1875*MT2*(-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.1875*(-1 + Cos(alp)*Csc(beta))*(pow(MU, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.0625*MW2*(5 + 9*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MB2, MB2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MD2, MD2)))/(Pi2*v2) + 
 (0.5*(-1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, ME2, ME2)))/(Pi2*v2) - 
 (0.5*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MHp2, MHp2)))/(Pi2*v2) + 
 (0.5*(-1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, ML2, ML2)))/(Pi2*v2) + 
 (0.5*(-1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MM2, MM2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MT2, MT2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2))*Re(B0i(bb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.5*(2 + 3*Cos(2*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.25*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, 0, ME2)))/(Pi2*v2) + 
 (0.25*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, 0, ML2)))/(Pi2*v2) + 
 (0.25*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, 0, MM2)))/(Pi2*v2) - 
 (1.*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*Re(B0i(bb00, MW2, 0, MW2)))/
  (Pi2*v2) - (0.125*Cos(alp)*Csc(beta)*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, MA02, MHp2)))/(Pi2*v2) + 
 (0.75*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, MB2, MT2)))/(Pi2*v2) + 
 (0.75*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.75*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb00, MW2, MD2, MU2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MW2, Mh2, MHp2)))/(Pi2*v2) - 
 (0.125*Cos(2*w)*(pow(Csc(w), 2))*
   (-1 + Cos(alp)*Csc(beta)*(pow(Sin(alp - beta), 2)))*
   Re(B0i(bb00, MW2, Mh2, MW2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Csc(w), 2))*
   (pow(Sin(alp - beta), 2))*Re(B0i(bb00, MW2, MHH2, MHp2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*Cos(2*w)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MW2, MHH2, MW2)))/(Pi2*v2) - 
 (0.125*(2 + 5*Cos(2*w) + 2*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MW2, MW2, MZ2)))/(Pi2*v2) - 
 (0.375*(-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, 0, 0)))/(Pi2*v2) + 
 (0.125*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MA02, Mh2)))/(Pi2*v2) + 
 (0.125*Cos(alp)*Csc(beta)*(pow(Cot(w), 2))*(pow(Sin(alp - beta), 2))*
   Re(B0i(bb00, MZ2, MA02, MHH2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Cot(w), 2))*Re(B0i(bb00, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.041666666666666664*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MC2, MC2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Cot(w), 2))*Re(B0i(bb00, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*(pow(Cot(w), 2))*(-1 + Cos(alp)*Csc(beta)*
     (pow(Sin(alp - beta), 2)))*Re(B0i(bb00, MZ2, Mh2, MZ2)))/(Pi2*v2) + 
 (0.125*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MHH2, MZ2)))/(Pi2*v2) + 
 (0.125*Cos(alp)*Csc(beta)*(pow(Cos(2*w), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MHp2, MHp2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Cot(w), 2))*Re(B0i(bb00, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.041666666666666664*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.041666666666666664*(-9 + 4*Cos(2*w) - 4*Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.0625*(7 + 8*Cos(2*w) + 3*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*
   (pow(Cot(w), 2))*Re(B0i(bb00, MZ2, MW2, MW2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MB2, MB2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MC2, MC2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MD2, MD2)))/(Pi2*v2) - 
 (0.5*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, ME2, ME2)))/(Pi2*v2) - 
 (0.5*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, ML2, ML2)))/(Pi2*v2) - 
 (0.5*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MM2, MM2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MS2, MS2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MT2, MT2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.25*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.03125*MW2*Cos(alp - beta)*(MHH2*Cos(2*alp - beta) + 
    (-2*Mh2 + MHH2)*Cos(beta))*Csc(beta)*(pow(Sec(w), 2))*
   Re(B0i(bb1, Mh2, MA02, MZ2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MB, 2))*(-2*Mh2 + 2*MHH2 - 
    Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*(pow(Csc(beta), 3)))*
   Re(B0i(bb1, Mh2, MB2, MB2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MC, 2))*(-2*Mh2 + 2*MHH2 - 
    Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*(pow(Csc(beta), 3)))*
   Re(B0i(bb1, Mh2, MC2, MC2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*(pow(MD, 2))*(-2*Mh2 + 2*MHH2 - 
    Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*(pow(Csc(beta), 3)))*
   Re(B0i(bb1, Mh2, MD2, MD2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.125*(pow(ME, 2))*(Mh2 - MHH2 + MHH2*Cos(alp)*Csc(beta)*
     (pow(Sec(beta), 2))*(pow(Sin(alp), 2)))*
   Re(B0i(bb1, Mh2, ME2, ME2)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0625*MW2*Cos(alp - beta)*(MHH2*Cos(2*alp - beta) + 
    (-2*Mh2 + MHH2)*Cos(beta))*Csc(beta)*Re(B0i(bb1, Mh2, MHp2, MW2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.125*(pow(ML, 2))*
   (Mh2 - MHH2 + MHH2*Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(bb1, Mh2, ML2, ML2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.125*(pow(MM, 2))*
   (Mh2 - MHH2 + MHH2*Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(bb1, Mh2, MM2, MM2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.1875*(pow(MS, 2))*
   (-2*Mh2 + 2*MHH2 - Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*
     (pow(Csc(beta), 3)))*Re(B0i(bb1, Mh2, MS2, MS2)))/
  ((Mh2 - MHH2)*Pi2*v2) + 
 (0.1875*MT2*(-2*Mh2 + 2*MHH2 - Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*
     (pow(Csc(beta), 3)))*Re(B0i(bb1, Mh2, MT2, MT2)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.1875*(pow(MU, 2))*
   (-2*Mh2 + 2*MHH2 - Cos(alp)*(-2*Mh2 + MHH2 + MHH2*Cos(2*alp))*
     (pow(Csc(beta), 3)))*Re(B0i(bb1, Mh2, MU2, MU2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.375*MHH2*Cos(alp)*(pow(MB, 2))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, MB2, MB2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*(pow(MC, 2))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb1, MHH2, MC2, MC2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.375*MHH2*Cos(alp)*(pow(MD, 2))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, MD2, MD2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.125*MHH2*Cos(alp)*Csc(beta)*(pow(ME, 2))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb1, MHH2, ME2, ME2)))/
  ((-Mh2 + MHH2)*Pi2*v2) - (0.125*MHH2*Cos(alp)*Csc(beta)*(pow(ML, 2))*
   (pow(Sec(beta), 2))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, ML2, ML2)))/((-Mh2 + MHH2)*Pi2*v2) - 
 (0.125*MHH2*Cos(alp)*Csc(beta)*(pow(MM, 2))*(pow(Sec(beta), 2))*
   (pow(Sin(alp), 2))*Re(B0i(bb1, MHH2, MM2, MM2)))/
  ((-Mh2 + MHH2)*Pi2*v2) - (0.375*MHH2*Cos(alp)*(pow(MS, 2))*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, MS2, MS2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*MT2*Cos(alp)*(pow(Csc(beta), 3))*(pow(Sin(alp), 2))*
   Re(B0i(bb1, MHH2, MT2, MT2)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.375*MHH2*Cos(alp)*(pow(MU, 2))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(bb1, MHH2, MU2, MU2)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.125*MW2*(-1 + Cos(alp)*Csc(beta))*
   (-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ME2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, 0, ML2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, 0, MM2)))/(Pi2*v2) - 
 (0.25*MW2*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*Re(B0i(bb1, MW2, 0, MW2)))/
  (Pi2*v2) + (0.375*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MB2, MT2)))/(Pi2*v2) - 
 (0.375*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + Cos(alp)*Csc(beta))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.25*MW2*Cos(2*w)*(-1 + Cos(alp)*Csc(beta))*(pow(Cot(w), 2))*
   Re(B0i(bb1, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.1875*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, 0, 0)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*
   (-1 + Cos(alp)*Csc(beta))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.25*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Cos(w), 2))*
   (pow(Cot(w), 2))*Re(B0i(bb1, MZ2, MW2, MW2)))/(Pi2*v2) - 
 (0.0009765625*Cos(alp)*(pow((2*MA02 - Mh2)*Cos(alp - 3*beta) + 
      (4*M2 - 2*MA02 - 3*Mh2)*Cos(alp + beta), 2))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*Re(B0i(dbb0, Mh2, MA02, MA02)))/(Pi2*v2) + 
 (0.75*(pow(MB, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.75*(pow(MC, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.75*(pow(MD, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.25*(pow(ME, 4))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb0, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.0087890625*(16*(pow(Mh, 4)) - 
    Cos(alp)*(pow(M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
        (2*M2 - 3*Mh2)*Cos(alp + beta), 2))*(pow(Csc(beta), 3))*
     (pow(Sec(beta), 2)))*Re(B0i(dbb0, Mh2, Mh2, Mh2)))/(Pi2*v2) - 
 (0.0078125*Cos(alp)*(pow(Cos(alp - beta), 2))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*(pow((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + 
      M2*Sin(2*beta), 2))*Re(B0i(dbb0, Mh2, Mh2, MHH2)))/(Pi2*v2) - 
 (0.00390625*Cos(alp)*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   (pow(Sin(alp - beta), 2))*
   (pow((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta), 2))*
   Re(B0i(dbb0, Mh2, MHH2, MHH2)))/(Pi2*v2) - 
 (0.001953125*Cos(alp)*(pow((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + 
      (-4*M2 + 3*Mh2 + 2*MHp2)*Cos(alp + beta), 2))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*Re(B0i(dbb0, Mh2, MHp2, MHp2)))/(Pi2*v2) - 
 (0.0625*Cos(alp)*Csc(beta)*(MHp2*(MHp2 - MW2) - Mh2*(2*MHp2 + MW2) + 
    (pow(Mh, 4)))*(pow(Cos(alp - beta), 2))*
   Re(B0i(dbb0, Mh2, MHp2, MW2)))/(Pi2*v2) + 
 (0.25*(pow(ML, 4))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb0, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.25*(pow(MM, 4))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb0, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.75*(pow(MS, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.75*(pow(MT, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.75*(pow(MU, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, Mh2, MU2, MU2)))/(Pi2*v2) + 
 (0.1111111111111111*MT2*(-1 + Cos(alp)*Csc(beta))*
   (12*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*
   Re(B0i(dbb0, MT2, 0, MT2)))/Pi2 + 
 (0.25*MT2*Cos(alp)*Csc(beta)*(pow(MB, 2))*(pow(Cot(beta), 2))*
   Re(B0i(dbb0, MT2, MB2, MHp2)))/(Pi2*v2) + 
 (0.25*MT2*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*
   Re(B0i(dbb0, MT2, MB2, MW2)))/(Pi2*v2) - 
 (0.25*(pow(MT, 4))*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb0, MT2, Mh2, MT2)))/(Pi2*v2) - 
 (0.25*Cos(alp)*(pow(MT, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(dbb0, MT2, MHH2, MT2)))/(Pi2*v2) + 
 (0.006944444444444444*MT2*(9*MT2 + (9*MT2 - 16*MW2)*Cos(2*w) + 
    16*MW2*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MT2, MT2, MZ2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MB2, MB2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MD2, MD2)))/(Pi2*v2) + 
 (1.*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, ME2, ME2)))/(Pi2*v2) - 
 (0.5*MW2*Cos(alp)*Csc(beta)*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MHp2, MHp2)))/(Pi2*v2) + 
 (1.*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, ML2, ML2)))/(Pi2*v2) + 
 (1.*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MM2, MM2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MT2, MT2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (1.5*MW2*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*Mh2*MW2*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   (pow(Sec(w), 2))*Re(B0i(dbb1, Mh2, MA02, MZ2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MB, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MC, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MD, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ME, 2))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb1, Mh2, ME2, ME2)))/(Pi2*v2) - 
 (0.125*Mh2*MW2*Cos(alp)*Csc(beta)*(pow(Cos(alp - beta), 2))*
   Re(B0i(dbb1, Mh2, MHp2, MW2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ML, 2))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(MM, 2))*(-1 + Cos(alp)*Csc(beta)*(pow(Sec(beta), 2))*
     (pow(Sin(alp), 2)))*Re(B0i(dbb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MS, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*Mh2*MT2*(-1 + (pow(Cos(alp), 3))*(pow(Csc(beta), 3)))*
   Re(B0i(dbb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MU, 2))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*Mh2*MW2*(-1 + Cos(alp)*Csc(beta)*(pow(Sin(alp - beta), 2)))*
   Re(B0i(dbb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*Mh2*MW2*(pow(Sec(w), 2))*
   (-1 + Cos(alp)*Csc(beta)*(pow(Sin(alp - beta), 2)))*
   Re(B0i(dbb1, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (0.1111111111111111*MT2*(-1 + Cos(alp)*Csc(beta))*
   (12*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*
   Re(B0i(dbb1, MT2, 0, MT2)))/Pi2 - 
 (0.125*Cos(alp)*Csc(beta)*(pow(MT, 4))*(pow(Cot(beta), 2))*
   Re(B0i(dbb1, MT2, MA02, MT2)))/(Pi2*v2) + 
 (0.125*MT2*Cos(alp)*Csc(beta)*(MT2 + (pow(MB, 2)))*
   (pow(Cot(beta), 2))*Re(B0i(dbb1, MT2, MB2, MHp2)))/(Pi2*v2) + 
 (0.125*MT2*(-1 + Cos(alp)*Csc(beta))*(MT2 + 2*MW2 + (pow(MB, 2)))*
   Re(B0i(dbb1, MT2, MB2, MW2)))/(Pi2*v2) - 
 (0.125*(pow(MT, 4))*(-1 + (pow(Cos(alp), 3))*
     (pow(Csc(beta), 3)))*Re(B0i(dbb1, MT2, Mh2, MT2)))/(Pi2*v2) - 
 (0.125*Cos(alp)*(pow(MT, 4))*(pow(Csc(beta), 3))*
   (pow(Sin(alp), 2))*Re(B0i(dbb1, MT2, MHH2, MT2)))/(Pi2*v2) + 
 (0.006944444444444444*MT2*(9*(MT2 + 2*MW2) + (9*MT2 - 8*MW2)*Cos(2*w) + 
    8*MW2*Cos(4*w))*(-1 + Cos(alp)*Csc(beta))*(pow(Sec(w), 2))*
   Re(B0i(dbb1, MT2, MT2, MZ2)))/(Pi2*v2) - 
 (0.046875*MT2*C0i(cc0, m12, m32, m22, Mh2, Mh2, MT2)*
   (-4*Mh2 - (M2*Cos(alp - 3*beta) + (M2 - Mh2)*Cos(3*alp - beta) + 
      (2*M2 - 3*Mh2)*Cos(alp + beta))*(pow(Cos(alp), 2))*
     (pow(Csc(beta), 3))*Sec(beta)))/(Pi2*v2) + 
 (0.0625*B0i(bb0, m32, MB2, MHp2)*
   (Cos(alp)*Cot(beta)*(MW2*Cos(beta) + 2*Cot(beta)*Csc(beta)*
       (pow(MB, 2))) + MW2*Cos(beta)*Sin(alp)))/(Pi2*v2) + 
 (0.03125*B0i(bb0, m32, MA02, MT2)*
   (Cos(alp)*Cot(beta)*(2*MT2*Cot(beta)*Csc(beta) + 
      MW2*Cos(beta)*(pow(Sec(w), 2))) + MW2*Cos(beta)*
     (pow(Sec(w), 2))*Sin(alp)))/(Pi2*v2) + 
 (0.00390625*Re(A0i(aa0, Mh2))*(8*MA02 + 5*M2*Cos(3*alp)*
     (pow(Csc(beta), 3)) + Cos(alp)*(5*M2 - 4*MA02 + 
      2*M2*Cos(2*(alp - 2*beta)) + (-Mh2 + MHH2)*Cos(4*alp - 2*beta) + 
      12*M2*Cos(2*beta) + 4*MA02*Cos(2*beta) - 11*Mh2*Cos(2*beta) - 
      MHH2*Cos(2*beta))*(pow(Csc(beta), 3)) - 
    3*Mh2*Csc(alp)*(pow(Csc(beta), 3))*Sin(4*alp)))/(MA02*Pi2*v2) - 
 (0.03125*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, MB2, MW2, MW2)*
   (1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.015625*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, MT2, MZ2, MZ2)*
   (pow(Sec(w), 2))*(1 + Sin(alp - beta)))/(Pi2*v2) - 
 (0.03125*(Mh2 - MHH2)*A0i(aa0, MW2)*Cos(alp)*Cot(beta)*Csc(beta)*
   Sin(2*(alp - beta)))/(MA02*Pi2*v2) - 
 (0.046875*(Mh2 - MHH2)*A0i(aa0, MZ2)*Cos(alp)*Cot(beta)*Csc(beta)*
   Sin(2*(alp - beta)))/(MA02*Pi2*v2) - 
 (0.03125*MW2*B0i(bb1, MA02, Mh2, MZ2)*Cos(alp)*Cot(beta)*Csc(beta)*
   (pow(Sec(w), 2))*Sin(2*(alp - beta)))/(Pi2*v2) + 
 (0.03125*MW2*B0i(bb1, MA02, MHH2, MZ2)*Cos(alp)*Cot(beta)*Csc(beta)*
   (pow(Sec(w), 2))*Sin(2*(alp - beta)))/(Pi2*v2) + 
 (0.0078125*Mh2*B0i(bb0, 0, Mh2, MZ2)*Cos(alp)*(MA02 - Mh2 + 2*MW2 + 
    (MA02 - Mh2)*Cos(2*w))*Cot(beta)*Csc(beta)*(pow(Sec(w), 2))*
   Sin(2*(alp - beta)))/(MA02*Pi2*v2) - 
 (0.0078125*B0i(bb0, MA02, MHH2, MZ2)*Cos(alp)*(-(MHH2*(MHH2 - 2*MW2)) + 
    MA02*(MHH2 + 2*MW2) + (MA02 - MHH2)*MHH2*Cos(2*w))*Cot(beta)*Csc(beta)*
   (pow(Sec(w), 2))*Sin(2*(alp - beta)))/(MA02*Pi2*v2) + 
 (0.0078125*MHH2*B0i(bb0, 0, MHH2, MZ2)*Cos(alp)*
   (-MA02 + MHH2 - 2*MW2 + (-MA02 + MHH2)*Cos(2*w))*Cot(beta)*Csc(beta)*
   (pow(Sec(w), 2))*Sin(2*(alp - beta)))/(MA02*Pi2*v2) + 
 (0.0078125*B0i(bb0, MA02, Mh2, MZ2)*Cos(alp)*Cot(beta)*Csc(beta)*
   (2*Mh2*MW2 + MA02*(Mh2 + 2*MW2) + (MA02 - Mh2)*Mh2*Cos(2*w) - 
    (pow(Mh, 4)))*(pow(Sec(w), 2))*Sin(2*(alp - beta)))/
  (MA02*Pi2*v2) + (0.015625*Csc(beta)*((MA02 - Mh2)*(MA02 - MHH2) - 
    (MA02 + MHH2)*MW2*(pow(Sec(w), 2)))*Re(B0i(bb0, MHH2, MA02, MZ2))*
   Sin(alp)*Sin(2*(alp - beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.03125*(Mh2*(MHH2 - MHp2) + MHp2*(MHp2 - MW2) - MHH2*(MHp2 + MW2))*
   Csc(beta)*Re(B0i(bb0, MHH2, MHp2, MW2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - 
 (0.015625*Csc(beta)*(Mh2*MHH2 - 2*MHH2*MW2 + 12*(pow(MW, 4)))*
   Re(B0i(bb0, MHH2, MW2, MW2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.0009765625*Csc(beta)*
   (3*Mh2*MHH2 - 8*MHH2*MW2 + 4*MHH2*(Mh2 - 2*MW2)*Cos(2*w) + 
    Mh2*MHH2*Cos(4*w) + 96*(pow(MW, 4)))*(pow(Sec(w), 4))*
   Re(B0i(bb0, MHH2, MZ2, MZ2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.03125*MHH2*MW2*Csc(beta)*(pow(Sec(w), 2))*
   Re(B0i(bb1, MHH2, MA02, MZ2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.0625*MHH2*MW2*Csc(beta)*
   Re(B0i(bb1, MHH2, MHp2, MW2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.0625*MHH2*MW2*Csc(beta)*
   Re(B0i(bb1, MHH2, MW2, MW2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) - (0.03125*MHH2*MW2*Csc(beta)*(pow(Sec(w), 2))*
   Re(B0i(bb1, MHH2, MZ2, MZ2))*Sin(alp)*Sin(2*(alp - beta)))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.03125*Csc(beta)*Re(B0i(bb0, Mh2, MHp2, MW2))*
   (2*(Mh2 - MHH2)*MW2*Cos(alp)*(pow(Cos(alp - beta), 2)) + 
    (MHp2*(MHH2 - MHp2 + MW2) + Mh2*(-MHH2 + MHp2 + MW2))*Sin(alp)*
     Sin(2*(alp - beta))))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.015625*MW2*Csc(beta)*Re(B0i(bb0, Mh2, MA02, MZ2))*
   (-2*Cos(alp)*(pow(Cos(alp - beta), 2))*(pow(Sec(w), 2)) + 
    (pow(Mh2 - MHH2, -1))*(pow(MW, -2))*((MA02 - Mh2)*(MA02 - MHH2) - 
      (MA02 + Mh2)*MW2*(pow(Sec(w), 2)))*Sin(alp)*Sin(2*(alp - beta))))/
  (Pi2*v2) + (0.015625*MW2*Re(B0i(bb0, Mh2, MZ2, MZ2))*
   (-8*(pow(Csc(2*w), 2)) - Csc(beta)*(pow(Csc(w), 2))*
     (pow(Sec(w), 2))*(-2*Cos(alp)*(pow(Sin(alp - beta), 2)) - 
      ((pow(Mh2 - MHH2, -1))*(pow(MW, -2))*(3*Mh2*MHH2 - 8*Mh2*MW2 + 
         4*Mh2*(MHH2 - 2*MW2)*Cos(2*w) + Mh2*MHH2*Cos(4*w) + 
         96*(pow(MW, 4)))*(pow(Sec(w), 2))*Sin(alp)*
        Sin(2*(alp - beta)))/16)))/(Pi2*v2*(pow(Csc(w), 2))) + 
 (0.0625*MW2*Re(B0i(bb1, Mh2, MW2, MW2))*(2*(Mh2 - MHH2) + 
    Sin(alp - beta)*(2*Mh2 - MHH2 + MHH2*Csc(beta)*Sin(2*alp - beta))))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.03125*MW2*(pow(Sec(w), 2))*
   Re(B0i(bb1, Mh2, MZ2, MZ2))*(2*(Mh2 - MHH2) + 
    Sin(alp - beta)*(2*Mh2 - MHH2 + MHH2*Csc(beta)*Sin(2*alp - beta))))/
  ((Mh2 - MHH2)*Pi2*v2) + (0.015625*Re(B0i(bb0, Mh2, MW2, MW2))*
   (-4*MW2 + (pow(Mh2 - MHH2, -1))*Sin(alp - beta)*
     (Mh2*(MHH2 - 4*MW2) + 2*MW2*(MHH2 + 6*MW2) + 
      Csc(beta)*(Mh2*MHH2 - 2*MHH2*MW2 + 12*(pow(MW, 4)))*
       Sin(2*alp - beta))))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ME2, ME2))*(Cos(alp) - Sin(beta)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ML2, ML2))*(Cos(alp) - Sin(beta)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MM2, MM2))*(Cos(alp) - Sin(beta)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, MC2, MC2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, ME2, ME2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, ML2, ML2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, MM2, MM2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*
   Csc(beta)*(pow(Csc(w), 2))*Re(B0i(bb1, MZ2, MT2, MT2))*
   (Cos(alp) - Sin(beta)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Csc(beta)*
   (pow(Csc(w), 2))*Re(B0i(bb1, MZ2, MU2, MU2))*(Cos(alp) - Sin(beta)))/
  (Pi2*v2) + (0.03125*Cos(alp)*(-(pow(MA02 - Mh2, 2)) + 
    (MA02 + Mh2)*MW2*(pow(Sec(w), 2)))*
   (pow(Cos(alp)*Cot(beta) + Sin(alp), 2))*Re(B0i(dbb0, Mh2, MA02, MZ2))*
   Sin(beta))/(Pi2*v2) + 
 (0.0078125*Csc(beta)*(-2*Mh2*MW2 + (pow(Mh, 4)) + 12*(pow(MW, 4)))*
   Re(B0i(dbb0, Mh2, MW2, MW2))*(-2*Cos(alp) + Cos(alp - 2*beta) + 
    Cos(3*alp - 2*beta) + 4*Sin(beta)))/(Pi2*v2) + 
 (0.00048828125*Csc(beta)*(-8*Mh2*MW2 + 4*Mh2*(Mh2 - 2*MW2)*Cos(2*w) + 
    3*(pow(Mh, 4)) + Cos(4*w)*(pow(Mh, 4)) + 96*(pow(MW, 4)))*
   (pow(Sec(w), 4))*Re(B0i(dbb0, Mh2, MZ2, MZ2))*
   (-2*Cos(alp) + Cos(alp - 2*beta) + Cos(3*alp - 2*beta) + 4*Sin(beta)))/
  (Pi2*v2) + (0.0078125*Re(A0i(aa0, MZ2))*
   (2*MA02*(-13 + 24*Cos(alp)*Csc(beta)*(pow(Cos(w), 2)))*
     (pow(Cot(w), 2)) - 2*MA02*(-13 + (5 + 6*Cos(4*w))*
       (pow(Csc(w), 2))) + Cos(alp)*(pow(Csc(beta), 2))*
     ((Mh2 - MHH2)*Sin(2*alp - 3*beta) + (Mh2 - MHH2)*Sin(2*alp - beta) - 
      4*MA02*(7 + 6*Cos(2*w))*Sin(beta))))/(MA02*Pi2*v2) + 
 (0.03125*MT2*C0i(cc0, m12, m32, m22, MHH2, MHH2, MT2)*
   (pow(Csc(beta), 3))*(pow(Sin(alp), 2))*Sec(beta)*Sin(alp - beta)*
   ((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) + 
 (0.005859375*Cos(alp - beta)*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh2)*Cos(3*alp - beta) + (2*M2 - 3*Mh2)*Cos(alp + beta))*
   (pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, Mh2, Mh2, Mh2))*Sin(alp)*((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.005859375*Cos(alp - beta)*(M2*Cos(alp - 3*beta) + 
    (M2 - Mh2)*Cos(3*alp - beta) + (2*M2 - 3*Mh2)*Cos(alp + beta))*
   (pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, Mh2, Mh2))*Sin(alp)*((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + 
    M2*Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.015625*MT2*C0i(cc0, m12, m22, m32, Mh2, MHH2, MT2)*Cos(alp - beta)*
   (pow(Csc(beta), 3))*Sec(beta)*Sin(2*alp)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) - 
 (0.015625*MT2*C0i(cc0, m12, m32, m22, Mh2, MHH2, MT2)*Cos(alp - beta)*
   (pow(Csc(beta), 3))*Sec(beta)*Sin(2*alp)*
   ((-3*M2 + 2*Mh2 + MHH2)*Sin(2*alp) + M2*Sin(2*beta)))/(Pi2*v2) - 
 (0.001953125*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, Mh2, Mh2, MHH2))*Sin(alp)*Sin(2*(alp - beta))*
   (-9*M2*Mh2 - 9*M2*MHH2 + 5*Mh2*MHH2 + 8*(pow(M2, 2)) + 
    2*(pow(Mh, 4)) + 2*(pow(MHH2, 2)) - 
    Cos(4*alp)*(5*Mh2*MHH2 - 9*M2*(Mh2 + MHH2) + 9*(pow(M2, 2)) + 
      2*(pow(Mh, 4)) + 2*(pow(MHH2, 2))) + 
    (pow(M2, 2))*(pow(Cos(beta), 4)) - 6*(pow(M2, 2))*
     (pow(Cos(beta), 2))*(pow(Sin(beta), 2)) + 
    (pow(M2, 2))*(pow(Sin(beta), 4)) - 2*M2*(Mh2 - MHH2)*Sin(2*alp)*
     Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.001953125*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, Mh2, MHH2))*Sin(alp)*Sin(2*(alp - beta))*
   (-9*M2*Mh2 - 9*M2*MHH2 + 5*Mh2*MHH2 + 8*(pow(M2, 2)) + 
    2*(pow(Mh, 4)) + 2*(pow(MHH2, 2)) - 
    Cos(4*alp)*(5*Mh2*MHH2 - 9*M2*(Mh2 + MHH2) + 9*(pow(M2, 2)) + 
      2*(pow(Mh, 4)) + 2*(pow(MHH2, 2))) + 
    (pow(M2, 2))*(pow(Cos(beta), 4)) - 6*(pow(M2, 2))*
     (pow(Cos(beta), 2))*(pow(Sin(beta), 2)) + 
    (pow(M2, 2))*(pow(Sin(beta), 4)) - 2*M2*(Mh2 - MHH2)*Sin(2*alp)*
     Sin(2*beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0234375*Csc(beta)*(pow(Csc(2*beta), 2))*Re(B0i(bb0, Mh2, MHH2, MHH2))*
   Sin(alp)*Sin(alp - beta)*((-3*M2 + Mh2 + 2*MHH2)*Sin(2*alp) - 
    M2*Sin(2*beta))*(M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0234375*Csc(beta)*(pow(Csc(2*beta), 2))*
   Re(B0i(bb0, MHH2, MHH2, MHH2))*Sin(alp)*Sin(alp - beta)*
   ((3*M2 - Mh2 - 2*MHH2)*Sin(2*alp) + M2*Sin(2*beta))*
   (M2*Sin(alp - 3*beta) + (-M2 + MHH2)*Sin(3*alp - beta) + 
    (2*M2 - 3*MHH2)*Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0009765625*((2*MA02 - Mh2)*Cos(alp - 3*beta) + 
    (4*M2 - 2*MA02 - 3*Mh2)*Cos(alp + beta))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, Mh2, MA02, MA02))*Sin(alp)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.0009765625*((-2*MA02 + Mh2)*Cos(alp - 3*beta) + 
    (-4*M2 + 2*MA02 + 3*Mh2)*Cos(alp + beta))*(pow(Csc(beta), 3))*
   (pow(Sec(beta), 2))*Re(B0i(bb0, MHH2, MA02, MA02))*Sin(alp)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) - 
 (0.0078125*(MA02 - MHH2)*B0i(bb0, 0, MA02, MHH2)*Cos(alp)*
   (pow(Csc(beta), 3))*Sin(alp - beta)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(MA02*Pi2*v2) - 
 (0.0078125*(MA02 - MHH2)*B0i(bb0, MA02, MA02, MHH2)*Cos(alp)*
   (pow(Csc(beta), 3))*Sin(alp - beta)*
   ((-2*MA02 + MHH2)*Sin(alp - 3*beta) + (-4*M2 + 2*MA02 + 3*MHH2)*
     Sin(alp + beta)))/(MA02*Pi2*v2) + 
 (0.0078125*A0i(aa0, MHH2)*Cos(alp)*(pow(Csc(beta), 3))*Sin(alp - beta)*
   (-2*(M2 - MA02)*Sin(alp - 3*beta) + (Mh2 - MHH2)*Sin(3*alp - beta) + 
    (-2*M2 - 2*MA02 + Mh2 + 3*MHH2)*Sin(alp + beta)))/(MA02*Pi2*v2) - 
 (0.001953125*((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, Mh2, MHp2, MHp2))*Sin(alp)*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.001953125*((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*(pow(Csc(beta), 3))*(pow(Sec(beta), 2))*
   Re(B0i(bb0, MHH2, MHp2, MHp2))*Sin(alp)*
   ((MHH2 - 2*MHp2)*Sin(alp - 3*beta) + (-4*M2 + 3*MHH2 + 2*MHp2)*
     Sin(alp + beta)))/((Mh2 - MHH2)*Pi2*v2) + 
 (0.00390625*(pow(Csc(w), 2))*Re(A0i(aa0, MW2))*
   (44 - 2*Cos(2*w) + 60*Cos(4*w) - 6*Cos(6*w) + Cos(alp)*(pow(MA02, -1))*
     (pow(Csc(beta), 2))*(2*(Mh2 - MHH2)*Sin(2*alp - 3*beta) + 
      2*(Mh2 - MHH2)*Sin(2*alp - beta) - 44*MA02*Sin(beta) + 
      3*MA02*Sin(beta - 6*w) - 30*MA02*Sin(beta - 4*w) - 
      Mh2*Sin(2*alp - 3*beta - 2*w) + MHH2*Sin(2*alp - 3*beta - 2*w) - 
      Mh2*Sin(2*alp - beta - 2*w) + MHH2*Sin(2*alp - beta - 2*w) + 
      MA02*Sin(beta - 2*w) - Mh2*Sin(2*alp - 3*beta + 2*w) + 
      MHH2*Sin(2*alp - 3*beta + 2*w) - Mh2*Sin(2*alp - beta + 2*w) + 
      MHH2*Sin(2*alp - beta + 2*w) + MA02*Sin(beta + 2*w) - 
      30*MA02*Sin(beta + 4*w) + 3*MA02*Sin(beta + 6*w))))/(Pi2*v2)
;
}

ComplexType dghgaga(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32)
{
return 0. + (2.6666666666666665*MT2*MW2*B0i(bb0, m32, MT2, MT2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (1.3333333333333333*(m12 + m22 - m32)*MT2*MW2*C0i(cc0, m22, m32, m12, MT2, 
    MT2, MT2)*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) - (10.666666666666666*MT2*MW2*
   C0i(cc00, m22, m32, m12, MT2, MT2, MT2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (5.333333333333333*m22*MT2*MW2*C0i(cc1, m22, m32, m12, MT2, MT2, MT2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.6666666666666665*(m12 + m22 - m32)*MT2*MW2*C0i(cc2, m22, m32, m12, MT2, 
    MT2, MT2)*(-1 + Cos(alp)*Csc(beta))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) - (0.25*MW2*B0i(bb0, m12, MHp2, MHp2)*
   ((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*Csc(2*beta)*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (1.*MW2*C0i(cc00, m22, m32, m12, MHp2, MHp2, MHp2)*
   ((Mh2 - 2*MHp2)*Cos(alp - 3*beta) + (-4*M2 + 3*Mh2 + 2*MHp2)*
     Cos(alp + beta))*Csc(2*beta)*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (0.6666666666666666*MW2*B0i(bb0, m32, MB2, MB2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (0.3333333333333333*(m12 + m22 - m32)*MW2*
   C0i(cc0, m22, m32, m12, MB2, MB2, MB2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MB, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (2.6666666666666665*MW2*C0i(cc00, m22, m32, m12, MB2, MB2, MB2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (1.3333333333333333*m22*MW2*
   C0i(cc1, m22, m32, m12, MB2, MB2, MB2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MB, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (0.6666666666666666*(m12 + m22 - m32)*MW2*C0i(cc2, m22, m32, m12, MB2, MB2, 
    MB2)*(-1 + Cos(alp)*Csc(beta))*(pow(MB, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (2.6666666666666665*MW2*B0i(bb0, m32, MC2, MC2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (1.3333333333333333*(m12 + m22 - m32)*MW2*
   C0i(cc0, m22, m32, m12, MC2, MC2, MC2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MC, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (10.666666666666666*MW2*C0i(cc00, m22, m32, m12, MC2, MC2, MC2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (5.333333333333333*m22*MW2*
   C0i(cc1, m22, m32, m12, MC2, MC2, MC2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MC, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.6666666666666665*(m12 + m22 - m32)*MW2*C0i(cc2, m22, m32, m12, MC2, MC2, 
    MC2)*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (0.6666666666666666*MW2*B0i(bb0, m32, MD2, MD2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MD, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (0.3333333333333333*(m12 + m22 - m32)*MW2*
   C0i(cc0, m22, m32, m12, MD2, MD2, MD2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MD, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (2.6666666666666665*MW2*C0i(cc00, m22, m32, m12, MD2, MD2, MD2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MD, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (1.3333333333333333*m22*MW2*
   C0i(cc1, m22, m32, m12, MD2, MD2, MD2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MD, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (0.6666666666666666*(m12 + m22 - m32)*MW2*C0i(cc2, m22, m32, m12, MD2, MD2, 
    MD2)*(-1 + Cos(alp)*Csc(beta))*(pow(MD, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (0.6666666666666666*MW2*B0i(bb0, m32, MS2, MS2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MS, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (0.3333333333333333*(m12 + m22 - m32)*MW2*
   C0i(cc0, m22, m32, m12, MS2, MS2, MS2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MS, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (2.6666666666666665*MW2*C0i(cc00, m22, m32, m12, MS2, MS2, MS2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MS, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (1.3333333333333333*m22*MW2*
   C0i(cc1, m22, m32, m12, MS2, MS2, MS2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MS, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (0.6666666666666666*(m12 + m22 - m32)*MW2*C0i(cc2, m22, m32, m12, MS2, MS2, 
    MS2)*(-1 + Cos(alp)*Csc(beta))*(pow(MS, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (2.6666666666666665*MW2*B0i(bb0, m32, MU2, MU2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MU, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (1.3333333333333333*(m12 + m22 - m32)*MW2*
   C0i(cc0, m22, m32, m12, MU2, MU2, MU2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MU, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (10.666666666666666*MW2*C0i(cc00, m22, m32, m12, MU2, MU2, MU2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MU, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (5.333333333333333*m22*MW2*
   C0i(cc1, m22, m32, m12, MU2, MU2, MU2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MU, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.6666666666666665*(m12 + m22 - m32)*MW2*C0i(cc2, m22, m32, m12, MU2, MU2, 
    MU2)*(-1 + Cos(alp)*Csc(beta))*(pow(MU, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (0.5*(6*m12 - 5*m22 - 5*m32 + Mh2)*
   C0i(cc0, m22, m32, m12, MW2, MW2, MW2)*(pow(MW, 4))*
   (pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (0.5*(m12 + m22 - m32)*C0i(cc1, m22, m32, m12, MW2, MW2, MW2)*
   (pow(MW, 4))*(pow(Cos((alp - beta)/2) + Sin((alp - beta)/2), 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (2.*MW2*B0i(bb0, m32, ME2, ME2)*(pow(ME, 2))*(pow(Sin(w), 2))*
   Sec(beta)*(Cos(beta) + Sin(alp)))/(Pi2*(pow(v, 3))) - 
 (1.*(m12 + m22 - m32)*MW2*C0i(cc0, m22, m32, m12, ME2, ME2, ME2)*
   (pow(ME, 2))*(pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*(pow(v, 3))) + (8.*MW2*C0i(cc00, m22, m32, m12, ME2, ME2, ME2)*
   (pow(ME, 2))*(pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*(pow(v, 3))) - (4.*m22*MW2*C0i(cc1, m22, m32, m12, ME2, ME2, ME2)*
   (pow(ME, 2))*(pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*(pow(v, 3))) - (2.*(m12 + m22 - m32)*MW2*
   C0i(cc2, m22, m32, m12, ME2, ME2, ME2)*(pow(ME, 2))*
   (pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*(pow(v, 3))) - (2.*MW2*B0i(bb0, m32, ML2, ML2)*(pow(ML, 2))*
   (pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*(pow(v, 3))) - (1.*(m12 + m22 - m32)*MW2*
   C0i(cc0, m22, m32, m12, ML2, ML2, ML2)*(pow(ML, 2))*
   (pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*(pow(v, 3))) + (8.*MW2*C0i(cc00, m22, m32, m12, ML2, ML2, ML2)*
   (pow(ML, 2))*(pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*(pow(v, 3))) - (4.*m22*MW2*C0i(cc1, m22, m32, m12, ML2, ML2, ML2)*
   (pow(ML, 2))*(pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*(pow(v, 3))) - (2.*(m12 + m22 - m32)*MW2*
   C0i(cc2, m22, m32, m12, ML2, ML2, ML2)*(pow(ML, 2))*
   (pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*(pow(v, 3))) - (2.*MW2*B0i(bb0, m32, MM2, MM2)*(pow(MM, 2))*
   (pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*(pow(v, 3))) - (1.*(m12 + m22 - m32)*MW2*
   C0i(cc0, m22, m32, m12, MM2, MM2, MM2)*(pow(MM, 2))*
   (pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*(pow(v, 3))) + (8.*MW2*C0i(cc00, m22, m32, m12, MM2, MM2, MM2)*
   (pow(MM, 2))*(pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*(pow(v, 3))) - (4.*m22*MW2*C0i(cc1, m22, m32, m12, MM2, MM2, MM2)*
   (pow(MM, 2))*(pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*(pow(v, 3))) - (2.*(m12 + m22 - m32)*MW2*
   C0i(cc2, m22, m32, m12, MM2, MM2, MM2)*(pow(MM, 2))*
   (pow(Sin(w), 2))*Sec(beta)*(Cos(beta) + Sin(alp)))/
  (Pi2*(pow(v, 3))) + (0.5*MW2*(Mh2 + 6*MW2)*B0i(bb0, m12, MW2, MW2)*
   (pow(Sin(w), 2))*(1 + Sin(alp - beta)))/(Pi2*(pow(v, 3))) - 
 (2.*MW2*(Mh2 + 6*MW2)*C0i(cc00, m22, m32, m12, MW2, MW2, MW2)*
   (pow(Sin(w), 2))*(1 + Sin(alp - beta)))/(Pi2*(pow(v, 3))) + 
 (0.5*B0i(bb0, m22, MW2, MW2)*(pow(MW, 4))*(pow(Sin(w), 2))*
   (1 + Sin(alp - beta)))/(Pi2*(pow(v, 3))) - 
 (0.5*B0i(bb0, m32, MW2, MW2)*(pow(MW, 4))*(pow(Sin(w), 2))*
   (1 + Sin(alp - beta)))/(Pi2*(pow(v, 3))) - 
 (1.*m12*C0i(cc2, m22, m32, m12, MW2, MW2, MW2)*(pow(MW, 4))*
   (pow(Sin(w), 2))*(1 + Sin(alp - beta)))/(Pi2*(pow(v, 3)))
;
}

ComplexType dghgg(double MHH2, double MA02, double MHp2, double M2, double beta, double alp, double m12, double m22, double m32)
{
return 0. + (0.6366197723675814*Alfashgg*MT2*B0i(bb0, m32, MT2, MT2)*
   (-1 + Cos(alp)*Csc(beta)))/v + (0.3183098861837907*Alfashgg*(m12 + m22 - m32)*
   MT2*C0i(cc0, m22, m32, m12, MT2, MT2, MT2)*(-1 + Cos(alp)*Csc(beta)))/v - 
 (2.5464790894703255*Alfashgg*MT2*C0i(cc00, m22, m32, m12, MT2, MT2, MT2)*
   (-1 + Cos(alp)*Csc(beta)))/v + 
 (1.2732395447351628*Alfashgg*m22*MT2*C0i(cc1, m22, m32, m12, MT2, MT2, MT2)*
   (-1 + Cos(alp)*Csc(beta)))/v + (0.6366197723675814*Alfashgg*(m12 + m22 - m32)*
   MT2*C0i(cc2, m22, m32, m12, MT2, MT2, MT2)*(-1 + Cos(alp)*Csc(beta)))/v + 
 (0.6366197723675814*Alfashgg*B0i(bb0, m32, MB2, MB2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MB, 2)))/v + (0.3183098861837907*Alfashgg*(m12 + m22 - m32)*
   C0i(cc0, m22, m32, m12, MB2, MB2, MB2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MB, 2)))/v - (2.5464790894703255*Alfashgg*
   C0i(cc00, m22, m32, m12, MB2, MB2, MB2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MB, 2)))/v + (1.2732395447351628*Alfashgg*m22*
   C0i(cc1, m22, m32, m12, MB2, MB2, MB2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MB, 2)))/v + (0.6366197723675814*Alfashgg*(m12 + m22 - m32)*
   C0i(cc2, m22, m32, m12, MB2, MB2, MB2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MB, 2)))/v + (0.6366197723675814*Alfashgg*B0i(bb0, m32, MC2, MC2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MC, 2)))/v + 
 (0.3183098861837907*Alfashgg*(m12 + m22 - m32)*C0i(cc0, m22, m32, m12, MC2, 
    MC2, MC2)*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2)))/v - 
 (2.5464790894703255*Alfashgg*C0i(cc00, m22, m32, m12, MC2, MC2, MC2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MC, 2)))/v + 
 (1.2732395447351628*Alfashgg*m22*C0i(cc1, m22, m32, m12, MC2, MC2, MC2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MC, 2)))/v + 
 (0.6366197723675814*Alfashgg*(m12 + m22 - m32)*C0i(cc2, m22, m32, m12, MC2, 
    MC2, MC2)*(-1 + Cos(alp)*Csc(beta))*(pow(MC, 2)))/v + 
 (0.6366197723675814*Alfashgg*B0i(bb0, m32, MD2, MD2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MD, 2)))/v + (0.3183098861837907*Alfashgg*(m12 + m22 - m32)*
   C0i(cc0, m22, m32, m12, MD2, MD2, MD2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MD, 2)))/v - (2.5464790894703255*Alfashgg*
   C0i(cc00, m22, m32, m12, MD2, MD2, MD2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MD, 2)))/v + (1.2732395447351628*Alfashgg*m22*
   C0i(cc1, m22, m32, m12, MD2, MD2, MD2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MD, 2)))/v + (0.6366197723675814*Alfashgg*(m12 + m22 - m32)*
   C0i(cc2, m22, m32, m12, MD2, MD2, MD2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MD, 2)))/v + (0.6366197723675814*Alfashgg*B0i(bb0, m32, MS2, MS2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MS, 2)))/v + 
 (0.3183098861837907*Alfashgg*(m12 + m22 - m32)*C0i(cc0, m22, m32, m12, MS2, 
    MS2, MS2)*(-1 + Cos(alp)*Csc(beta))*(pow(MS, 2)))/v - 
 (2.5464790894703255*Alfashgg*C0i(cc00, m22, m32, m12, MS2, MS2, MS2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MS, 2)))/v + 
 (1.2732395447351628*Alfashgg*m22*C0i(cc1, m22, m32, m12, MS2, MS2, MS2)*
   (-1 + Cos(alp)*Csc(beta))*(pow(MS, 2)))/v + 
 (0.6366197723675814*Alfashgg*(m12 + m22 - m32)*C0i(cc2, m22, m32, m12, MS2, 
    MS2, MS2)*(-1 + Cos(alp)*Csc(beta))*(pow(MS, 2)))/v + 
 (0.6366197723675814*Alfashgg*B0i(bb0, m32, MU2, MU2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MU, 2)))/v + (0.3183098861837907*Alfashgg*(m12 + m22 - m32)*
   C0i(cc0, m22, m32, m12, MU2, MU2, MU2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MU, 2)))/v - (2.5464790894703255*Alfashgg*
   C0i(cc00, m22, m32, m12, MU2, MU2, MU2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MU, 2)))/v + (1.2732395447351628*Alfashgg*m22*
   C0i(cc1, m22, m32, m12, MU2, MU2, MU2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MU, 2)))/v + (0.6366197723675814*Alfashgg*(m12 + m22 - m32)*
   C0i(cc2, m22, m32, m12, MU2, MU2, MU2)*(-1 + Cos(alp)*Csc(beta))*
   (pow(MU, 2)))/v
;
}

} //end namespace TypeLS
