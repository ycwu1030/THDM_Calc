#include "CouplingFunctionSM.h"

namespace SM{
ComplexType dKappahWW(double m12, double m22, double m32)
{
return 1. - (0.09375*Mh2*B0i(bb0, m12, Mh2, Mh2))/(Pi2*v2) + 
 (0.0625*(-Mh2 - 12*MW2)*B0i(bb0, m12, MW2, MW2))/(Pi2*v2) + 
 (0.03125*(-Mh2 - 24*MW2)*B0i(bb0, m12, MZ2, MZ2))/(Pi2*v2) - 
 (0.125*MW2*B0i(bb0, m22, Mh2, MW2))/(Pi2*v2) - 
 (0.125*MW2*B0i(bb0, m32, Mh2, MW2))/(Pi2*v2) - 
 (0.375*Mh2*MW2*C0i(cc0, m12, m32, m22, Mh2, Mh2, MW2))/(Pi2*v2) + 
 (0.75*MT2*(-m22 + MT2)*C0i(cc0, m22, m12, m32, MB2, MT2, MT2))/(Pi2*v2) + 
 (0.375*Mh2*C0i(cc00, m12, m32, m22, Mh2, Mh2, MW2))/(Pi2*v2) - 
 (1.5*MT2*C0i(cc00, m22, m12, m32, MB2, MT2, MT2))/(Pi2*v2) + 
 (0.125*(Mh2 + 2*MW2)*C0i(cc00, m22, m12, m32, Mh2, MW2, MW2))/(Pi2*v2) - 
 (0.1875*(m12 + 5*m22 - m32)*MT2*C0i(cc1, m22, m12, m32, MB2, MT2, MT2))/
  (Pi2*v2) - (0.125*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, MW2, MZ2, 
    MZ2))/(Pi2*v2) + (0.1875*(3*m12 - 3*m22 - m32)*MT2*
   C0i(cc2, m22, m12, m32, MB2, MT2, MT2))/(Pi2*v2) + 
 (0.125*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, MW2, MZ2, MZ2))/
  (Pi2*v2) + (0.125*MW2*C0i(cc1, m12, m32, m22, MW2, MW2, MZ2)*
   (m12 + (3*m12 + 2*m22 - 2*m32)*Cos(2*w)))/(Pi2*v2) + 
 (0.0625*MW2*C0i(cc2, m12, m32, m22, MW2, MW2, MZ2)*
   (m12 + m22 - m32 + (m12 + 9*m22 - m32)*Cos(2*w)))/(Pi2*v2) + 
 (0.125*C0i(cc00, m12, m32, m22, MW2, MW2, MZ2)*
   (Mh2 + 10*MW2 + 8*MW2*Cos(2*w)))/(Pi2*v2) - 
 (1.5*C0i(cc00, m12, m32, m22, MB2, MB2, MT2)*(pow(MB, 2)))/(Pi2*v2) + 
 (0.375*(2*m12 + m22 - m32)*C0i(cc1, m12, m32, m22, MB2, MB2, MT2)*
   (pow(MB, 2)))/(Pi2*v2) + 
 (0.1875*(m12 + 5*m22 - m32)*C0i(cc2, m12, m32, m22, MB2, MB2, MT2)*
   (pow(MB, 2)))/(Pi2*v2) + 
 (0.75*B0i(bb0, m32, MB2, MT2)*(MT2 + (pow(MB, 2))))/(Pi2*v2) + 
 (0.1875*C0i(cc0, m12, m32, m22, MB2, MB2, MT2)*(pow(MB, 2))*
   (m12 + m22 - m32 + 4*(pow(MB, 2))))/(Pi2*v2) - 
 (1.5*C0i(cc00, m12, m32, m22, MC2, MC2, MS2)*(pow(MC, 2)))/(Pi2*v2) + 
 (0.375*(2*m12 + m22 - m32)*C0i(cc1, m12, m32, m22, MC2, MC2, MS2)*
   (pow(MC, 2)))/(Pi2*v2) + 
 (0.1875*(m12 + 5*m22 - m32)*C0i(cc2, m12, m32, m22, MC2, MC2, MS2)*
   (pow(MC, 2)))/(Pi2*v2) + 
 (0.1875*C0i(cc0, m12, m32, m22, MC2, MC2, MS2)*(pow(MC, 2))*
   (m12 + m22 - m32 + 4*(pow(MC, 2))))/(Pi2*v2) - 
 (1.5*C0i(cc00, m12, m32, m22, MD2, MD2, MU2)*(pow(MD, 2)))/(Pi2*v2) + 
 (0.375*(2*m12 + m22 - m32)*C0i(cc1, m12, m32, m22, MD2, MD2, MU2)*
   (pow(MD, 2)))/(Pi2*v2) + 
 (0.1875*(m12 + 5*m22 - m32)*C0i(cc2, m12, m32, m22, MD2, MD2, MU2)*
   (pow(MD, 2)))/(Pi2*v2) + 
 (0.1875*C0i(cc0, m12, m32, m22, MD2, MD2, MU2)*(pow(MD, 2))*
   (m12 + m22 - m32 + 4*(pow(MD, 2))))/(Pi2*v2) + 
 (0.25*B0i(bb0, m32, 0, ME2)*(pow(ME, 2)))/(Pi2*v2) - 
 (0.5*C0i(cc00, m22, m12, m32, 0, ME2, ME2)*(pow(ME, 2)))/(Pi2*v2) - 
 (0.0625*(m12 + 5*m22 - m32)*C0i(cc1, m22, m12, m32, 0, ME2, ME2)*
   (pow(ME, 2)))/(Pi2*v2) + 
 (0.0625*(3*m12 - 3*m22 - m32)*C0i(cc2, m22, m12, m32, 0, ME2, ME2)*
   (pow(ME, 2)))/(Pi2*v2) + (0.25*C0i(cc0, m22, m12, m32, 0, ME2, ME2)*
   (pow(ME, 2))*(-m22 + (pow(ME, 2))))/(Pi2*v2) + 
 (0.25*B0i(bb0, m32, 0, ML2)*(pow(ML, 2)))/(Pi2*v2) - 
 (0.5*C0i(cc00, m22, m12, m32, 0, ML2, ML2)*(pow(ML, 2)))/(Pi2*v2) - 
 (0.0625*(m12 + 5*m22 - m32)*C0i(cc1, m22, m12, m32, 0, ML2, ML2)*
   (pow(ML, 2)))/(Pi2*v2) + 
 (0.0625*(3*m12 - 3*m22 - m32)*C0i(cc2, m22, m12, m32, 0, ML2, ML2)*
   (pow(ML, 2)))/(Pi2*v2) + (0.25*C0i(cc0, m22, m12, m32, 0, ML2, ML2)*
   (pow(ML, 2))*(-m22 + (pow(ML, 2))))/(Pi2*v2) + 
 (0.25*B0i(bb0, m32, 0, MM2)*(pow(MM, 2)))/(Pi2*v2) - 
 (0.5*C0i(cc00, m22, m12, m32, 0, MM2, MM2)*(pow(MM, 2)))/(Pi2*v2) - 
 (0.0625*(m12 + 5*m22 - m32)*C0i(cc1, m22, m12, m32, 0, MM2, MM2)*
   (pow(MM, 2)))/(Pi2*v2) + 
 (0.0625*(3*m12 - 3*m22 - m32)*C0i(cc2, m22, m12, m32, 0, MM2, MM2)*
   (pow(MM, 2)))/(Pi2*v2) + (0.25*C0i(cc0, m22, m12, m32, 0, MM2, MM2)*
   (pow(MM, 2))*(-m22 + (pow(MM, 2))))/(Pi2*v2) - 
 (1.5*C0i(cc00, m22, m12, m32, MC2, MS2, MS2)*(pow(MS, 2)))/(Pi2*v2) - 
 (0.1875*(m12 + 5*m22 - m32)*C0i(cc1, m22, m12, m32, MC2, MS2, MS2)*
   (pow(MS, 2)))/(Pi2*v2) + 
 (0.1875*(3*m12 - 3*m22 - m32)*C0i(cc2, m22, m12, m32, MC2, MS2, MS2)*
   (pow(MS, 2)))/(Pi2*v2) + (0.75*C0i(cc0, m22, m12, m32, MC2, MS2, MS2)*
   (pow(MS, 2))*(-m22 + (pow(MS, 2))))/(Pi2*v2) + 
 (0.75*B0i(bb0, m32, MC2, MS2)*((pow(MC, 2)) + (pow(MS, 2))))/
  (Pi2*v2) - (1.5*C0i(cc00, m22, m12, m32, MD2, MU2, MU2)*(pow(MU, 2)))/
  (Pi2*v2) - (0.1875*(m12 + 5*m22 - m32)*C0i(cc1, m22, m12, m32, MD2, MU2, 
    MU2)*(pow(MU, 2)))/(Pi2*v2) + 
 (0.1875*(3*m12 - 3*m22 - m32)*C0i(cc2, m22, m12, m32, MD2, MU2, MU2)*
   (pow(MU, 2)))/(Pi2*v2) + (0.75*C0i(cc0, m22, m12, m32, MD2, MU2, MU2)*
   (pow(MU, 2))*(-m22 + (pow(MU, 2))))/(Pi2*v2) + 
 (0.75*B0i(bb0, m32, MD2, MU2)*((pow(MD, 2)) + (pow(MU, 2))))/
  (Pi2*v2) - (0.25*C0i(cc0, m22, m12, m32, Mh2, MW2, MW2)*(pow(MW, 4)))/
  (Pi2*v2) + (0.015625*MW2*B0i(bb0, m32, MW2, MZ2)*
   (23 + 36*Cos(2*w) + 5*Cos(4*w))*(pow(Sec(w), 2)))/(Pi2*v2) - 
 (0.125*MW2*B0i(bb0, m22, 0, MW2)*(pow(Sin(w), 2)))/(Pi2*v2) + 
 (0.625*MW2*B0i(bb0, m32, 0, MW2)*(pow(Sin(w), 2)))/(Pi2*v2) - 
 (0.125*(3*m12 + 2*m22 - 4*m32 + Mh2 - 6*MW2)*MW2*
   C0i(cc0, m22, m12, m32, 0, MW2, MW2)*(pow(Sin(w), 2)))/(Pi2*v2) + 
 (2.*MW2*C0i(cc00, m22, m12, m32, 0, MW2, MW2)*(pow(Sin(w), 2)))/
  (Pi2*v2) - (0.125*(m12 + 9*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, 0, MW2, 
    MW2)*(pow(Sin(w), 2)))/(Pi2*v2) + 
 (0.125*(5*m12 - 5*m22 - 3*m32)*MW2*C0i(cc2, m22, m12, m32, 0, MW2, MW2)*
   (pow(Sin(w), 2)))/(Pi2*v2) - 
 (0.125*MW2*B0i(bb0, m22, MW2, MZ2)*(pow(Sin(w), 2))*
   (pow(Tan(w), 2)))/(Pi2*v2) + 
 (0.125*C0i(cc00, m22, m12, m32, MW2, MZ2, MZ2)*
   (Mh2 + 18*MW2 + 2*MW2*(pow(Tan(w), 2))))/(Pi2*v2) - 
 (0.125*MW2*C0i(cc0, m12, m32, m22, MW2, MW2, MZ2)*(pow(Cos(w), 2))*
   (4*m12 - 6*m22 - 4*m32 - 4*MW2 + (2*m12 + m22 - m32 + 2*MW2)*
     (pow(Tan(w), 2)) + Mh2*(pow(Tan(w), 4))))/(Pi2*v2) + 
 (0.125*MW2*C0i(cc0, m22, m12, m32, MW2, MZ2, MZ2)*
   (-5*m12 + m22 + 5*m32 + 4*MW2*(pow(Sec(w), 2)) - 
    2*MW2*(pow(Tan(w), 4))))/(Pi2*v2) - 
 (0.005208333333333333*(56 - 11*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MB2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MC2)))/(Pi2*v2) - 
 (0.005208333333333333*(56 - 11*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MD2)))/(Pi2*v2) - 
 (0.015625*(-8 + 21*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ME2)))/(Pi2*v2) + (0.03125*Re(A0i(aa0, Mh2)))/(Pi2*v2) - 
 (0.015625*(-8 + 21*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ML2)))/(Pi2*v2) - 
 (0.015625*(-8 + 21*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MM2)))/(Pi2*v2) - 
 (0.005208333333333333*(56 - 11*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MS2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MT2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MU2)))/(Pi2*v2) + 
 (0.0078125*(2 - 23*Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MW2)))/(Pi2*v2) + 
 (0.015625*(13 + 11*Cos(2*w))*(pow(Csc(w), 2))*Re(A0i(aa0, MZ2)))/
  (Pi2*v2) - (1.125*MW2*(pow(Sin(w), 2))*Re(B0i(bb0, 0, MW2, MW2)))/
  (Pi2*v2) + (0.0625*MW2*Re(B0i(bb0, Mh2, MW2, MW2)))/(Pi2*v2) + 
 (0.03125*MW2*(pow(Sec(w), 2))*Re(B0i(bb0, Mh2, MZ2, MZ2)))/(Pi2*v2) + 
 (0.125*MW2*(1 - 5*Cos(2*w))*Re(B0i(bb0, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.375*(-MT2 + MW2*Cos(2*w))*(pow(Csc(w), 2))*
   Re(B0i(bb0, MW2, MB2, MT2)))/(Pi2*v2) - 
 (0.375*(pow(MC, 2))*(pow(Csc(w), 2))*Re(B0i(bb0, MW2, MC2, MS2)))/
  (Pi2*v2) + (0.375*(MW2*Cos(2*w) - (pow(MU, 2)))*(pow(Csc(w), 2))*
   Re(B0i(bb0, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.125*MW2*(pow(Csc(w), 2))*Re(B0i(bb0, MW2, Mh2, MW2)))/(Pi2*v2) - 
 (0.0078125*MW2*(22 + 45*Cos(2*w) + 10*Cos(4*w) + 3*Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Sec(w), 2))*Re(B0i(bb0, MW2, MW2, MZ2)))/
  (Pi2*v2) + (0.1875*(pow(MB, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.1875*(pow(MC, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MC2, MC2)))/
  (Pi2*v2) + (0.1875*(pow(MD, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*(pow(ME, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, ME2, ME2)))/
  (Pi2*v2) - (0.125*MW2*(pow(Csc(w), 2))*Re(B0i(bb0, MZ2, Mh2, MZ2)))/
  (Pi2*v2) + (0.0625*(pow(ML, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*(pow(MM, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MM2, MM2)))/
  (Pi2*v2) + (0.1875*(pow(MS, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.1875*MT2*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.1875*(pow(MU, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MU2, MU2)))/
  (Pi2*v2) + (0.0625*MW2*(5 + 9*Cos(2*w))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MB2, MB2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MD2, MD2)))/(Pi2*v2) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, ME2, ME2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, ML2, ML2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, MM2, MM2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MT2, MT2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.5*(2 + 3*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, MW2, MW2)))/
  (Pi2*v2) + (0.25*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, 0, ME2)))/
  (Pi2*v2) + (0.25*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, 0, ML2)))/
  (Pi2*v2) + (0.25*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, 0, MM2)))/
  (Pi2*v2) - (1.*Re(B0i(bb00, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.75*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, MB2, MT2)))/(Pi2*v2) + 
 (0.75*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.75*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, MD2, MU2)))/(Pi2*v2) - 
 (0.125*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, Mh2, MW2)))/(Pi2*v2) - 
 (0.125*(5 + 4*Cos(2*w))*(pow(Csc(w), 2))*Re(B0i(bb00, MW2, MW2, MZ2)))/
  (Pi2*v2) - (0.375*(pow(Cot(w), 2))*Re(B0i(bb00, MZ2, 0, 0)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MB2, MB2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MC2, MC2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MD2, MD2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*(pow(Cot(w), 2))*Re(B0i(bb00, MZ2, Mh2, MZ2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ML2, ML2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MM2, MM2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MS2, MS2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MT2, MT2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.0625*(7 + 8*Cos(2*w) + 3*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MW2, MW2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MB2, MB2)))/
  (Pi2*v2) - (0.6666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MC2, MC2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MD2, MD2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ME2, ME2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ML2, ML2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MM2, MM2)))/
  (Pi2*v2) - (0.16666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MS2, MS2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MT2, MT2)))/
  (Pi2*v2) - (0.6666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.25*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.375*(pow(MB, 2))*Re(B0i(bb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*(pow(MC, 2))*Re(B0i(bb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*(pow(MD, 2))*Re(B0i(bb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*(pow(ME, 2))*Re(B0i(bb1, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*(pow(ML, 2))*Re(B0i(bb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*(pow(MM, 2))*Re(B0i(bb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*(pow(MS, 2))*Re(B0i(bb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*MT2*Re(B0i(bb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*(pow(MU, 2))*Re(B0i(bb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*MW2*Re(B0i(bb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*MW2*(pow(Sec(w), 2))*Re(B0i(bb1, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ME2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ML2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, MM2)))/(Pi2*v2) - 
 (0.25*MW2*Cos(2*w)*Re(B0i(bb1, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MB2, MT2)))/
  (Pi2*v2) - (0.375*MW2*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MD2, MU2)))/
  (Pi2*v2) + (0.25*MW2*Cos(2*w)*(pow(Cot(w), 2))*
   Re(B0i(bb1, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.1875*MW2*(pow(Csc(w), 2))*Re(B0i(bb1, MZ2, 0, 0)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.25*MW2*(pow(Cos(w), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb1, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.75*(pow(MB, 4))*Re(B0i(dbb0, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.75*(pow(MC, 4))*Re(B0i(dbb0, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.75*(pow(MD, 4))*Re(B0i(dbb0, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.25*(pow(ME, 4))*Re(B0i(dbb0, Mh2, ME2, ME2)))/(Pi2*v2) - 
 (0.140625*(pow(Mh, 4))*Re(B0i(dbb0, Mh2, Mh2, Mh2)))/(Pi2*v2) + 
 (0.25*(pow(ML, 4))*Re(B0i(dbb0, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.25*(pow(MM, 4))*Re(B0i(dbb0, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.75*(pow(MS, 4))*Re(B0i(dbb0, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.75*(pow(MT, 4))*Re(B0i(dbb0, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.75*(pow(MU, 4))*Re(B0i(dbb0, Mh2, MU2, MU2)))/(Pi2*v2) + 
 (0.03125*(2*Mh2*MW2 - (pow(Mh, 4)) - 12*(pow(MW, 4)))*
   Re(B0i(dbb0, Mh2, MW2, MW2)))/(Pi2*v2) + 
 (0.015625*(-(pow(Mh, 4)) + 2*Mh2*MW2*(pow(Sec(w), 2)) - 
    12*(pow(MW, 4))*(pow(Sec(w), 4)))*Re(B0i(dbb0, Mh2, MZ2, MZ2)))/
  (Pi2*v2) + (1.*(pow(MW, 4))*(pow(Sin(w), 2))*
   Re(B0i(dbb0, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.75*(MT2 - MW2)*MW2*Re(B0i(dbb0, MW2, MB2, MT2)))/(Pi2*v2) + 
 (0.75*MW2*(pow(MC, 2))*Re(B0i(dbb0, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.75*MW2*(-MW2 + (pow(MU, 2)))*Re(B0i(dbb0, MW2, MD2, MU2)))/(Pi2*v2) - 
 (0.25*(pow(MW, 4))*Re(B0i(dbb0, MW2, Mh2, MW2)))/(Pi2*v2) + 
 (0.0625*(7 + 12*Cos(2*w) + Cos(4*w))*(pow(MW, 4))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MB2, MB2)))/
  (Pi2*v2) + (1.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MD2, MD2)))/
  (Pi2*v2) + (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ME2, ME2)))/
  (Pi2*v2) + (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ML2, ML2)))/
  (Pi2*v2) + (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MM2, MM2)))/
  (Pi2*v2) + (0.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MT2, MT2)))/
  (Pi2*v2) + (1.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (1.5*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MW2, MW2)))/(Pi2*v2) - 
 (0.5*MW2*Re(B0i(dbb00, MW2, 0, ME2)))/(Pi2*v2) - 
 (0.5*MW2*Re(B0i(dbb00, MW2, 0, ML2)))/(Pi2*v2) - 
 (0.5*MW2*Re(B0i(dbb00, MW2, 0, MM2)))/(Pi2*v2) + 
 (2.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, MW2, 0, MW2)))/(Pi2*v2) - 
 (1.5*MW2*Re(B0i(dbb00, MW2, MB2, MT2)))/(Pi2*v2) - 
 (1.5*MW2*Re(B0i(dbb00, MW2, MC2, MS2)))/(Pi2*v2) - 
 (1.5*MW2*Re(B0i(dbb00, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.25*MW2*Re(B0i(dbb00, MW2, Mh2, MW2)))/(Pi2*v2) + 
 (0.25*MW2*(5 + 4*Cos(2*w))*Re(B0i(dbb00, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MB, 2))*Re(B0i(dbb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MC, 2))*Re(B0i(dbb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MD, 2))*Re(B0i(dbb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ME, 2))*Re(B0i(dbb1, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ML, 2))*Re(B0i(dbb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(MM, 2))*Re(B0i(dbb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MS, 2))*Re(B0i(dbb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*Mh2*MT2*Re(B0i(dbb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MU, 2))*Re(B0i(dbb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*Mh2*MW2*Re(B0i(dbb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*Mh2*MW2*(pow(Sec(w), 2))*Re(B0i(dbb1, Mh2, MZ2, MZ2)))/
  (Pi2*v2) + (0.25*(pow(MW, 4))*Re(B0i(dbb1, MW2, 0, ME2)))/(Pi2*v2) + 
 (0.25*(pow(MW, 4))*Re(B0i(dbb1, MW2, 0, ML2)))/(Pi2*v2) + 
 (0.25*(pow(MW, 4))*Re(B0i(dbb1, MW2, 0, MM2)))/(Pi2*v2) + 
 (0.5*(pow(MW, 4))*(pow(Sin(w), 2))*Re(B0i(dbb1, MW2, 0, MW2)))/
  (Pi2*v2) - (0.75*(pow(MW, 4))*Re(B0i(dbb1, MW2, MB2, MT2)))/(Pi2*v2) + 
 (0.75*(pow(MW, 4))*Re(B0i(dbb1, MW2, MC2, MS2)))/(Pi2*v2) - 
 (0.75*(pow(MW, 4))*Re(B0i(dbb1, MW2, MD2, MU2)))/(Pi2*v2) - 
 (0.5*(pow(MW, 4))*(pow(Cos(w), 2))*Re(B0i(dbb1, MW2, MW2, MZ2)))/
  (Pi2*v2)
;
}

ComplexType dKappahZZ(double m12, double m22, double m32)
{
return 1. - (0.09375*Mh2*B0i(bb0, m12, Mh2, Mh2))/(Pi2*v2) - 
 (0.03125*Mh2*B0i(bb0, m12, MZ2, MZ2))/(Pi2*v2) + 
 (0.375*Mh2*C0i(cc00, m12, m32, m22, Mh2, Mh2, MZ2))/(Pi2*v2) + 
 (0.03125*MW2*B0i(bb0, m32, MW2, MW2)*(11 + 20*Cos(2*w) + Cos(4*w)))/
  (Pi2*v2) + (0.08333333333333333*MT2*B0i(bb0, m32, MT2, MT2)*
   (9 - 2*Cos(2*w) + 2*Cos(4*w)))/(Pi2*v2) - 
 (0.16666666666666666*MT2*C0i(cc00, m22, m32, m12, MT2, MT2, MT2)*
   (9 - 4*Cos(2*w) + 4*Cos(4*w)))/(Pi2*v2) + 
 (0.020833333333333332*MT2*C0i(cc1, m22, m32, m12, MT2, MT2, MT2)*
   (9*(m12 + 5*m22 - m32) - 16*m22*Cos(2*w) + 16*m22*Cos(4*w)))/(Pi2*v2) + 
 (0.041666666666666664*MT2*C0i(cc2, m22, m32, m12, MT2, MT2, MT2)*
   (9*(2*m12 + m22 - m32) - 4*(m12 + m22 - m32)*Cos(2*w) + 
    4*(m12 + m22 - m32)*Cos(4*w)))/(Pi2*v2) + 
 (0.020833333333333332*MT2*C0i(cc0, m22, m32, m12, MT2, MT2, MT2)*
   (9*(m12 + m22 - m32 + 4*MT2) - 4*(m12 + m22 - m32)*Cos(2*w) + 
    4*(m12 + m22 - m32)*Cos(4*w)))/(Pi2*v2) + 
 (0.125*C0i(cc00, m22, m32, m12, MW2, MW2, MW2)*
   (Mh2 + 14*MW2 + 16*MW2*Cos(2*w) + (Mh2 + 6*MW2)*Cos(4*w)))/(Pi2*v2) - 
 (0.16666666666666666*C0i(cc00, m22, m32, m12, MB2, MB2, MB2)*
   (6 + 2*Cos(2*w) + Cos(4*w))*(pow(MB, 2)))/(Pi2*v2) + 
 (0.041666666666666664*B0i(bb0, m32, MB2, MB2)*(15 + 2*Cos(2*w) + Cos(4*w))*
   (pow(MB, 2)))/(Pi2*v2) + 
 (0.020833333333333332*C0i(cc1, m22, m32, m12, MB2, MB2, MB2)*
   (9*m12 + 33*m22 - 9*m32 + 8*m22*Cos(2*w) + 4*m22*Cos(4*w))*
   (pow(MB, 2)))/(Pi2*v2) + 
 (0.041666666666666664*C0i(cc2, m22, m32, m12, MB2, MB2, MB2)*
   (15*m12 + 6*m22 - 6*m32 + 2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w))*(pow(MB, 2)))/(Pi2*v2) + 
 (0.020833333333333332*C0i(cc0, m22, m32, m12, MB2, MB2, MB2)*
   (pow(MB, 2))*(2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w) + 6*(m12 + m22 - m32 + 6*(pow(MB, 2)))))/
  (Pi2*v2) + (0.08333333333333333*B0i(bb0, m32, MC2, MC2)*
   (9 - 2*Cos(2*w) + 2*Cos(4*w))*(pow(MC, 2)))/(Pi2*v2) - 
 (0.16666666666666666*C0i(cc00, m22, m32, m12, MC2, MC2, MC2)*
   (9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(MC, 2)))/(Pi2*v2) + 
 (0.020833333333333332*C0i(cc1, m22, m32, m12, MC2, MC2, MC2)*
   (9*(m12 + 5*m22 - m32) - 16*m22*Cos(2*w) + 16*m22*Cos(4*w))*
   (pow(MC, 2)))/(Pi2*v2) + 
 (0.041666666666666664*C0i(cc2, m22, m32, m12, MC2, MC2, MC2)*
   (9*(2*m12 + m22 - m32) - 4*(m12 + m22 - m32)*Cos(2*w) + 
    4*(m12 + m22 - m32)*Cos(4*w))*(pow(MC, 2)))/(Pi2*v2) + 
 (0.020833333333333332*C0i(cc0, m22, m32, m12, MC2, MC2, MC2)*
   (pow(MC, 2))*(-4*(m12 + m22 - m32)*Cos(2*w) + 
    4*(m12 + m22 - m32)*Cos(4*w) + 9*(m12 + m22 - m32 + 4*(pow(MC, 2)))))/
  (Pi2*v2) - (0.16666666666666666*C0i(cc00, m22, m32, m12, MD2, MD2, MD2)*
   (6 + 2*Cos(2*w) + Cos(4*w))*(pow(MD, 2)))/(Pi2*v2) + 
 (0.041666666666666664*B0i(bb0, m32, MD2, MD2)*(15 + 2*Cos(2*w) + Cos(4*w))*
   (pow(MD, 2)))/(Pi2*v2) + 
 (0.020833333333333332*C0i(cc1, m22, m32, m12, MD2, MD2, MD2)*
   (9*m12 + 33*m22 - 9*m32 + 8*m22*Cos(2*w) + 4*m22*Cos(4*w))*
   (pow(MD, 2)))/(Pi2*v2) + 
 (0.041666666666666664*C0i(cc2, m22, m32, m12, MD2, MD2, MD2)*
   (15*m12 + 6*m22 - 6*m32 + 2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w))*(pow(MD, 2)))/(Pi2*v2) + 
 (0.020833333333333332*C0i(cc0, m22, m32, m12, MD2, MD2, MD2)*
   (pow(MD, 2))*(2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w) + 6*(m12 + m22 - m32 + 6*(pow(MD, 2)))))/
  (Pi2*v2) - (0.5*C0i(cc00, m22, m32, m12, ME2, ME2, ME2)*
   (2 - 2*Cos(2*w) + Cos(4*w))*(pow(ME, 2)))/(Pi2*v2) + 
 (0.125*B0i(bb0, m32, ME2, ME2)*(3 - 2*Cos(2*w) + Cos(4*w))*(pow(ME, 2)))/
  (Pi2*v2) + (0.0625*C0i(cc1, m22, m32, m12, ME2, ME2, ME2)*
   (m12 + 9*m22 - m32 - 8*m22*Cos(2*w) + 4*m22*Cos(4*w))*(pow(ME, 2)))/
  (Pi2*v2) + (0.125*C0i(cc2, m22, m32, m12, ME2, ME2, ME2)*
   (3*m12 + 2*m22 - 2*m32 - 2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w))*(pow(ME, 2)))/(Pi2*v2) + 
 (0.0625*C0i(cc0, m22, m32, m12, ME2, ME2, ME2)*(pow(ME, 2))*
   (-2*(m12 + m22 - m32)*Cos(2*w) + (m12 + m22 - m32)*Cos(4*w) + 
    2*(m12 + m22 - m32 + 2*(pow(ME, 2)))))/(Pi2*v2) - 
 (0.5*C0i(cc00, m22, m32, m12, ML2, ML2, ML2)*(2 - 2*Cos(2*w) + Cos(4*w))*
   (pow(ML, 2)))/(Pi2*v2) + 
 (0.125*B0i(bb0, m32, ML2, ML2)*(3 - 2*Cos(2*w) + Cos(4*w))*(pow(ML, 2)))/
  (Pi2*v2) + (0.0625*C0i(cc1, m22, m32, m12, ML2, ML2, ML2)*
   (m12 + 9*m22 - m32 - 8*m22*Cos(2*w) + 4*m22*Cos(4*w))*(pow(ML, 2)))/
  (Pi2*v2) + (0.125*C0i(cc2, m22, m32, m12, ML2, ML2, ML2)*
   (3*m12 + 2*m22 - 2*m32 - 2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w))*(pow(ML, 2)))/(Pi2*v2) + 
 (0.0625*C0i(cc0, m22, m32, m12, ML2, ML2, ML2)*(pow(ML, 2))*
   (-2*(m12 + m22 - m32)*Cos(2*w) + (m12 + m22 - m32)*Cos(4*w) + 
    2*(m12 + m22 - m32 + 2*(pow(ML, 2)))))/(Pi2*v2) - 
 (0.5*C0i(cc00, m22, m32, m12, MM2, MM2, MM2)*(2 - 2*Cos(2*w) + Cos(4*w))*
   (pow(MM, 2)))/(Pi2*v2) + 
 (0.125*B0i(bb0, m32, MM2, MM2)*(3 - 2*Cos(2*w) + Cos(4*w))*(pow(MM, 2)))/
  (Pi2*v2) + (0.0625*C0i(cc1, m22, m32, m12, MM2, MM2, MM2)*
   (m12 + 9*m22 - m32 - 8*m22*Cos(2*w) + 4*m22*Cos(4*w))*(pow(MM, 2)))/
  (Pi2*v2) + (0.125*C0i(cc2, m22, m32, m12, MM2, MM2, MM2)*
   (3*m12 + 2*m22 - 2*m32 - 2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w))*(pow(MM, 2)))/(Pi2*v2) + 
 (0.0625*C0i(cc0, m22, m32, m12, MM2, MM2, MM2)*(pow(MM, 2))*
   (-2*(m12 + m22 - m32)*Cos(2*w) + (m12 + m22 - m32)*Cos(4*w) + 
    2*(m12 + m22 - m32 + 2*(pow(MM, 2)))))/(Pi2*v2) - 
 (0.16666666666666666*C0i(cc00, m22, m32, m12, MS2, MS2, MS2)*
   (6 + 2*Cos(2*w) + Cos(4*w))*(pow(MS, 2)))/(Pi2*v2) + 
 (0.041666666666666664*B0i(bb0, m32, MS2, MS2)*(15 + 2*Cos(2*w) + Cos(4*w))*
   (pow(MS, 2)))/(Pi2*v2) + 
 (0.020833333333333332*C0i(cc1, m22, m32, m12, MS2, MS2, MS2)*
   (9*m12 + 33*m22 - 9*m32 + 8*m22*Cos(2*w) + 4*m22*Cos(4*w))*
   (pow(MS, 2)))/(Pi2*v2) + 
 (0.041666666666666664*C0i(cc2, m22, m32, m12, MS2, MS2, MS2)*
   (15*m12 + 6*m22 - 6*m32 + 2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w))*(pow(MS, 2)))/(Pi2*v2) + 
 (0.020833333333333332*C0i(cc0, m22, m32, m12, MS2, MS2, MS2)*
   (pow(MS, 2))*(2*(m12 + m22 - m32)*Cos(2*w) + 
    (m12 + m22 - m32)*Cos(4*w) + 6*(m12 + m22 - m32 + 6*(pow(MS, 2)))))/
  (Pi2*v2) + (0.08333333333333333*B0i(bb0, m32, MU2, MU2)*
   (9 - 2*Cos(2*w) + 2*Cos(4*w))*(pow(MU, 2)))/(Pi2*v2) - 
 (0.16666666666666666*C0i(cc00, m22, m32, m12, MU2, MU2, MU2)*
   (9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(MU, 2)))/(Pi2*v2) + 
 (0.020833333333333332*C0i(cc1, m22, m32, m12, MU2, MU2, MU2)*
   (9*(m12 + 5*m22 - m32) - 16*m22*Cos(2*w) + 16*m22*Cos(4*w))*
   (pow(MU, 2)))/(Pi2*v2) + 
 (0.041666666666666664*C0i(cc2, m22, m32, m12, MU2, MU2, MU2)*
   (9*(2*m12 + m22 - m32) - 4*(m12 + m22 - m32)*Cos(2*w) + 
    4*(m12 + m22 - m32)*Cos(4*w))*(pow(MU, 2)))/(Pi2*v2) + 
 (0.020833333333333332*C0i(cc0, m22, m32, m12, MU2, MU2, MU2)*
   (pow(MU, 2))*(-4*(m12 + m22 - m32)*Cos(2*w) + 
    4*(m12 + m22 - m32)*Cos(4*w) + 9*(m12 + m22 - m32 + 4*(pow(MU, 2)))))/
  (Pi2*v2) + (0.25*MW2*C0i(cc2, m22, m32, m12, MW2, MW2, MW2)*
   (3*m12 + 2*m22 - 2*m32 + m12*Cos(2*w))*(pow(Cos(w), 2)))/(Pi2*v2) + 
 (0.125*MW2*C0i(cc1, m22, m32, m12, MW2, MW2, MW2)*
   (m12 + 9*m22 - m32 + (m12 + m22 - m32)*Cos(2*w))*(pow(Cos(w), 2)))/
  (Pi2*v2) - (0.125*MW2*B0i(bb0, m22, Mh2, MZ2)*(pow(Sec(w), 2)))/
  (Pi2*v2) - (0.125*MW2*B0i(bb0, m32, Mh2, MZ2)*(pow(Sec(w), 2)))/
  (Pi2*v2) - (0.375*Mh2*MW2*C0i(cc0, m12, m32, m22, Mh2, Mh2, MZ2)*
   (pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.125*C0i(cc00, m22, m12, m32, Mh2, MZ2, MZ2)*
   (Mh2 + 2*MW2*(pow(Sec(w), 2))))/(Pi2*v2) - 
 (0.25*C0i(cc0, m22, m12, m32, Mh2, MZ2, MZ2)*(pow(MW, 4))*
   (pow(Sec(w), 4)))/(Pi2*v2) - 
 (0.25*MW2*B0i(bb0, m22, MW2, MW2)*(pow(Sin(w), 4)))/(Pi2*v2) - 
 (0.0625*B0i(bb0, m12, MW2, MW2)*(Mh2 - 2*Mh2*(pow(Cot(w), 2)) + 
    (Mh2 + 24*MW2)*(pow(Cot(w), 4)))*(pow(Sin(w), 4)))/(Pi2*v2) + 
 (0.005208333333333333*(-4 - 43*Cos(2*w) + 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MB2)))/(Pi2*v2) + 
 (0.010416666666666666*(7 + Cos(2*w) + 8*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MC2)))/(Pi2*v2) + 
 (0.005208333333333333*(-4 - 43*Cos(2*w) + 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MD2)))/(Pi2*v2) + 
 (0.015625*(4 - 11*Cos(2*w) + 2*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ME2)))/(Pi2*v2) + (0.03125*Re(A0i(aa0, Mh2)))/(Pi2*v2) + 
 (0.015625*(4 - 11*Cos(2*w) + 2*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ML2)))/(Pi2*v2) + 
 (0.015625*(4 - 11*Cos(2*w) + 2*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MM2)))/(Pi2*v2) + 
 (0.005208333333333333*(-4 - 43*Cos(2*w) + 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MS2)))/(Pi2*v2) + 
 (0.010416666666666666*(7 + Cos(2*w) + 8*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MT2)))/(Pi2*v2) + 
 (0.010416666666666666*(7 + Cos(2*w) + 8*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MU2)))/(Pi2*v2) - 
 (0.0078125*(10 + Cos(2*w) + 34*Cos(4*w) + 3*Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MW2)))/(Pi2*v2) + (0.015625*(7 + 11*Cos(2*w) + 6*Cos(4*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MZ2)))/(Pi2*v2) - 
 (1.125*MW2*(pow(Sin(w), 2))*Re(B0i(bb0, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.0625*MW2*Re(B0i(bb0, Mh2, MW2, MW2)))/(Pi2*v2) + 
 (0.03125*MW2*(pow(Sec(w), 2))*Re(B0i(bb0, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (0.5*MW2*Cos(2*w)*Re(B0i(bb0, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.375*(-MT2 + MW2)*Cos(2*w)*(pow(Csc(w), 2))*
   Re(B0i(bb0, MW2, MB2, MT2)))/(Pi2*v2) - 
 (0.375*(pow(MC, 2))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*Cos(2*w)*(MW2 - (pow(MU, 2)))*(pow(Csc(w), 2))*
   Re(B0i(bb0, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb0, MW2, Mh2, MW2)))/
  (Pi2*v2) - (0.015625*MW2*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Sec(w), 2))*Re(B0i(bb0, MW2, MW2, MZ2)))/
  (Pi2*v2) + (0.1875*(pow(MB, 2))*(-2 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.1875*(pow(MC, 2))*(-2 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.1875*(pow(MD, 2))*(-2 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*(pow(ME, 2))*(-2 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MZ2, ME2, ME2)))/(Pi2*v2) - 
 (0.0625*MW2*(-1 + 3*Cos(2*w))*(pow(Csc(w), 2))*(pow(Sec(w), 2))*
   Re(B0i(bb0, MZ2, Mh2, MZ2)))/(Pi2*v2) + 
 (0.0625*(pow(ML, 2))*(-2 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*(pow(MM, 2))*(-2 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.1875*(pow(MS, 2))*(-2 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.1875*MT2*(-2 + (pow(Cot(w), 2)))*Re(B0i(bb0, MZ2, MT2, MT2)))/
  (Pi2*v2) + (0.1875*(pow(MU, 2))*(-2 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.015625*MW2*(27 + 12*Cos(2*w) + 17*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb0, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MB2, MB2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MD2, MD2)))/(Pi2*v2) + 
 (0.5*(-1 + 2*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, ME2, ME2)))/
  (Pi2*v2) + (0.5*(-1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, ML2, ML2)))/(Pi2*v2) + 
 (0.5*(-1 + 2*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, MM2, MM2)))/
  (Pi2*v2) + (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MT2, MT2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.5*(2 + 3*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, MW2, MW2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, ME2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, ML2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, MM2)))/
  (Pi2*v2) - (1.*Cos(2*w)*Re(B0i(bb00, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MB2, MT2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.125*(1 - (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, Mh2, MW2)))/(Pi2*v2) - 
 (0.125*(2 + 5*Cos(2*w) + 2*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MW2, MW2, MZ2)))/(Pi2*v2) - 
 (0.375*(-2 + (pow(Cot(w), 2)))*Re(B0i(bb00, MZ2, 0, 0)))/(Pi2*v2) - 
 (0.010416666666666666*(-6 + 35*Cos(2*w) + 4*Cos(4*w) + 3*Cos(6*w))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MZ2, MB2, MB2)))/(Pi2*v2) - 
 (0.020833333333333332*(-15 + 37*Cos(2*w) - 10*Cos(4*w) + 6*Cos(6*w))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MZ2, MC2, MC2)))/(Pi2*v2) - 
 (0.010416666666666666*(-6 + 35*Cos(2*w) + 4*Cos(4*w) + 3*Cos(6*w))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MZ2, MD2, MD2)))/(Pi2*v2) - 
 (0.03125*(-10 + 19*Cos(2*w) - 8*Cos(4*w) + 3*Cos(6*w))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*(-2 + (pow(Cot(w), 2)))*Re(B0i(bb00, MZ2, Mh2, MZ2)))/(Pi2*v2) - 
 (0.03125*(-10 + 19*Cos(2*w) - 8*Cos(4*w) + 3*Cos(6*w))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MZ2, ML2, ML2)))/(Pi2*v2) - 
 (0.03125*(-10 + 19*Cos(2*w) - 8*Cos(4*w) + 3*Cos(6*w))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MZ2, MM2, MM2)))/(Pi2*v2) - 
 (0.010416666666666666*(-6 + 35*Cos(2*w) + 4*Cos(4*w) + 3*Cos(6*w))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MZ2, MS2, MS2)))/(Pi2*v2) - 
 (0.020833333333333332*(-15 + 37*Cos(2*w) - 10*Cos(4*w) + 6*Cos(6*w))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MZ2, MT2, MT2)))/(Pi2*v2) - 
 (0.020833333333333332*(-15 + 37*Cos(2*w) - 10*Cos(4*w) + 6*Cos(6*w))*
   (pow(Csc(w), 2))*Re(B0i(bb00, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.015625*(10 + 35*Cos(2*w) + 18*Cos(4*w) + 9*Cos(6*w))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MZ2, MW2, MW2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MB2, MB2)))/
  (Pi2*v2) - (0.6666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MC2, MC2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MD2, MD2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ME2, ME2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ML2, ML2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MM2, MM2)))/
  (Pi2*v2) - (0.16666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MS2, MS2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MT2, MT2)))/
  (Pi2*v2) - (0.6666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.25*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.375*(pow(MB, 2))*Re(B0i(bb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*(pow(MC, 2))*Re(B0i(bb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*(pow(MD, 2))*Re(B0i(bb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*(pow(ME, 2))*Re(B0i(bb1, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*(pow(ML, 2))*Re(B0i(bb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*(pow(MM, 2))*Re(B0i(bb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*(pow(MS, 2))*Re(B0i(bb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*MT2*Re(B0i(bb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*(pow(MU, 2))*Re(B0i(bb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*MW2*Re(B0i(bb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*MW2*(pow(Sec(w), 2))*Re(B0i(bb1, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ME2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ML2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, MM2)))/(Pi2*v2) - 
 (0.25*MW2*Cos(2*w)*Re(B0i(bb1, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MB2, MT2)))/
  (Pi2*v2) - (0.375*MW2*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MD2, MU2)))/
  (Pi2*v2) + (0.25*MW2*Cos(2*w)*(pow(Cot(w), 2))*
   Re(B0i(bb1, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.1875*MW2*(pow(Csc(w), 2))*Re(B0i(bb1, MZ2, 0, 0)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.25*MW2*(pow(Cos(w), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb1, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.75*(pow(MB, 4))*Re(B0i(dbb0, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.75*(pow(MC, 4))*Re(B0i(dbb0, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.75*(pow(MD, 4))*Re(B0i(dbb0, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.25*(pow(ME, 4))*Re(B0i(dbb0, Mh2, ME2, ME2)))/(Pi2*v2) - 
 (0.140625*(pow(Mh, 4))*Re(B0i(dbb0, Mh2, Mh2, Mh2)))/(Pi2*v2) + 
 (0.25*(pow(ML, 4))*Re(B0i(dbb0, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.25*(pow(MM, 4))*Re(B0i(dbb0, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.75*(pow(MS, 4))*Re(B0i(dbb0, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.75*(pow(MT, 4))*Re(B0i(dbb0, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.75*(pow(MU, 4))*Re(B0i(dbb0, Mh2, MU2, MU2)))/(Pi2*v2) + 
 (0.03125*(2*Mh2*MW2 - (pow(Mh, 4)) - 12*(pow(MW, 4)))*
   Re(B0i(dbb0, Mh2, MW2, MW2)))/(Pi2*v2) + 
 (0.015625*(-(pow(Mh, 4)) + 2*Mh2*MW2*(pow(Sec(w), 2)) - 
    12*(pow(MW, 4))*(pow(Sec(w), 4)))*Re(B0i(dbb0, Mh2, MZ2, MZ2)))/
  (Pi2*v2) + (0.375*MW2*(pow(MB, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*MW2*(pow(MC, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*MW2*(pow(MD, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*MW2*(pow(ME, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, ME2, ME2)))/(Pi2*v2) - 
 (0.25*(pow(MW, 4))*(pow(Sec(w), 4))*Re(B0i(dbb0, MZ2, Mh2, MZ2)))/
  (Pi2*v2) + (0.125*MW2*(pow(ML, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*MW2*(pow(MM, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*MW2*(pow(MS, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*MT2*MW2*(pow(Sec(w), 2))*Re(B0i(dbb0, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*MW2*(pow(MU, 2))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.125*(5 + 9*Cos(2*w))*(pow(MW, 4))*(pow(Sec(w), 2))*
   Re(B0i(dbb0, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MB2, MB2)))/
  (Pi2*v2) + (1.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MD2, MD2)))/
  (Pi2*v2) + (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ME2, ME2)))/
  (Pi2*v2) + (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ML2, ML2)))/
  (Pi2*v2) + (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MM2, MM2)))/
  (Pi2*v2) + (0.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MT2, MT2)))/
  (Pi2*v2) + (1.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (1.5*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MW2, MW2)))/(Pi2*v2) - 
 (0.75*MW2*(pow(Sec(w), 2))*Re(B0i(dbb00, MZ2, 0, 0)))/(Pi2*v2) - 
 (0.08333333333333333*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Sec(w), 2))*
   Re(B0i(dbb00, MZ2, MB2, MB2)))/(Pi2*v2) - 
 (0.08333333333333333*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Sec(w), 2))*
   Re(B0i(dbb00, MZ2, MC2, MC2)))/(Pi2*v2) - 
 (0.08333333333333333*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Sec(w), 2))*
   Re(B0i(dbb00, MZ2, MD2, MD2)))/(Pi2*v2) - 
 (0.25*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Sec(w), 2))*
   Re(B0i(dbb00, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.25*MW2*(pow(Sec(w), 2))*Re(B0i(dbb00, MZ2, Mh2, MZ2)))/(Pi2*v2) - 
 (0.25*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Sec(w), 2))*
   Re(B0i(dbb00, MZ2, ML2, ML2)))/(Pi2*v2) - 
 (0.25*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Sec(w), 2))*
   Re(B0i(dbb00, MZ2, MM2, MM2)))/(Pi2*v2) - 
 (0.08333333333333333*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Sec(w), 2))*
   Re(B0i(dbb00, MZ2, MS2, MS2)))/(Pi2*v2) - 
 (0.08333333333333333*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Sec(w), 2))*
   Re(B0i(dbb00, MZ2, MT2, MT2)))/(Pi2*v2) - 
 (0.08333333333333333*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Sec(w), 2))*
   Re(B0i(dbb00, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.125*MW2*(7 + 8*Cos(2*w) + 3*Cos(4*w))*(pow(Sec(w), 2))*
   Re(B0i(dbb00, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MB, 2))*Re(B0i(dbb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MC, 2))*Re(B0i(dbb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MD, 2))*Re(B0i(dbb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ME, 2))*Re(B0i(dbb1, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ML, 2))*Re(B0i(dbb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(MM, 2))*Re(B0i(dbb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MS, 2))*Re(B0i(dbb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*Mh2*MT2*Re(B0i(dbb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MU, 2))*Re(B0i(dbb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*Mh2*MW2*Re(B0i(dbb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*Mh2*MW2*(pow(Sec(w), 2))*Re(B0i(dbb1, Mh2, MZ2, MZ2)))/
  (Pi2*v2) + (0.375*(pow(MW, 4))*(pow(Sec(w), 4))*
   Re(B0i(dbb1, MZ2, 0, 0)))/(Pi2*v2) + 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(MW, 4))*
   (pow(Sec(w), 4))*Re(B0i(dbb1, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(MW, 4))*
   (pow(Sec(w), 4))*Re(B0i(dbb1, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(MW, 4))*
   (pow(Sec(w), 4))*Re(B0i(dbb1, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(MW, 4))*(pow(Sec(w), 4))*
   Re(B0i(dbb1, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(MW, 4))*(pow(Sec(w), 4))*
   Re(B0i(dbb1, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(MW, 4))*(pow(Sec(w), 4))*
   Re(B0i(dbb1, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(MW, 4))*
   (pow(Sec(w), 4))*Re(B0i(dbb1, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(MW, 4))*
   (pow(Sec(w), 4))*Re(B0i(dbb1, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(MW, 4))*
   (pow(Sec(w), 4))*Re(B0i(dbb1, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.5*(pow(MW, 4))*Re(B0i(dbb1, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.25*MW2*C0i(cc0, m22, m32, m12, MW2, MW2, MW2)*(pow(Sin(w), 3))*
   (Cos(w)*Cot(w)*(2*m12 + m22 - m32 + 2*MW2 + 
      (-4*m12 + 6*m22 + 4*(m32 + MW2))*(pow(Cot(w), 2))) - 
    (Mh2 + 2*MW2)*Sin(w)))/(Pi2*v2)
;
}

ComplexType dKappahbb(double m12, double m22, double m32)
{
return 1. + (0.125*(MT2 + MW2)*B0i(bb0, m32, MT2, MW2))/(Pi2*v2) + 
 (0.125*m12*MT2*C0i(cc1, m12, m32, m22, MT2, MT2, MW2))/(Pi2*v2) - 
 (0.03125*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, MT2, MW2, MW2))/
  (Pi2*v2) + (0.0625*(m12 + m22 - m32)*MT2*C0i(cc2, m12, m32, m22, MT2, MT2, 
    MW2))/(Pi2*v2) + (0.03125*(3*m12 - 3*m22 - m32)*MW2*
   C0i(cc2, m22, m12, m32, MT2, MW2, MW2))/(Pi2*v2) - 
 (0.0625*B0i(bb0, m32, MB2, Mh2)*(pow(MB, 2)))/(Pi2*v2) - 
 (0.1875*Mh2*C0i(cc0, m22, m12, m32, MB2, Mh2, Mh2)*(pow(MB, 2)))/
  (Pi2*v2) - (0.0625*m12*C0i(cc1, m12, m32, m22, MB2, MB2, Mh2)*
   (pow(MB, 2)))/(Pi2*v2) - 
 (0.03125*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MB2, MB2, Mh2)*
   (pow(MB, 2)))/(Pi2*v2) - (0.125*C0i(cc0, m12, m32, m22, MB2, MB2, Mh2)*
   (pow(MB, 4)))/(Pi2*v2) + (0.25*C0i(cc0, m12, m32, m22, MT2, MT2, MW2)*
   (pow(MT, 4)))/(Pi2*v2) + (0.125*C0i(cc0, m22, m12, m32, MT2, MW2, MW2)*
   (Mh2*MT2 - m22*MW2 + (pow(MW, 4))))/(Pi2*v2) - 
 (0.015625*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, MB2, MZ2, MZ2)*
   (pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.015625*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, MB2, MZ2, MZ2)*
   (pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.001736111111111111*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MB2, MB2, 
    MZ2)*(-12*MW2 + 4*MW2*Cos(4*w) + 9*(pow(MB, 2)) + 
    Cos(2*w)*(8*MW2 + 9*(pow(MB, 2))))*(pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.006944444444444444*B0i(bb0, m32, MB2, MZ2)*(9*(pow(MB, 2)) + 
    MW2*(pow(1 + 2*Cos(2*w), 2))*(pow(Sec(w), 2))))/(Pi2*v2) + 
 (0.006944444444444444*C0i(cc0, m22, m12, m32, MB2, MZ2, MZ2)*
   (9*(pow(Csc(w), 2))*(Mh2*(pow(MB, 2)) + (pow(MW, 4))*
       (pow(Sec(w), 4))) - 4*(9*m22*MW2*(pow(Csc(2*w), 2)) + 
      4*(2 + Cos(2*w))*(pow(MW, 4))*(pow(Sec(w), 4)))))/
  (Pi2*v2*(pow(Csc(w), 2))) + (0.027777777777777776*B0i(bb0, m32, 0, MB2)*
   (48*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2))))/Pi2 - 
 (0.013888888888888888*(m12 + m22 - m32)*C0i(cc1, m22, m12, m32, 0, MB2, MB2)*
   (48*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2))))/Pi2 + 
 (0.013888888888888888*(m12 - m22 + m32)*C0i(cc2, m22, m12, m32, 0, MB2, MB2)*
   (48*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2))))/Pi2 + 
 (0.013888888888888888*C0i(cc0, m22, m12, m32, 0, MB2, MB2)*
   (-m12 - m22 + m32 + 4*(pow(MB, 2)))*
   (48*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2))))/Pi2 + 
 (0.006944444444444444*m12*C0i(cc1, m12, m32, m22, MB2, MB2, MZ2)*
   (9*(pow(MB, 2)) - 8*MW2*(2 + Cos(2*w))*(pow(Tan(w), 2))))/
  (Pi2*v2) + (0.013888888888888888*C0i(cc0, m12, m32, m22, MB2, MB2, MZ2)*
   (pow(MB, 2))*(9*(pow(MB, 2)) - 8*MW2*(2 + Cos(2*w))*
     (pow(Tan(w), 2))))/(Pi2*v2) - 
 (0.005208333333333333*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MB2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MC2)))/(Pi2*v2) - 
 (0.005208333333333333*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MD2)))/(Pi2*v2) - 
 (0.015625*(-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ME2)))/(Pi2*v2) - (0.03125*Re(A0i(aa0, Mh2)))/(Pi2*v2) - 
 (0.015625*(-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ML2)))/(Pi2*v2) - 
 (0.015625*(-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MM2)))/(Pi2*v2) - 
 (0.005208333333333333*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MS2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MT2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MU2)))/(Pi2*v2) + 
 (0.0078125*(-22 + Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MW2)))/(Pi2*v2) + (0.015625*(5 + 13*Cos(2*w) + 6*Cos(4*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MZ2)))/(Pi2*v2) - 
 (1.125*MW2*(pow(Sin(w), 2))*Re(B0i(bb0, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.027777777777777776*(-48*Alfas*Pi - 4*MW2*(pow(v, -2))*
     (pow(Sin(w), 2)))*Re(B0i(bb0, MB2, 0, MB2)))/Pi2 + 
 (0.0625*(pow(MB, 2))*Re(B0i(bb0, MB2, MB2, Mh2)))/(Pi2*v2) + 
 (0.006944444444444444*(-9*(pow(MB, 2)) + 8*MW2*(2 + Cos(2*w))*
     (pow(Tan(w), 2)))*Re(B0i(bb0, MB2, MB2, MZ2)))/(Pi2*v2) - 
 (0.125*MT2*Re(B0i(bb0, MB2, MT2, MW2)))/(Pi2*v2) + 
 (0.0625*MW2*Re(B0i(bb0, Mh2, MW2, MW2)))/(Pi2*v2) + 
 (0.03125*MW2*(pow(Sec(w), 2))*Re(B0i(bb0, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (0.5*MW2*Cos(2*w)*Re(B0i(bb0, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.375*(-MT2 + MW2)*Cos(2*w)*(pow(Csc(w), 2))*
   Re(B0i(bb0, MW2, MB2, MT2)))/(Pi2*v2) - 
 (0.375*(pow(MC, 2))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*Cos(2*w)*(MW2 - (pow(MU, 2)))*(pow(Csc(w), 2))*
   Re(B0i(bb0, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb0, MW2, Mh2, MW2)))/
  (Pi2*v2) - (0.015625*MW2*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Sec(w), 2))*Re(B0i(bb0, MW2, MW2, MZ2)))/
  (Pi2*v2) + (0.1875*(pow(MB, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.1875*(pow(MC, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MC2, MC2)))/
  (Pi2*v2) + (0.1875*(pow(MD, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*(pow(ME, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, ME2, ME2)))/
  (Pi2*v2) - (0.125*MW2*(pow(Csc(w), 2))*Re(B0i(bb0, MZ2, Mh2, MZ2)))/
  (Pi2*v2) + (0.0625*(pow(ML, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*(pow(MM, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MM2, MM2)))/
  (Pi2*v2) + (0.1875*(pow(MS, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.1875*MT2*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.1875*(pow(MU, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MU2, MU2)))/
  (Pi2*v2) + (0.0625*MW2*(5 + 9*Cos(2*w))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MB2, MB2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MD2, MD2)))/(Pi2*v2) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, ME2, ME2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, ML2, ML2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, MM2, MM2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MT2, MT2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.5*(2 + 3*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, MW2, MW2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, ME2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, ML2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, MM2)))/
  (Pi2*v2) - (1.*Cos(2*w)*Re(B0i(bb00, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MB2, MT2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.125*(1 - (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, Mh2, MW2)))/(Pi2*v2) - 
 (0.125*(2 + 5*Cos(2*w) + 2*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MW2, MW2, MZ2)))/(Pi2*v2) - 
 (0.375*(pow(Cot(w), 2))*Re(B0i(bb00, MZ2, 0, 0)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MB2, MB2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MC2, MC2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MD2, MD2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*(pow(Cot(w), 2))*Re(B0i(bb00, MZ2, Mh2, MZ2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ML2, ML2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MM2, MM2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MS2, MS2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MT2, MT2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.0625*(7 + 8*Cos(2*w) + 3*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MW2, MW2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MB2, MB2)))/
  (Pi2*v2) - (0.6666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MC2, MC2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MD2, MD2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ME2, ME2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ML2, ML2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MM2, MM2)))/
  (Pi2*v2) - (0.16666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MS2, MS2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MT2, MT2)))/
  (Pi2*v2) - (0.6666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.25*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.375*(pow(MB, 2))*Re(B0i(bb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*(pow(MC, 2))*Re(B0i(bb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*(pow(MD, 2))*Re(B0i(bb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*(pow(ME, 2))*Re(B0i(bb1, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*(pow(ML, 2))*Re(B0i(bb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*(pow(MM, 2))*Re(B0i(bb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*(pow(MS, 2))*Re(B0i(bb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*MT2*Re(B0i(bb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*(pow(MU, 2))*Re(B0i(bb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*MW2*Re(B0i(bb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*MW2*(pow(Sec(w), 2))*Re(B0i(bb1, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ME2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ML2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, MM2)))/(Pi2*v2) - 
 (0.25*MW2*Cos(2*w)*Re(B0i(bb1, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MB2, MT2)))/
  (Pi2*v2) - (0.375*MW2*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MD2, MU2)))/
  (Pi2*v2) + (0.25*MW2*Cos(2*w)*(pow(Cot(w), 2))*
   Re(B0i(bb1, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.1875*MW2*(pow(Csc(w), 2))*Re(B0i(bb1, MZ2, 0, 0)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.25*MW2*(pow(Cos(w), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb1, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.027777777777777776*(pow(MB, 2))*(48*Alfas*Pi + 
    4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*Re(B0i(dbb0, MB2, 0, MB2)))/
  Pi2 - (0.125*(pow(MB, 4))*Re(B0i(dbb0, MB2, MB2, Mh2)))/(Pi2*v2) + 
 (0.013888888888888888*(pow(MB, 2))*(9*(pow(MB, 2)) - 
    8*MW2*(2 + Cos(2*w))*(pow(Tan(w), 2)))*Re(B0i(dbb0, MB2, MB2, MZ2)))/
  (Pi2*v2) + (0.25*MT2*(pow(MB, 2))*Re(B0i(dbb0, MB2, MT2, MW2)))/
  (Pi2*v2) + (0.75*(pow(MB, 4))*Re(B0i(dbb0, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.75*(pow(MC, 4))*Re(B0i(dbb0, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.75*(pow(MD, 4))*Re(B0i(dbb0, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.25*(pow(ME, 4))*Re(B0i(dbb0, Mh2, ME2, ME2)))/(Pi2*v2) - 
 (0.140625*(pow(Mh, 4))*Re(B0i(dbb0, Mh2, Mh2, Mh2)))/(Pi2*v2) + 
 (0.25*(pow(ML, 4))*Re(B0i(dbb0, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.25*(pow(MM, 4))*Re(B0i(dbb0, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.75*(pow(MS, 4))*Re(B0i(dbb0, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.75*(pow(MT, 4))*Re(B0i(dbb0, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.75*(pow(MU, 4))*Re(B0i(dbb0, Mh2, MU2, MU2)))/(Pi2*v2) + 
 (0.03125*(2*Mh2*MW2 - (pow(Mh, 4)) - 12*(pow(MW, 4)))*
   Re(B0i(dbb0, Mh2, MW2, MW2)))/(Pi2*v2) + 
 (0.015625*(-(pow(Mh, 4)) + 2*Mh2*MW2*(pow(Sec(w), 2)) - 
    12*(pow(MW, 4))*(pow(Sec(w), 4)))*Re(B0i(dbb0, Mh2, MZ2, MZ2)))/
  (Pi2*v2) + (0.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MB2, MB2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MC2, MC2)))/
  (Pi2*v2) + (0.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MD2, MD2)))/(Pi2*v2) + 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ME2, ME2)))/(Pi2*v2) + 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ML2, ML2)))/(Pi2*v2) + 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MM2, MM2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MS2, MS2)))/
  (Pi2*v2) + (1.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MT2, MT2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MU2, MU2)))/
  (Pi2*v2) - (1.5*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MW2, MW2)))/
  (Pi2*v2) - (0.027777777777777776*(pow(MB, 2))*
   (48*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*
   Re(B0i(dbb1, MB2, 0, MB2)))/Pi2 + 
 (0.125*(pow(MB, 4))*Re(B0i(dbb1, MB2, MB2, Mh2)))/(Pi2*v2) + 
 (0.006944444444444444*(pow(MB, 2))*(12*MW2 + 2*MW2*Cos(4*w) + 
    9*(pow(MB, 2)) + Cos(2*w)*(4*MW2 + 9*(pow(MB, 2))))*
   (pow(Sec(w), 2))*Re(B0i(dbb1, MB2, MB2, MZ2)))/(Pi2*v2) + 
 (0.125*(pow(MB, 2))*(MT2 + 2*MW2 + (pow(MB, 2)))*
   Re(B0i(dbb1, MB2, MT2, MW2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MB, 2))*Re(B0i(dbb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MC, 2))*Re(B0i(dbb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MD, 2))*Re(B0i(dbb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ME, 2))*Re(B0i(dbb1, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ML, 2))*Re(B0i(dbb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(MM, 2))*Re(B0i(dbb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MS, 2))*Re(B0i(dbb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*Mh2*MT2*Re(B0i(dbb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MU, 2))*Re(B0i(dbb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*Mh2*MW2*Re(B0i(dbb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*Mh2*MW2*(pow(Sec(w), 2))*Re(B0i(dbb1, Mh2, MZ2, MZ2)))/(Pi2*v2)
;
}

ComplexType dKappahcc(double m12, double m22, double m32)
{
return 1. - (0.03125*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, MS2, MW2, MW2))/
  (Pi2*v2) + (0.03125*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, MS2, 
    MW2, MW2))/(Pi2*v2) - (0.0625*B0i(bb0, m32, MC2, Mh2)*(pow(MC, 2)))/
  (Pi2*v2) - (0.1875*Mh2*C0i(cc0, m22, m12, m32, MC2, Mh2, Mh2)*
   (pow(MC, 2)))/(Pi2*v2) - 
 (0.0625*m12*C0i(cc1, m12, m32, m22, MC2, MC2, Mh2)*(pow(MC, 2)))/
  (Pi2*v2) - (0.03125*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MC2, MC2, 
    Mh2)*(pow(MC, 2)))/(Pi2*v2) - 
 (0.125*C0i(cc0, m12, m32, m22, MC2, MC2, Mh2)*(pow(MC, 4)))/(Pi2*v2) + 
 (0.125*m12*C0i(cc1, m12, m32, m22, MS2, MS2, MW2)*(pow(MS, 2)))/
  (Pi2*v2) + (0.0625*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MS2, MS2, MW2)*
   (pow(MS, 2)))/(Pi2*v2) + 
 (0.125*B0i(bb0, m32, MS2, MW2)*(MW2 + (pow(MS, 2))))/(Pi2*v2) + 
 (0.25*C0i(cc0, m12, m32, m22, MS2, MS2, MW2)*(pow(MS, 4)))/(Pi2*v2) + 
 (0.125*C0i(cc0, m22, m12, m32, MS2, MW2, MW2)*
   (-(m22*MW2) + Mh2*(pow(MS, 2)) + (pow(MW, 4))))/(Pi2*v2) - 
 (0.015625*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, MC2, MZ2, MZ2)*
   (pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.015625*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, MC2, MZ2, MZ2)*
   (pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.001736111111111111*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MC2, MC2, 
    MZ2)*(16*MW2*Cos(4*w) + 9*(pow(MC, 2)) + 
    Cos(2*w)*(-16*MW2 + 9*(pow(MC, 2))))*(pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.006944444444444444*B0i(bb0, m32, MC2, MZ2)*(9*(pow(MC, 2)) + 
    MW2*(pow(1 - 4*Cos(2*w), 2))*(pow(Sec(w), 2))))/(Pi2*v2) + 
 (0.006944444444444444*C0i(cc0, m22, m12, m32, MC2, MZ2, MZ2)*
   (9*(pow(Csc(w), 2))*(Mh2*(pow(MC, 2)) + (pow(MW, 4))*
       (pow(Sec(w), 4))) - 4*(9*m22*MW2*(pow(Csc(2*w), 2)) + 
      8*(1 + 2*Cos(2*w))*(pow(MW, 4))*(pow(Sec(w), 4)))))/
  (Pi2*v2*(pow(Csc(w), 2))) - 
 (0.2222222222222222*(m12 + m22 - m32)*C0i(cc1, m22, m12, m32, 0, MC2, MC2)*
   (3*Alfas*Pi*v2 + MW2*(pow(Sin(w), 2))))/(Pi2*v2) + 
 (0.2222222222222222*(m12 - m22 + m32)*C0i(cc2, m22, m12, m32, 0, MC2, MC2)*
   (3*Alfas*Pi*v2 + MW2*(pow(Sin(w), 2))))/(Pi2*v2) - 
 (0.2222222222222222*C0i(cc0, m22, m12, m32, 0, MC2, MC2)*
   (m12 + m22 - m32 - 4*(pow(MC, 2)))*(3*Alfas*Pi*v2 + 
    MW2*(pow(Sin(w), 2))))/(Pi2*v2) + 
 (0.1111111111111111*B0i(bb0, m32, 0, MC2)*
   (12*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2))))/Pi2 + 
 (0.006944444444444444*m12*C0i(cc1, m12, m32, m22, MC2, MC2, MZ2)*
   (9*(pow(MC, 2)) - 16*MW2*(1 + 2*Cos(2*w))*(pow(Tan(w), 2))))/
  (Pi2*v2) + (0.013888888888888888*C0i(cc0, m12, m32, m22, MC2, MC2, MZ2)*
   (pow(MC, 2))*(9*(pow(MC, 2)) - 16*MW2*(1 + 2*Cos(2*w))*
     (pow(Tan(w), 2))))/(Pi2*v2) - 
 (0.005208333333333333*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MB2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MC2)))/(Pi2*v2) - 
 (0.005208333333333333*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MD2)))/(Pi2*v2) - 
 (0.015625*(-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ME2)))/(Pi2*v2) - (0.03125*Re(A0i(aa0, Mh2)))/(Pi2*v2) - 
 (0.015625*(-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ML2)))/(Pi2*v2) - 
 (0.015625*(-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MM2)))/(Pi2*v2) - 
 (0.005208333333333333*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MS2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MT2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MU2)))/(Pi2*v2) + 
 (0.0078125*(-22 + Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MW2)))/(Pi2*v2) + (0.015625*(5 + 13*Cos(2*w) + 6*Cos(4*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MZ2)))/(Pi2*v2) - 
 (1.125*MW2*(pow(Sin(w), 2))*Re(B0i(bb0, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.1111111111111111*(-12*Alfas*Pi - 4*MW2*(pow(v, -2))*
     (pow(Sin(w), 2)))*Re(B0i(bb0, MC2, 0, MC2)))/Pi2 + 
 (0.0625*(pow(MC, 2))*Re(B0i(bb0, MC2, MC2, Mh2)))/(Pi2*v2) + 
 (0.006944444444444444*(-9*(pow(MC, 2)) + 16*MW2*(1 + 2*Cos(2*w))*
     (pow(Tan(w), 2)))*Re(B0i(bb0, MC2, MC2, MZ2)))/(Pi2*v2) - 
 (0.125*(pow(MS, 2))*Re(B0i(bb0, MC2, MS2, MW2)))/(Pi2*v2) + 
 (0.0625*MW2*Re(B0i(bb0, Mh2, MW2, MW2)))/(Pi2*v2) + 
 (0.03125*MW2*(pow(Sec(w), 2))*Re(B0i(bb0, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (0.5*MW2*Cos(2*w)*Re(B0i(bb0, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.375*(-MT2 + MW2)*Cos(2*w)*(pow(Csc(w), 2))*
   Re(B0i(bb0, MW2, MB2, MT2)))/(Pi2*v2) - 
 (0.375*(pow(MC, 2))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*Cos(2*w)*(MW2 - (pow(MU, 2)))*(pow(Csc(w), 2))*
   Re(B0i(bb0, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb0, MW2, Mh2, MW2)))/
  (Pi2*v2) - (0.015625*MW2*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Sec(w), 2))*Re(B0i(bb0, MW2, MW2, MZ2)))/
  (Pi2*v2) + (0.1875*(pow(MB, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.1875*(pow(MC, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MC2, MC2)))/
  (Pi2*v2) + (0.1875*(pow(MD, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*(pow(ME, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, ME2, ME2)))/
  (Pi2*v2) - (0.125*MW2*(pow(Csc(w), 2))*Re(B0i(bb0, MZ2, Mh2, MZ2)))/
  (Pi2*v2) + (0.0625*(pow(ML, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*(pow(MM, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MM2, MM2)))/
  (Pi2*v2) + (0.1875*(pow(MS, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.1875*MT2*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.1875*(pow(MU, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MU2, MU2)))/
  (Pi2*v2) + (0.0625*MW2*(5 + 9*Cos(2*w))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MB2, MB2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MD2, MD2)))/(Pi2*v2) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, ME2, ME2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, ML2, ML2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, MM2, MM2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MT2, MT2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.5*(2 + 3*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, MW2, MW2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, ME2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, ML2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, MM2)))/
  (Pi2*v2) - (1.*Cos(2*w)*Re(B0i(bb00, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MB2, MT2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.125*(1 - (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, Mh2, MW2)))/(Pi2*v2) - 
 (0.125*(2 + 5*Cos(2*w) + 2*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MW2, MW2, MZ2)))/(Pi2*v2) - 
 (0.375*(pow(Cot(w), 2))*Re(B0i(bb00, MZ2, 0, 0)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MB2, MB2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MC2, MC2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MD2, MD2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*(pow(Cot(w), 2))*Re(B0i(bb00, MZ2, Mh2, MZ2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ML2, ML2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MM2, MM2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MS2, MS2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MT2, MT2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.0625*(7 + 8*Cos(2*w) + 3*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MW2, MW2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MB2, MB2)))/
  (Pi2*v2) - (0.6666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MC2, MC2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MD2, MD2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ME2, ME2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ML2, ML2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MM2, MM2)))/
  (Pi2*v2) - (0.16666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MS2, MS2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MT2, MT2)))/
  (Pi2*v2) - (0.6666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.25*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.375*(pow(MB, 2))*Re(B0i(bb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*(pow(MC, 2))*Re(B0i(bb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*(pow(MD, 2))*Re(B0i(bb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*(pow(ME, 2))*Re(B0i(bb1, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*(pow(ML, 2))*Re(B0i(bb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*(pow(MM, 2))*Re(B0i(bb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*(pow(MS, 2))*Re(B0i(bb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*MT2*Re(B0i(bb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*(pow(MU, 2))*Re(B0i(bb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*MW2*Re(B0i(bb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*MW2*(pow(Sec(w), 2))*Re(B0i(bb1, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ME2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ML2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, MM2)))/(Pi2*v2) - 
 (0.25*MW2*Cos(2*w)*Re(B0i(bb1, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MB2, MT2)))/
  (Pi2*v2) - (0.375*MW2*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MD2, MU2)))/
  (Pi2*v2) + (0.25*MW2*Cos(2*w)*(pow(Cot(w), 2))*
   Re(B0i(bb1, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.1875*MW2*(pow(Csc(w), 2))*Re(B0i(bb1, MZ2, 0, 0)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.25*MW2*(pow(Cos(w), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb1, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.1111111111111111*(pow(MC, 2))*(12*Alfas*Pi + 
    4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*Re(B0i(dbb0, MC2, 0, MC2)))/
  Pi2 - (0.125*(pow(MC, 4))*Re(B0i(dbb0, MC2, MC2, Mh2)))/(Pi2*v2) + 
 (0.013888888888888888*(pow(MC, 2))*(9*(pow(MC, 2)) - 
    16*MW2*(1 + 2*Cos(2*w))*(pow(Tan(w), 2)))*
   Re(B0i(dbb0, MC2, MC2, MZ2)))/(Pi2*v2) + 
 (0.25*(pow(MC, 2))*(pow(MS, 2))*Re(B0i(dbb0, MC2, MS2, MW2)))/
  (Pi2*v2) + (0.75*(pow(MB, 4))*Re(B0i(dbb0, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.75*(pow(MC, 4))*Re(B0i(dbb0, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.75*(pow(MD, 4))*Re(B0i(dbb0, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.25*(pow(ME, 4))*Re(B0i(dbb0, Mh2, ME2, ME2)))/(Pi2*v2) - 
 (0.140625*(pow(Mh, 4))*Re(B0i(dbb0, Mh2, Mh2, Mh2)))/(Pi2*v2) + 
 (0.25*(pow(ML, 4))*Re(B0i(dbb0, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.25*(pow(MM, 4))*Re(B0i(dbb0, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.75*(pow(MS, 4))*Re(B0i(dbb0, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.75*(pow(MT, 4))*Re(B0i(dbb0, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.75*(pow(MU, 4))*Re(B0i(dbb0, Mh2, MU2, MU2)))/(Pi2*v2) + 
 (0.03125*(2*Mh2*MW2 - (pow(Mh, 4)) - 12*(pow(MW, 4)))*
   Re(B0i(dbb0, Mh2, MW2, MW2)))/(Pi2*v2) + 
 (0.015625*(-(pow(Mh, 4)) + 2*Mh2*MW2*(pow(Sec(w), 2)) - 
    12*(pow(MW, 4))*(pow(Sec(w), 4)))*Re(B0i(dbb0, Mh2, MZ2, MZ2)))/
  (Pi2*v2) + (0.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MB2, MB2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MC2, MC2)))/
  (Pi2*v2) + (0.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MD2, MD2)))/(Pi2*v2) + 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ME2, ME2)))/(Pi2*v2) + 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ML2, ML2)))/(Pi2*v2) + 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MM2, MM2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MS2, MS2)))/
  (Pi2*v2) + (1.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MT2, MT2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MU2, MU2)))/
  (Pi2*v2) - (1.5*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MW2, MW2)))/
  (Pi2*v2) - (0.1111111111111111*(pow(MC, 2))*
   (12*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*
   Re(B0i(dbb1, MC2, 0, MC2)))/Pi2 + 
 (0.125*(pow(MC, 4))*Re(B0i(dbb1, MC2, MC2, Mh2)))/(Pi2*v2) + 
 (0.006944444444444444*(pow(MC, 2))*(8*MW2*Cos(4*w) + 
    9*(2*MW2 + (pow(MC, 2))) + Cos(2*w)*(-8*MW2 + 9*(pow(MC, 2))))*
   (pow(Sec(w), 2))*Re(B0i(dbb1, MC2, MC2, MZ2)))/(Pi2*v2) + 
 (0.125*(pow(MC, 2))*(2*MW2 + (pow(MC, 2)) + (pow(MS, 2)))*
   Re(B0i(dbb1, MC2, MS2, MW2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MB, 2))*Re(B0i(dbb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MC, 2))*Re(B0i(dbb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MD, 2))*Re(B0i(dbb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ME, 2))*Re(B0i(dbb1, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ML, 2))*Re(B0i(dbb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(MM, 2))*Re(B0i(dbb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MS, 2))*Re(B0i(dbb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*Mh2*MT2*Re(B0i(dbb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MU, 2))*Re(B0i(dbb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*Mh2*MW2*Re(B0i(dbb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*Mh2*MW2*(pow(Sec(w), 2))*Re(B0i(dbb1, Mh2, MZ2, MZ2)))/(Pi2*v2)
;
}

ComplexType dKappahtautau(double m12, double m22, double m32)
{
return 1. + (0.125*MW2*B0i(bb0, m32, 0, MW2))/(Pi2*v2) + 
 (0.125*MW2*(-m22 + MW2)*C0i(cc0, m22, m12, m32, 0, MW2, MW2))/(Pi2*v2) - 
 (0.03125*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, 0, MW2, MW2))/
  (Pi2*v2) + (0.03125*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, 0, 
    MW2, MW2))/(Pi2*v2) - (0.0625*B0i(bb0, m32, Mh2, ML2)*(pow(ML, 2)))/
  (Pi2*v2) - (0.1875*Mh2*C0i(cc0, m12, m32, m22, Mh2, Mh2, ML2)*
   (pow(ML, 2)))/(Pi2*v2) + 
 (0.03125*(m12 + m22 - m32)*C0i(cc1, m22, m12, m32, Mh2, ML2, ML2)*
   (pow(ML, 2)))/(Pi2*v2) - 
 (0.03125*(m12 - m22 + m32)*C0i(cc2, m22, m12, m32, Mh2, ML2, ML2)*
   (pow(ML, 2)))/(Pi2*v2) + 
 (0.03125*C0i(cc0, m22, m12, m32, Mh2, ML2, ML2)*
   (m12 + m22 - m32 - 4*(pow(ML, 2)))*(pow(ML, 2)))/(Pi2*v2) - 
 (0.015625*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, ML2, MZ2, MZ2)*
   (pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.015625*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, ML2, MZ2, MZ2)*
   (pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.015625*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, ML2, ML2, MZ2)*
   (4*MW2 + 4*MW2*Cos(4*w) + (pow(ML, 2)) + 
    Cos(2*w)*(-8*MW2 + (pow(ML, 2))))*(pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.0625*B0i(bb0, m32, ML2, MZ2)*((pow(ML, 2)) + 
    MW2*(pow(1 - 2*Cos(2*w), 2))*(pow(Sec(w), 2))))/(Pi2*v2) + 
 (1.*MW2*B0i(bb0, m32, 0, ML2)*(pow(Sin(w), 2)))/(Pi2*v2) - 
 (0.5*(m12 + m22 - m32)*MW2*C0i(cc1, m22, m12, m32, 0, ML2, ML2)*
   (pow(Sin(w), 2)))/(Pi2*v2) + 
 (0.5*(m12 - m22 + m32)*MW2*C0i(cc2, m22, m12, m32, 0, ML2, ML2)*
   (pow(Sin(w), 2)))/(Pi2*v2) - 
 (0.5*MW2*C0i(cc0, m22, m12, m32, 0, ML2, ML2)*
   (m12 + m22 - m32 - 4*(pow(ML, 2)))*(pow(Sin(w), 2)))/(Pi2*v2) + 
 (0.0625*m12*C0i(cc1, m12, m32, m22, ML2, ML2, MZ2)*
   ((pow(ML, 2)) - 8*MW2*Cos(2*w)*(pow(Tan(w), 2))))/(Pi2*v2) + 
 (0.125*C0i(cc0, m12, m32, m22, ML2, ML2, MZ2)*((pow(ML, 4)) - 
    8*MW2*Cos(2*w)*(pow(ML, 2))*(pow(Tan(w), 2))))/(Pi2*v2) + 
 (0.0625*C0i(cc0, m22, m12, m32, ML2, MZ2, MZ2)*(Mh2*(pow(ML, 2)) - 
    m22*MW2*(pow(Sec(w), 2)) + (pow(MW, 4))*(pow(Sec(w), 4)) + 
    16*(pow(MW, 4))*(-1 + (pow(Tan(w), 2)))*(pow(Tan(w), 2))))/
  (Pi2*v2) - (0.005208333333333333*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + 
    Cos(6*w))*(pow(Csc(w), 2))*Re(A0i(aa0, MB2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MC2)))/(Pi2*v2) - 
 (0.005208333333333333*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MD2)))/(Pi2*v2) - 
 (0.015625*(-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ME2)))/(Pi2*v2) - (0.03125*Re(A0i(aa0, Mh2)))/(Pi2*v2) - 
 (0.015625*(-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ML2)))/(Pi2*v2) - 
 (0.015625*(-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MM2)))/(Pi2*v2) - 
 (0.005208333333333333*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MS2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MT2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MU2)))/(Pi2*v2) + 
 (0.0078125*(-22 + Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MW2)))/(Pi2*v2) + (0.015625*(5 + 13*Cos(2*w) + 6*Cos(4*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MZ2)))/(Pi2*v2) - 
 (1.125*MW2*(pow(Sin(w), 2))*Re(B0i(bb0, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.0625*MW2*Re(B0i(bb0, Mh2, MW2, MW2)))/(Pi2*v2) + 
 (0.03125*MW2*(pow(Sec(w), 2))*Re(B0i(bb0, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (1.*MW2*(pow(Sin(w), 2))*Re(B0i(bb0, ML2, 0, ML2)))/(Pi2*v2) + 
 (0.0625*(pow(ML, 2))*Re(B0i(bb0, ML2, Mh2, ML2)))/(Pi2*v2) + 
 (0.0625*(-(pow(ML, 2)) + 8*MW2*Cos(2*w)*(pow(Tan(w), 2)))*
   Re(B0i(bb0, ML2, ML2, MZ2)))/(Pi2*v2) - 
 (0.5*MW2*Cos(2*w)*Re(B0i(bb0, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.375*(-MT2 + MW2)*Cos(2*w)*(pow(Csc(w), 2))*
   Re(B0i(bb0, MW2, MB2, MT2)))/(Pi2*v2) - 
 (0.375*(pow(MC, 2))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*Cos(2*w)*(MW2 - (pow(MU, 2)))*(pow(Csc(w), 2))*
   Re(B0i(bb0, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb0, MW2, Mh2, MW2)))/
  (Pi2*v2) - (0.015625*MW2*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Sec(w), 2))*Re(B0i(bb0, MW2, MW2, MZ2)))/
  (Pi2*v2) + (0.1875*(pow(MB, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.1875*(pow(MC, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MC2, MC2)))/
  (Pi2*v2) + (0.1875*(pow(MD, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*(pow(ME, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, ME2, ME2)))/
  (Pi2*v2) - (0.125*MW2*(pow(Csc(w), 2))*Re(B0i(bb0, MZ2, Mh2, MZ2)))/
  (Pi2*v2) + (0.0625*(pow(ML, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*(pow(MM, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MM2, MM2)))/
  (Pi2*v2) + (0.1875*(pow(MS, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.1875*MT2*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.1875*(pow(MU, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MU2, MU2)))/
  (Pi2*v2) + (0.0625*MW2*(5 + 9*Cos(2*w))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MB2, MB2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MD2, MD2)))/(Pi2*v2) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, ME2, ME2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, ML2, ML2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, MM2, MM2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MT2, MT2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.5*(2 + 3*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, MW2, MW2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, ME2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, ML2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, MM2)))/
  (Pi2*v2) - (1.*Cos(2*w)*Re(B0i(bb00, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MB2, MT2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.125*(1 - (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, Mh2, MW2)))/(Pi2*v2) - 
 (0.125*(2 + 5*Cos(2*w) + 2*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MW2, MW2, MZ2)))/(Pi2*v2) - 
 (0.375*(pow(Cot(w), 2))*Re(B0i(bb00, MZ2, 0, 0)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MB2, MB2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MC2, MC2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MD2, MD2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*(pow(Cot(w), 2))*Re(B0i(bb00, MZ2, Mh2, MZ2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ML2, ML2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MM2, MM2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MS2, MS2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MT2, MT2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.0625*(7 + 8*Cos(2*w) + 3*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MW2, MW2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MB2, MB2)))/
  (Pi2*v2) - (0.6666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MC2, MC2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MD2, MD2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ME2, ME2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ML2, ML2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MM2, MM2)))/
  (Pi2*v2) - (0.16666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MS2, MS2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MT2, MT2)))/
  (Pi2*v2) - (0.6666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.25*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.375*(pow(MB, 2))*Re(B0i(bb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*(pow(MC, 2))*Re(B0i(bb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*(pow(MD, 2))*Re(B0i(bb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*(pow(ME, 2))*Re(B0i(bb1, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*(pow(ML, 2))*Re(B0i(bb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*(pow(MM, 2))*Re(B0i(bb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*(pow(MS, 2))*Re(B0i(bb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*MT2*Re(B0i(bb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*(pow(MU, 2))*Re(B0i(bb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*MW2*Re(B0i(bb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*MW2*(pow(Sec(w), 2))*Re(B0i(bb1, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ME2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ML2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, MM2)))/(Pi2*v2) - 
 (0.25*MW2*Cos(2*w)*Re(B0i(bb1, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MB2, MT2)))/
  (Pi2*v2) - (0.375*MW2*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MD2, MU2)))/
  (Pi2*v2) + (0.25*MW2*Cos(2*w)*(pow(Cot(w), 2))*
   Re(B0i(bb1, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.1875*MW2*(pow(Csc(w), 2))*Re(B0i(bb1, MZ2, 0, 0)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.25*MW2*(pow(Cos(w), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb1, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.75*(pow(MB, 4))*Re(B0i(dbb0, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.75*(pow(MC, 4))*Re(B0i(dbb0, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.75*(pow(MD, 4))*Re(B0i(dbb0, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.25*(pow(ME, 4))*Re(B0i(dbb0, Mh2, ME2, ME2)))/(Pi2*v2) - 
 (0.140625*(pow(Mh, 4))*Re(B0i(dbb0, Mh2, Mh2, Mh2)))/(Pi2*v2) + 
 (0.25*(pow(ML, 4))*Re(B0i(dbb0, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.25*(pow(MM, 4))*Re(B0i(dbb0, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.75*(pow(MS, 4))*Re(B0i(dbb0, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.75*(pow(MT, 4))*Re(B0i(dbb0, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.75*(pow(MU, 4))*Re(B0i(dbb0, Mh2, MU2, MU2)))/(Pi2*v2) + 
 (0.03125*(2*Mh2*MW2 - (pow(Mh, 4)) - 12*(pow(MW, 4)))*
   Re(B0i(dbb0, Mh2, MW2, MW2)))/(Pi2*v2) + 
 (0.015625*(-(pow(Mh, 4)) + 2*Mh2*MW2*(pow(Sec(w), 2)) - 
    12*(pow(MW, 4))*(pow(Sec(w), 4)))*Re(B0i(dbb0, Mh2, MZ2, MZ2)))/
  (Pi2*v2) + (1.*MW2*(pow(ML, 2))*(pow(Sin(w), 2))*
   Re(B0i(dbb0, ML2, 0, ML2)))/(Pi2*v2) - 
 (0.25*(pow(ML, 4))*Re(B0i(dbb0, ML2, Mh2, ML2)))/(Pi2*v2) + 
 (0.125*((pow(ML, 4)) - 8*MW2*Cos(2*w)*(pow(ML, 2))*
     (pow(Tan(w), 2)))*Re(B0i(dbb0, ML2, ML2, MZ2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MB2, MB2)))/
  (Pi2*v2) + (1.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MD2, MD2)))/
  (Pi2*v2) + (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ME2, ME2)))/
  (Pi2*v2) + (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ML2, ML2)))/
  (Pi2*v2) + (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MM2, MM2)))/
  (Pi2*v2) + (0.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MT2, MT2)))/
  (Pi2*v2) + (1.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (1.5*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MB, 2))*Re(B0i(dbb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MC, 2))*Re(B0i(dbb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MD, 2))*Re(B0i(dbb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ME, 2))*Re(B0i(dbb1, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ML, 2))*Re(B0i(dbb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(MM, 2))*Re(B0i(dbb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MS, 2))*Re(B0i(dbb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*Mh2*MT2*Re(B0i(dbb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MU, 2))*Re(B0i(dbb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*Mh2*MW2*Re(B0i(dbb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*Mh2*MW2*(pow(Sec(w), 2))*Re(B0i(dbb1, Mh2, MZ2, MZ2)))/
  (Pi2*v2) - (1.*MW2*(pow(ML, 2))*(pow(Sin(w), 2))*
   Re(B0i(dbb1, ML2, 0, ML2)))/(Pi2*v2) + 
 (0.125*(2*MW2*(pow(ML, 2)) + (pow(ML, 4)))*
   Re(B0i(dbb1, ML2, 0, MW2)))/(Pi2*v2) - 
 (0.125*(pow(ML, 4))*Re(B0i(dbb1, ML2, Mh2, ML2)))/(Pi2*v2) + 
 (0.0625*(pow(ML, 2))*(4*MW2 + 2*MW2*Cos(4*w) + (pow(ML, 2)) + 
    Cos(2*w)*(-4*MW2 + (pow(ML, 2))))*(pow(Sec(w), 2))*
   Re(B0i(dbb1, ML2, ML2, MZ2)))/(Pi2*v2)
;
}

ComplexType dKappahtt(double m12, double m22, double m32)
{
return 1. - (0.0625*MT2*B0i(bb0, m32, Mh2, MT2))/(Pi2*v2) - 
 (0.1875*Mh2*MT2*C0i(cc0, m12, m32, m22, Mh2, Mh2, MT2))/(Pi2*v2) + 
 (0.03125*(m12 + m22 - m32 - 4*MT2)*MT2*C0i(cc0, m22, m12, m32, Mh2, MT2, 
    MT2))/(Pi2*v2) - (0.03125*(m12 + 5*m22 - m32)*MW2*
   C0i(cc1, m22, m12, m32, MB2, MW2, MW2))/(Pi2*v2) + 
 (0.03125*(m12 + m22 - m32)*MT2*C0i(cc1, m22, m12, m32, Mh2, MT2, MT2))/
  (Pi2*v2) + (0.03125*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, MB2, 
    MW2, MW2))/(Pi2*v2) - (0.03125*(m12 - m22 + m32)*MT2*
   C0i(cc2, m22, m12, m32, Mh2, MT2, MT2))/(Pi2*v2) + 
 (0.125*m12*C0i(cc1, m12, m32, m22, MB2, MB2, MW2)*(pow(MB, 2)))/
  (Pi2*v2) + (0.0625*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MB2, MB2, MW2)*
   (pow(MB, 2)))/(Pi2*v2) + 
 (0.125*B0i(bb0, m32, MB2, MW2)*(MW2 + (pow(MB, 2))))/(Pi2*v2) + 
 (0.25*C0i(cc0, m12, m32, m22, MB2, MB2, MW2)*(pow(MB, 4)))/(Pi2*v2) + 
 (0.125*C0i(cc0, m22, m12, m32, MB2, MW2, MW2)*
   (-(m22*MW2) + Mh2*(pow(MB, 2)) + (pow(MW, 4))))/(Pi2*v2) - 
 (0.015625*(m12 + 5*m22 - m32)*MW2*C0i(cc1, m22, m12, m32, MT2, MZ2, MZ2)*
   (pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.015625*(3*m12 - 3*m22 - m32)*MW2*C0i(cc2, m22, m12, m32, MT2, MZ2, MZ2)*
   (pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.001736111111111111*(m12 + m22 - m32)*C0i(cc2, m12, m32, m22, MT2, MT2, 
    MZ2)*(9*MT2 + (9*MT2 - 16*MW2)*Cos(2*w) + 16*MW2*Cos(4*w))*
   (pow(Sec(w), 2)))/(Pi2*v2) + 
 (0.006944444444444444*B0i(bb0, m32, MT2, MZ2)*
   (9*MT2 + MW2*(pow(1 - 4*Cos(2*w), 2))*(pow(Sec(w), 2))))/
  (Pi2*v2) + (0.006944444444444444*C0i(cc0, m22, m12, m32, MT2, MZ2, MZ2)*
   (9*(pow(Csc(w), 2))*(Mh2*MT2 + (pow(MW, 4))*
       (pow(Sec(w), 4))) - 4*(9*m22*MW2*(pow(Csc(2*w), 2)) + 
      8*(1 + 2*Cos(2*w))*(pow(MW, 4))*(pow(Sec(w), 4)))))/
  (Pi2*v2*(pow(Csc(w), 2))) - 
 (0.2222222222222222*(m12 + m22 - m32 - 4*MT2)*C0i(cc0, m22, m12, m32, 0, 
    MT2, MT2)*(3*Alfas*Pi*v2 + MW2*(pow(Sin(w), 2))))/(Pi2*v2) - 
 (0.2222222222222222*(m12 + m22 - m32)*C0i(cc1, m22, m12, m32, 0, MT2, MT2)*
   (3*Alfas*Pi*v2 + MW2*(pow(Sin(w), 2))))/(Pi2*v2) + 
 (0.2222222222222222*(m12 - m22 + m32)*C0i(cc2, m22, m12, m32, 0, MT2, MT2)*
   (3*Alfas*Pi*v2 + MW2*(pow(Sin(w), 2))))/(Pi2*v2) + 
 (0.1111111111111111*B0i(bb0, m32, 0, MT2)*
   (12*Alfas*Pi + 4*MW2*(pow(v, -2))*(pow(Sin(w), 2))))/Pi2 + 
 (0.013888888888888888*MT2*C0i(cc0, m12, m32, m22, MT2, MT2, MZ2)*
   (9*MT2 - 16*MW2*(1 + 2*Cos(2*w))*(pow(Tan(w), 2))))/(Pi2*v2) + 
 (0.006944444444444444*m12*C0i(cc1, m12, m32, m22, MT2, MT2, MZ2)*
   (9*MT2 - 16*MW2*(1 + 2*Cos(2*w))*(pow(Tan(w), 2))))/(Pi2*v2) - 
 (0.005208333333333333*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MB2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MC2)))/(Pi2*v2) - 
 (0.005208333333333333*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MD2)))/(Pi2*v2) - 
 (0.015625*(-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ME2)))/(Pi2*v2) - (0.03125*Re(A0i(aa0, Mh2)))/(Pi2*v2) - 
 (0.015625*(-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, ML2)))/(Pi2*v2) - 
 (0.015625*(-16 + 29*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MM2)))/(Pi2*v2) - 
 (0.005208333333333333*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MS2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MT2)))/(Pi2*v2) - 
 (0.010416666666666666*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MU2)))/(Pi2*v2) + 
 (0.0078125*(-22 + Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*(pow(Csc(w), 2))*
   Re(A0i(aa0, MW2)))/(Pi2*v2) + (0.015625*(5 + 13*Cos(2*w) + 6*Cos(4*w))*
   (pow(Csc(w), 2))*Re(A0i(aa0, MZ2)))/(Pi2*v2) - 
 (1.125*MW2*(pow(Sin(w), 2))*Re(B0i(bb0, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.0625*MW2*Re(B0i(bb0, Mh2, MW2, MW2)))/(Pi2*v2) + 
 (0.03125*MW2*(pow(Sec(w), 2))*Re(B0i(bb0, Mh2, MZ2, MZ2)))/(Pi2*v2) + 
 (0.1111111111111111*(-12*Alfas*Pi - 4*MW2*(pow(v, -2))*
     (pow(Sin(w), 2)))*Re(B0i(bb0, MT2, 0, MT2)))/Pi2 - 
 (0.125*(pow(MB, 2))*Re(B0i(bb0, MT2, MB2, MW2)))/(Pi2*v2) + 
 (0.0625*MT2*Re(B0i(bb0, MT2, Mh2, MT2)))/(Pi2*v2) + 
 (0.006944444444444444*(-9*MT2 + 16*MW2*(1 + 2*Cos(2*w))*
     (pow(Tan(w), 2)))*Re(B0i(bb0, MT2, MT2, MZ2)))/(Pi2*v2) - 
 (0.5*MW2*Cos(2*w)*Re(B0i(bb0, MW2, 0, MW2)))/(Pi2*v2) - 
 (0.375*(MT2 - MW2)*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb0, MW2, MB2, MT2)))/
  (Pi2*v2) - (0.375*(pow(MC, 2))*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb0, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*Cos(2*w)*(MW2 - (pow(MU, 2)))*(pow(Csc(w), 2))*
   Re(B0i(bb0, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb0, MW2, Mh2, MW2)))/
  (Pi2*v2) - (0.015625*MW2*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*
   (pow(Csc(w), 2))*(pow(Sec(w), 2))*Re(B0i(bb0, MW2, MW2, MZ2)))/
  (Pi2*v2) + (0.1875*(pow(MB, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.1875*(pow(MC, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MC2, MC2)))/
  (Pi2*v2) + (0.1875*(pow(MD, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*(pow(ME, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, ME2, ME2)))/
  (Pi2*v2) - (0.125*MW2*(pow(Csc(w), 2))*Re(B0i(bb0, MZ2, Mh2, MZ2)))/
  (Pi2*v2) + (0.0625*(pow(ML, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*(pow(MM, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MM2, MM2)))/
  (Pi2*v2) + (0.1875*(pow(MS, 2))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.1875*MT2*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.1875*(pow(MU, 2))*(pow(Cot(w), 2))*Re(B0i(bb0, MZ2, MU2, MU2)))/
  (Pi2*v2) + (0.0625*MW2*(5 + 9*Cos(2*w))*(pow(Cot(w), 2))*
   Re(B0i(bb0, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MB2, MB2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MD2, MD2)))/(Pi2*v2) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, ME2, ME2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, ML2, ML2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.5*(-4 + (pow(Csc(w), 2)))*Re(B0i(bb00, 0, MM2, MM2)))/
  (Pi2*v2*(pow(Csc(w), 4))) + 
 (0.16666666666666666*(1 + 2*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MT2, MT2)))/(Pi2*v2) + 
 (0.3333333333333333*(-1 + 4*Cos(2*w))*(pow(Sin(w), 2))*
   Re(B0i(bb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.5*(2 + 3*Cos(2*w))*(pow(Sin(w), 2))*Re(B0i(bb00, 0, MW2, MW2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, ME2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, ML2)))/
  (Pi2*v2) + (0.25*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, 0, MM2)))/
  (Pi2*v2) - (1.*Cos(2*w)*Re(B0i(bb00, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MB2, MT2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.75*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, MD2, MU2)))/(Pi2*v2) + 
 (0.125*(1 - (pow(Cot(w), 2)))*Re(B0i(bb00, MW2, Mh2, MW2)))/(Pi2*v2) - 
 (0.125*(2 + 5*Cos(2*w) + 2*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb00, MW2, MW2, MZ2)))/(Pi2*v2) - 
 (0.375*(pow(Cot(w), 2))*Re(B0i(bb00, MZ2, 0, 0)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MB2, MB2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MC2, MC2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MD2, MD2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*(pow(Cot(w), 2))*Re(B0i(bb00, MZ2, Mh2, MZ2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, ML2, ML2)))/(Pi2*v2) - 
 (0.125*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MM2, MM2)))/(Pi2*v2) - 
 (0.041666666666666664*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MS2, MS2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MT2, MT2)))/(Pi2*v2) - 
 (0.041666666666666664*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.0625*(7 + 8*Cos(2*w) + 3*Cos(4*w))*(pow(Cot(w), 2))*
   Re(B0i(bb00, MZ2, MW2, MW2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MB2, MB2)))/
  (Pi2*v2) - (0.6666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MC2, MC2)))/(Pi2*v2) - 
 (0.16666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MD2, MD2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ME2, ME2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, ML2, ML2)))/
  (Pi2*v2) - (0.5*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MM2, MM2)))/
  (Pi2*v2) - (0.16666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MS2, MS2)))/(Pi2*v2) - 
 (0.6666666666666666*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MT2, MT2)))/
  (Pi2*v2) - (0.6666666666666666*MW2*(pow(Sin(w), 2))*
   Re(B0i(bb1, 0, MU2, MU2)))/(Pi2*v2) - 
 (0.25*MW2*(pow(Sin(w), 2))*Re(B0i(bb1, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.375*(pow(MB, 2))*Re(B0i(bb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*(pow(MC, 2))*Re(B0i(bb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*(pow(MD, 2))*Re(B0i(bb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*(pow(ME, 2))*Re(B0i(bb1, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*(pow(ML, 2))*Re(B0i(bb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*(pow(MM, 2))*Re(B0i(bb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*(pow(MS, 2))*Re(B0i(bb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*MT2*Re(B0i(bb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*(pow(MU, 2))*Re(B0i(bb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*MW2*Re(B0i(bb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*MW2*(pow(Sec(w), 2))*Re(B0i(bb1, Mh2, MZ2, MZ2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ME2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, ML2)))/(Pi2*v2) - 
 (0.125*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, 0, MM2)))/(Pi2*v2) - 
 (0.25*MW2*Cos(2*w)*Re(B0i(bb1, MW2, 0, MW2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MB2, MT2)))/
  (Pi2*v2) - (0.375*MW2*(-1 + (pow(Cot(w), 2)))*
   Re(B0i(bb1, MW2, MC2, MS2)))/(Pi2*v2) + 
 (0.375*MW2*(-1 + (pow(Cot(w), 2)))*Re(B0i(bb1, MW2, MD2, MU2)))/
  (Pi2*v2) + (0.25*MW2*Cos(2*w)*(pow(Cot(w), 2))*
   Re(B0i(bb1, MW2, MW2, MZ2)))/(Pi2*v2) + 
 (0.1875*MW2*(pow(Csc(w), 2))*Re(B0i(bb1, MZ2, 0, 0)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MB2, MB2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MC2, MC2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MD2, MD2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, ME2, ME2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, ML2, ML2)))/(Pi2*v2) + 
 (0.0625*MW2*(2 - 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MM2, MM2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MS2, MS2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MT2, MT2)))/(Pi2*v2) + 
 (0.020833333333333332*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*(pow(Csc(w), 2))*
   Re(B0i(bb1, MZ2, MU2, MU2)))/(Pi2*v2) + 
 (0.25*MW2*(pow(Cos(w), 2))*(pow(Cot(w), 2))*
   Re(B0i(bb1, MZ2, MW2, MW2)))/(Pi2*v2) + 
 (0.75*(pow(MB, 4))*Re(B0i(dbb0, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.75*(pow(MC, 4))*Re(B0i(dbb0, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.75*(pow(MD, 4))*Re(B0i(dbb0, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.25*(pow(ME, 4))*Re(B0i(dbb0, Mh2, ME2, ME2)))/(Pi2*v2) - 
 (0.140625*(pow(Mh, 4))*Re(B0i(dbb0, Mh2, Mh2, Mh2)))/(Pi2*v2) + 
 (0.25*(pow(ML, 4))*Re(B0i(dbb0, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.25*(pow(MM, 4))*Re(B0i(dbb0, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.75*(pow(MS, 4))*Re(B0i(dbb0, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.75*(pow(MT, 4))*Re(B0i(dbb0, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.75*(pow(MU, 4))*Re(B0i(dbb0, Mh2, MU2, MU2)))/(Pi2*v2) + 
 (0.03125*(2*Mh2*MW2 - (pow(Mh, 4)) - 12*(pow(MW, 4)))*
   Re(B0i(dbb0, Mh2, MW2, MW2)))/(Pi2*v2) + 
 (0.015625*(-(pow(Mh, 4)) + 2*Mh2*MW2*(pow(Sec(w), 2)) - 
    12*(pow(MW, 4))*(pow(Sec(w), 4)))*Re(B0i(dbb0, Mh2, MZ2, MZ2)))/
  (Pi2*v2) + (0.1111111111111111*MT2*(12*Alfas*Pi + 
    4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*Re(B0i(dbb0, MT2, 0, MT2)))/
  Pi2 + (0.25*MT2*(pow(MB, 2))*Re(B0i(dbb0, MT2, MB2, MW2)))/(Pi2*v2) - 
 (0.25*(pow(MT, 4))*Re(B0i(dbb0, MT2, Mh2, MT2)))/(Pi2*v2) + 
 (0.013888888888888888*MT2*(9*MT2 - 16*MW2*(1 + 2*Cos(2*w))*
     (pow(Tan(w), 2)))*Re(B0i(dbb0, MT2, MT2, MZ2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MB2, MB2)))/
  (Pi2*v2) + (1.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MC2, MC2)))/(Pi2*v2) + 
 (0.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MD2, MD2)))/
  (Pi2*v2) + (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ME2, ME2)))/
  (Pi2*v2) + (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, ML2, ML2)))/
  (Pi2*v2) + (1.*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MM2, MM2)))/
  (Pi2*v2) + (0.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MS2, MS2)))/(Pi2*v2) + 
 (1.3333333333333333*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MT2, MT2)))/
  (Pi2*v2) + (1.3333333333333333*MW2*(pow(Sin(w), 2))*
   Re(B0i(dbb00, 0, MU2, MU2)))/(Pi2*v2) - 
 (1.5*MW2*(pow(Sin(w), 2))*Re(B0i(dbb00, 0, MW2, MW2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MB, 2))*Re(B0i(dbb1, Mh2, MB2, MB2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MC, 2))*Re(B0i(dbb1, Mh2, MC2, MC2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MD, 2))*Re(B0i(dbb1, Mh2, MD2, MD2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ME, 2))*Re(B0i(dbb1, Mh2, ME2, ME2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(ML, 2))*Re(B0i(dbb1, Mh2, ML2, ML2)))/(Pi2*v2) + 
 (0.125*Mh2*(pow(MM, 2))*Re(B0i(dbb1, Mh2, MM2, MM2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MS, 2))*Re(B0i(dbb1, Mh2, MS2, MS2)))/(Pi2*v2) + 
 (0.375*Mh2*MT2*Re(B0i(dbb1, Mh2, MT2, MT2)))/(Pi2*v2) + 
 (0.375*Mh2*(pow(MU, 2))*Re(B0i(dbb1, Mh2, MU2, MU2)))/(Pi2*v2) - 
 (0.125*Mh2*MW2*Re(B0i(dbb1, Mh2, MW2, MW2)))/(Pi2*v2) - 
 (0.0625*Mh2*MW2*(pow(Sec(w), 2))*Re(B0i(dbb1, Mh2, MZ2, MZ2)))/
  (Pi2*v2) - (0.1111111111111111*MT2*(12*Alfas*Pi + 
    4*MW2*(pow(v, -2))*(pow(Sin(w), 2)))*Re(B0i(dbb1, MT2, 0, MT2)))/
  Pi2 + (0.125*MT2*(MT2 + 2*MW2 + (pow(MB, 2)))*
   Re(B0i(dbb1, MT2, MB2, MW2)))/(Pi2*v2) - 
 (0.125*(pow(MT, 4))*Re(B0i(dbb1, MT2, Mh2, MT2)))/(Pi2*v2) + 
 (0.006944444444444444*MT2*(9*(MT2 + 2*MW2) + (9*MT2 - 8*MW2)*Cos(2*w) + 
    8*MW2*Cos(4*w))*(pow(Sec(w), 2))*Re(B0i(dbb1, MT2, MT2, MZ2)))/
  (Pi2*v2)
;
}

ComplexType dghgaga(double m12, double m22, double m32)
{
return 0. - (0.5*MW2*(Mh2 + 6*MW2)*B0i(bb0, m12, MW2, MW2)*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (2.6666666666666665*MT2*MW2*B0i(bb0, m32, MT2, MT2)*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (1.3333333333333333*(m12 + m22 - m32)*MT2*MW2*C0i(cc0, m22, m32, m12, MT2, 
    MT2, MT2)*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (10.666666666666666*MT2*MW2*C0i(cc00, m22, m32, m12, MT2, MT2, MT2)*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.*MW2*(Mh2 + 6*MW2)*C0i(cc00, m22, m32, m12, MW2, MW2, MW2)*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (5.333333333333333*m22*MT2*MW2*C0i(cc1, m22, m32, m12, MT2, MT2, MT2)*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.6666666666666665*(m12 + m22 - m32)*MT2*MW2*C0i(cc2, m22, m32, m12, MT2, 
    MT2, MT2)*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (0.6666666666666666*MW2*B0i(bb0, m32, MB2, MB2)*(pow(MB, 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (0.3333333333333333*(m12 + m22 - m32)*MW2*C0i(cc0, m22, m32, m12, MB2, MB2, 
    MB2)*(pow(MB, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (2.6666666666666665*MW2*C0i(cc00, m22, m32, m12, MB2, MB2, MB2)*
   (pow(MB, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (1.3333333333333333*m22*MW2*C0i(cc1, m22, m32, m12, MB2, MB2, MB2)*
   (pow(MB, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (0.6666666666666666*(m12 + m22 - m32)*MW2*C0i(cc2, m22, m32, m12, MB2, MB2, 
    MB2)*(pow(MB, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.6666666666666665*MW2*B0i(bb0, m32, MC2, MC2)*(pow(MC, 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (1.3333333333333333*(m12 + m22 - m32)*MW2*C0i(cc0, m22, m32, m12, MC2, MC2, 
    MC2)*(pow(MC, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (10.666666666666666*MW2*C0i(cc00, m22, m32, m12, MC2, MC2, MC2)*
   (pow(MC, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (5.333333333333333*m22*MW2*C0i(cc1, m22, m32, m12, MC2, MC2, MC2)*
   (pow(MC, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.6666666666666665*(m12 + m22 - m32)*MW2*C0i(cc2, m22, m32, m12, MC2, MC2, 
    MC2)*(pow(MC, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (0.6666666666666666*MW2*B0i(bb0, m32, MD2, MD2)*(pow(MD, 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (0.3333333333333333*(m12 + m22 - m32)*MW2*C0i(cc0, m22, m32, m12, MD2, MD2, 
    MD2)*(pow(MD, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (2.6666666666666665*MW2*C0i(cc00, m22, m32, m12, MD2, MD2, MD2)*
   (pow(MD, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (1.3333333333333333*m22*MW2*C0i(cc1, m22, m32, m12, MD2, MD2, MD2)*
   (pow(MD, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (0.6666666666666666*(m12 + m22 - m32)*MW2*C0i(cc2, m22, m32, m12, MD2, MD2, 
    MD2)*(pow(MD, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.*MW2*B0i(bb0, m32, ME2, ME2)*(pow(ME, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (1.*(m12 + m22 - m32)*MW2*
   C0i(cc0, m22, m32, m12, ME2, ME2, ME2)*(pow(ME, 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (8.*MW2*C0i(cc00, m22, m32, m12, ME2, ME2, ME2)*(pow(ME, 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (4.*m22*MW2*C0i(cc1, m22, m32, m12, ME2, ME2, ME2)*(pow(ME, 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.*(m12 + m22 - m32)*MW2*C0i(cc2, m22, m32, m12, ME2, ME2, ME2)*
   (pow(ME, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.*MW2*B0i(bb0, m32, ML2, ML2)*(pow(ML, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (1.*(m12 + m22 - m32)*MW2*
   C0i(cc0, m22, m32, m12, ML2, ML2, ML2)*(pow(ML, 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (8.*MW2*C0i(cc00, m22, m32, m12, ML2, ML2, ML2)*(pow(ML, 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (4.*m22*MW2*C0i(cc1, m22, m32, m12, ML2, ML2, ML2)*(pow(ML, 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.*(m12 + m22 - m32)*MW2*C0i(cc2, m22, m32, m12, ML2, ML2, ML2)*
   (pow(ML, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.*MW2*B0i(bb0, m32, MM2, MM2)*(pow(MM, 2))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (1.*(m12 + m22 - m32)*MW2*
   C0i(cc0, m22, m32, m12, MM2, MM2, MM2)*(pow(MM, 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (8.*MW2*C0i(cc00, m22, m32, m12, MM2, MM2, MM2)*(pow(MM, 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (4.*m22*MW2*C0i(cc1, m22, m32, m12, MM2, MM2, MM2)*(pow(MM, 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.*(m12 + m22 - m32)*MW2*C0i(cc2, m22, m32, m12, MM2, MM2, MM2)*
   (pow(MM, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (0.6666666666666666*MW2*B0i(bb0, m32, MS2, MS2)*(pow(MS, 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (0.3333333333333333*(m12 + m22 - m32)*MW2*C0i(cc0, m22, m32, m12, MS2, MS2, 
    MS2)*(pow(MS, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (2.6666666666666665*MW2*C0i(cc00, m22, m32, m12, MS2, MS2, MS2)*
   (pow(MS, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (1.3333333333333333*m22*MW2*C0i(cc1, m22, m32, m12, MS2, MS2, MS2)*
   (pow(MS, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (0.6666666666666666*(m12 + m22 - m32)*MW2*C0i(cc2, m22, m32, m12, MS2, MS2, 
    MS2)*(pow(MS, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.6666666666666665*MW2*B0i(bb0, m32, MU2, MU2)*(pow(MU, 2))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (1.3333333333333333*(m12 + m22 - m32)*MW2*C0i(cc0, m22, m32, m12, MU2, MU2, 
    MU2)*(pow(MU, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (10.666666666666666*MW2*C0i(cc00, m22, m32, m12, MU2, MU2, MU2)*
   (pow(MU, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (5.333333333333333*m22*MW2*C0i(cc1, m22, m32, m12, MU2, MU2, MU2)*
   (pow(MU, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (2.6666666666666665*(m12 + m22 - m32)*MW2*C0i(cc2, m22, m32, m12, MU2, MU2, 
    MU2)*(pow(MU, 2))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (0.5*B0i(bb0, m22, MW2, MW2)*(pow(MW, 4))*(pow(Sin(w), 2)))/
  (Pi2*(pow(v, 3))) + (0.5*B0i(bb0, m32, MW2, MW2)*(pow(MW, 4))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) - 
 (0.5*(6*m12 - 5*m22 - 5*m32 + Mh2)*C0i(cc0, m22, m32, m12, MW2, MW2, MW2)*
   (pow(MW, 4))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (0.5*(m12 + m22 - m32)*C0i(cc1, m22, m32, m12, MW2, MW2, MW2)*
   (pow(MW, 4))*(pow(Sin(w), 2)))/(Pi2*(pow(v, 3))) + 
 (1.*m12*C0i(cc2, m22, m32, m12, MW2, MW2, MW2)*(pow(MW, 4))*
   (pow(Sin(w), 2)))/(Pi2*(pow(v, 3)))
;
}

ComplexType dghgg(double m12, double m22, double m32)
{
return 0. + (0.6366197723675814*Alfashgg*MT2*B0i(bb0, m32, MT2, MT2))/v + 
 (0.3183098861837907*Alfashgg*(m12 + m22 - m32)*MT2*C0i(cc0, m22, m32, m12, MT2, 
    MT2, MT2))/v - (2.5464790894703255*Alfashgg*MT2*
   C0i(cc00, m22, m32, m12, MT2, MT2, MT2))/v + 
 (1.2732395447351628*Alfashgg*m22*MT2*C0i(cc1, m22, m32, m12, MT2, MT2, MT2))/
  v + (0.6366197723675814*Alfashgg*(m12 + m22 - m32)*MT2*
   C0i(cc2, m22, m32, m12, MT2, MT2, MT2))/v + 
 (0.6366197723675814*Alfashgg*B0i(bb0, m32, MB2, MB2)*(pow(MB, 2)))/v + 
 (0.3183098861837907*Alfashgg*(m12 + m22 - m32)*C0i(cc0, m22, m32, m12, MB2, 
    MB2, MB2)*(pow(MB, 2)))/v - 
 (2.5464790894703255*Alfashgg*C0i(cc00, m22, m32, m12, MB2, MB2, MB2)*
   (pow(MB, 2)))/v + (1.2732395447351628*Alfashgg*m22*
   C0i(cc1, m22, m32, m12, MB2, MB2, MB2)*(pow(MB, 2)))/v + 
 (0.6366197723675814*Alfashgg*(m12 + m22 - m32)*C0i(cc2, m22, m32, m12, MB2, 
    MB2, MB2)*(pow(MB, 2)))/v + 
 (0.6366197723675814*Alfashgg*B0i(bb0, m32, MC2, MC2)*(pow(MC, 2)))/v + 
 (0.3183098861837907*Alfashgg*(m12 + m22 - m32)*C0i(cc0, m22, m32, m12, MC2, 
    MC2, MC2)*(pow(MC, 2)))/v - 
 (2.5464790894703255*Alfashgg*C0i(cc00, m22, m32, m12, MC2, MC2, MC2)*
   (pow(MC, 2)))/v + (1.2732395447351628*Alfashgg*m22*
   C0i(cc1, m22, m32, m12, MC2, MC2, MC2)*(pow(MC, 2)))/v + 
 (0.6366197723675814*Alfashgg*(m12 + m22 - m32)*C0i(cc2, m22, m32, m12, MC2, 
    MC2, MC2)*(pow(MC, 2)))/v + 
 (0.6366197723675814*Alfashgg*B0i(bb0, m32, MD2, MD2)*(pow(MD, 2)))/v + 
 (0.3183098861837907*Alfashgg*(m12 + m22 - m32)*C0i(cc0, m22, m32, m12, MD2, 
    MD2, MD2)*(pow(MD, 2)))/v - 
 (2.5464790894703255*Alfashgg*C0i(cc00, m22, m32, m12, MD2, MD2, MD2)*
   (pow(MD, 2)))/v + (1.2732395447351628*Alfashgg*m22*
   C0i(cc1, m22, m32, m12, MD2, MD2, MD2)*(pow(MD, 2)))/v + 
 (0.6366197723675814*Alfashgg*(m12 + m22 - m32)*C0i(cc2, m22, m32, m12, MD2, 
    MD2, MD2)*(pow(MD, 2)))/v + 
 (0.6366197723675814*Alfashgg*B0i(bb0, m32, MS2, MS2)*(pow(MS, 2)))/v + 
 (0.3183098861837907*Alfashgg*(m12 + m22 - m32)*C0i(cc0, m22, m32, m12, MS2, 
    MS2, MS2)*(pow(MS, 2)))/v - 
 (2.5464790894703255*Alfashgg*C0i(cc00, m22, m32, m12, MS2, MS2, MS2)*
   (pow(MS, 2)))/v + (1.2732395447351628*Alfashgg*m22*
   C0i(cc1, m22, m32, m12, MS2, MS2, MS2)*(pow(MS, 2)))/v + 
 (0.6366197723675814*Alfashgg*(m12 + m22 - m32)*C0i(cc2, m22, m32, m12, MS2, 
    MS2, MS2)*(pow(MS, 2)))/v + 
 (0.6366197723675814*Alfashgg*B0i(bb0, m32, MU2, MU2)*(pow(MU, 2)))/v + 
 (0.3183098861837907*Alfashgg*(m12 + m22 - m32)*C0i(cc0, m22, m32, m12, MU2, 
    MU2, MU2)*(pow(MU, 2)))/v - 
 (2.5464790894703255*Alfashgg*C0i(cc00, m22, m32, m12, MU2, MU2, MU2)*
   (pow(MU, 2)))/v + (1.2732395447351628*Alfashgg*m22*
   C0i(cc1, m22, m32, m12, MU2, MU2, MU2)*(pow(MU, 2)))/v + 
 (0.6366197723675814*Alfashgg*(m12 + m22 - m32)*C0i(cc2, m22, m32, m12, MU2, 
    MU2, MU2)*(pow(MU, 2)))/v
;
}

} //end namespace SM
