#include "ModelParameters.h"

double lam1(double MHH, double MHA, double MHp, double M2, double beta, double alpha)
{
    double _MHH2 = MHH*MHH;
    double _MHL2 = Mh2;
    double _sb = sin(beta);
    double _cb = cos(beta);
    double _m122 = M2*_sb*_cb;
    double _cba = cos(beta-alpha);
    double _sba = sin(beta-alpha);
    double _vev2 = v2;
    double l1 = (_cb*((_MHL2-_MHH2)*pow(_cba*_sb-_cb*_sba,2)+_MHH2)-_m122*_sb)/pow(_cb,3)/_vev2;
    return l1;
}
double lam2(double MHH, double MHA, double MHp, double M2, double beta, double alpha)
{
    double _MHH2 = MHH*MHH;
    double _MHL2 = Mh2;
    double _sb = sin(beta);
    double _cb = cos(beta);
    double _m122 = M2*_sb*_cb;
    double _cba = cos(beta-alpha);
    double _sba = sin(beta-alpha);
    double _vev2 = v2;
    double l2 = (_sb*((_MHH2-_MHL2)*pow(_cba*_sb-_cb*_sba,2)+_MHL2)-_m122*_cb)/pow(_sb,3)/_vev2;
    return l2;
}
double lam3(double MHH, double MHA, double MHp, double M2, double beta, double alpha)
{
    double _MHH2 = MHH*MHH;
    double _MHL2 = Mh2;
    double _MHpm2 = MHp*MHp;
    double _sb = sin(beta);
    double _cb = cos(beta);
    double _m122 = M2*_sb*_cb;
    double _cba = cos(beta-alpha);
    double _sba = sin(beta-alpha);
    double _vev2 = v2;
    double l3 = ((_MHH2-_MHL2)*(_cba*_sb-_cb*_sba)*(_cb*_cba+_sb*_sba)+2.0*_cb*_sb*_MHpm2-_m122)/_cb/_sb/_vev2;
    return l3;
}
double lam4(double MHH, double MHA, double MHp, double M2, double beta, double alpha)
{
    double _MHA2 = MHA*MHA;
    double _MHpm2 = MHp*MHp;
    double _sb = sin(beta);
    double _cb = cos(beta);
    double _m122 = M2*_sb*_cb;
    double _cba = cos(beta-alpha);
    double _sba = sin(beta-alpha);
    double _vev2 = v2;
    double l4 = (_cb*_sb*(_MHA2-2*_MHpm2)+_m122)/_cb/_sb/_vev2;
    return l4;
}
double lam5(double MHH, double MHA, double MHp, double M2, double beta, double alpha)
{
    double _MHA2 = MHA*MHA;
    double _MHpm2 = MHp*MHp;
    double _sb = sin(beta);
    double _cb = cos(beta);
    double _m122 = M2*_sb*_cb;
    double _cba = cos(beta-alpha);
    double _sba = sin(beta-alpha);
    double _vev2 = v2;
    double l5 = (_m122-_cb*_sb*_MHA2)/_sb/_cb/_vev2;
    return l5;
}

bool Check_Unitarity(double MHH, double MHA, double MHp, double M2, double beta, double alpha)
{
    double MAX = 0.5;
    double _lam1 = lam1(MHH,MHA,MHp,M2,beta,alpha);
    double _lam2 = lam2(MHH,MHA,MHp,M2,beta,alpha);
    double _lam3 = lam3(MHH,MHA,MHp,M2,beta,alpha);
    double _lam4 = lam4(MHH,MHA,MHp,M2,beta,alpha);
    double _lam5 = lam5(MHH,MHA,MHp,M2,beta,alpha);
    double a0Eigens[12];
    a0Eigens[0] = _lam3+_lam4;
    a0Eigens[1] = _lam3-_lam4;
    a0Eigens[2] = _lam3+_lam5;
    a0Eigens[3] = _lam3-_lam5;
    a0Eigens[4] = _lam3+2*_lam4+3*_lam5;
    a0Eigens[5] = _lam3+2*_lam4-3*_lam5;
    a0Eigens[6] = (_lam1+_lam2+sqrt(pow(_lam1-_lam2,2)+4*pow(_lam5,2)))/2;
    a0Eigens[7] = (_lam1+_lam2-sqrt(pow(_lam1-_lam2,2)+4*pow(_lam5,2)))/2;
    a0Eigens[8] = (_lam1+_lam2+sqrt(pow(_lam1-_lam2,2)+4*pow(_lam4,2)))/2;
    a0Eigens[9] = (_lam1+_lam2-sqrt(pow(_lam1-_lam2,2)+4*pow(_lam4,2)))/2;
    a0Eigens[10] = (3*(_lam1+_lam2)+sqrt(9*pow(_lam1-_lam2,2)+4*pow(2*_lam3+_lam4,2)))/2;
    a0Eigens[11] = (3*(_lam1+_lam2)-sqrt(9*pow(_lam1-_lam2,2)+4*pow(2*_lam3+_lam4,2)))/2;
    bool good = true;
    for (int i = 0; i < 12; i++)
    {
        good*=(abs(a0Eigens[i])<MAX*16.0*Pi);
        if (!good)
        {
            return good;
        }
    }
    return good;
}

bool Check_Stability(double MHH, double MHA, double MHp, double M2, double beta, double alpha)
{
    double _lam1 = lam1(MHH,MHA,MHp,M2,beta,alpha);
    double _lam2 = lam2(MHH,MHA,MHp,M2,beta,alpha);
    double _lam3 = lam3(MHH,MHA,MHp,M2,beta,alpha);
    double _lam4 = lam4(MHH,MHA,MHp,M2,beta,alpha);
    double _lam5 = lam5(MHH,MHA,MHp,M2,beta,alpha);
    if (_lam1 <= 0 || _lam2 <= 0)
    {
        return false;
    }
    
    if (_lam3 <= -sqrt(_lam1*_lam2) || _lam3+_lam4-_lam5 <= -sqrt(_lam1*_lam2))
    {
        return false;
    }

    return true;
}
// THDMParameter::THDMParameter()
// {
//     MH=500;
//     MA0=500;
//     MHp=500;
//     M2=500*500;
//     beta=atan(1);
//     alp=beta-PiHalf;

//     MH2=MH*MH;
//     MA02=MA0*MA0;
//     MHp2=MHp*MHp;
// }