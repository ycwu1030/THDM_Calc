#include <iostream>
#include <fstream>
#include <vector>
#include "ModelParameters.h"
#include "clooptools.h"
#include "CouplingFunctionSM.h"
#include "CouplingFunctionTypeI.h"
#include "CouplingFunctionTypeII.h"
#include "CouplingFunctionTypeLS.h"
#include "CouplingFunctionTypeFL.h"
#include "CouplingFunctionSTU.h"
#include "HSS_STU_ExpData.h"

using namespace std;
using namespace TypeII;

double RandomReal(double min, double max)
{
    return (rand()/(double)RAND_MAX)*(max-min)+min;
}
int main(int argc, char const *argv[])
{
    ofstream output("THDM_SCAN_TYPEII.dat");
    ostream_iterator<string> oss(output,"    ");
    ostream_iterator<double> od(output,"    ");
    ltini();
    ComplexType dkappab,dkappac,dkappatau;
    ComplexType dkappag,dkappaga,dkappaW,dkappaZ;
    ComplexType SMkappab = SM::dKappahbb(Mh2,MB2,MB2);
    ComplexType SMkappac = SM::dKappahcc(Mh2,MC2,MC2);
    ComplexType SMkappatau = SM::dKappahtautau(Mh2,ML2,ML2);
    ComplexType SMkappaW = SM::dKappahWW(Mh2,30*30,30*30);
    ComplexType SMkappaZ = SM::dKappahZZ(Mh2,30*30,30*30);
    ComplexType SMkappag = SM::dghgg(Mh2,0.01,0.01);
    ComplexType SMkappaga = SM::dghgaga(Mh2,0.01,0.01);

    double kappau, kappad, kappac, kappas, kappat, kappab, kappaele, kappatau, kappamu, kappaW, kappaZ, kappag, kappaga, kappaZga;
    

    double MHH,MHA,MHp;
    double beta,alpha;
    double M2;
    double tb;
    double m122;

    double S,T,U;

    MHH = 80;
    MHA = 400;
    MHp = 400;
    int GOT = 0;
    double logtb;
    double MHH2,MHA2,MHp2;
    int DOFLHC8,DOFLHC13,DOFLHC3000,DOFCEPC;
    double chi2muLHC8, chi2muLHC13, chi2muLHC3000, chi2muCEPC;
    double chi2STULHC8, chi2STULHC13, chi2STULHC3000, chi2STUCEPC;
    bool goodLHC8, goodLHC13, goodLHC3000, goodCEPC;
    vector<string> title = {"GOT", "MHH", "MHA", "MHp", "m122", "beta", "alpha", "kb", "kc", "ktau", "kW", "kZ", "kg", "kga", "S", "T", "U", "chi2LHC8", "DOFLHC8", "GOODLHC8", "chi2LHC13", "DOFLHC13", "GOODLHC13","chi2LHC3000","DOFLHC3000","GOODLHC3000","chi2CEPC","DOFCEPC","GOODCEPC"};
    for (auto e :title)
    {
        oss = e;
    }
    output<<endl;
    while (GOT<10000)
    {
        MHH2 = MHH*MHH;
        MHA2 = MHA*MHA;
        MHp2= MHp*MHp;
        logtb=RandomReal(-1,1.7);
        tb = pow(10,logtb);
        m122 = RandomReal(0,10000);
        beta = atan(tb);
        alpha = beta - Pi/2.0;
        M2 = m122/sin(beta)/cos(beta);
        
        if ((!Check_Unitarity(MHH,MHA,MHp,M2,beta,alpha))||(!Check_Stability(MHH,MHA,MHp,M2,beta,alpha)))
        {
            continue;
        }
        
        dkappab = dKappahbb(MHH2,MHA2,MHp2,M2,beta,alpha,Mh2,MB2,MB2);
        dkappac = dKappahcc(MHH2,MHA2,MHp2,M2,beta,alpha,Mh2,MC2,MC2);
        dkappatau = dKappahtautau(MHH2,MHA2,MHp2,M2,beta,alpha,Mh2,ML2,ML2);
        dkappag = dghgg(MHH2,MHA2,MHp2,M2,beta,alpha,Mh2,0.01,0.01);
        dkappaga = dghgaga(MHH2,MHA2,MHp2,M2,beta,alpha,Mh2,0.01,0.01);
        dkappaW = dKappahWW(MHH2,MHA2,MHp2,M2,beta,alpha,Mh2,30*30,30*30);
        dkappaZ = dKappahZZ(MHH2,MHA2,MHp2,M2,beta,alpha,Mh2,30*30,30*30);

        S = Re(STU::SinTHDM(MHH2,MHA2,MHp2,M2,beta,alpha));
        T = Re(STU::TinTHDM(MHH2,MHA2,MHp2,M2,beta,alpha));
        U = Re(STU::UinTHDM(MHH2,MHA2,MHp2,M2,beta,alpha));

        kappau = 1.0;
        kappad = 1.0;
        kappac = Re(dkappac/SMkappac)+1.0;
        kappas = 1.0;
        kappat = 1.0;
        kappab = Re(dkappab/SMkappab)+1.0;
        kappatau = Re(dkappatau/SMkappatau)+1.0;
        kappamu = 1.0;
        kappaele = 1.0;
        kappaW = Re(dkappaW/SMkappaW)+1.0;
        kappaZ = Re(dkappaZ/SMkappaZ)+1.0;
        kappag = Re(dkappag/SMkappag)+1.0;
        kappaga = Re(dkappaga/SMkappaga)+1.0;
        kappaZga = 1.0;

        KAPPAS kappainput = {kappau,kappad,kappac,kappas,kappat,kappab,kappaele,kappamu,kappatau,kappaW,kappaZ,kappag,kappaga,kappaZga};

        HiggsSignalStrengthSTU_Test(muLHC8,STULHC,kappainput,S,T,U,chi2muLHC8,chi2STULHC8,DOFLHC8,goodLHC8);
        HiggsSignalStrengthSTU_Test(muATLAS13+muCMS13,STULHC,kappainput,S,T,U,chi2muLHC13,chi2STULHC13,DOFLHC13,goodLHC13);
        HiggsSignalStrengthSTU_Test(muHLLHC3000,STULHC,kappainput,S,T,U,chi2muLHC3000,chi2STULHC3000,DOFLHC3000,goodLHC3000);
        HiggsSignalStrengthSTU_Test(muCEPC,STUCEPC,kappainput,S,T,U,chi2muCEPC,chi2STUCEPC,DOFCEPC,goodCEPC);

        od = GOT;
        od = MHH;
        od = MHA;
        od = MHp;
        od = m122;
        od = beta;
        od = alpha;
        od = kappab;
        od = kappac;
        od = kappatau;
        od = kappaW;
        od = kappaZ;
        od = kappag;
        od = kappaga;
        od = S;
        od = T;
        od = U;
        od = chi2muLHC8 + chi2STULHC8;
        od = DOFLHC8;
        od = goodLHC8;
        od = chi2muLHC13 + chi2STULHC13;
        od = DOFLHC13;
        od = goodLHC13;
        od = chi2muLHC3000 + chi2STULHC3000;
        od = DOFLHC3000;
        od = goodLHC3000;
        od = chi2muCEPC + chi2STUCEPC;
        od = DOFCEPC;
        od = goodCEPC;
        output << endl;
        ++GOT;
    }
    


    return 0;
}
