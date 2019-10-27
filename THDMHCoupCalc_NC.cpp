#include <iostream>
#include <fstream>
#include "ModelParameters.h"
#include "clooptools.h"
#include "CouplingFunctionSM.h"
#include "CouplingFunctionTypeI.h"
#include "CouplingFunctionTypeII.h"
#include "CouplingFunctionTypeLS.h"
#include "CouplingFunctionTypeFL.h"

using namespace std;
using namespace TypeFL;

int main()
{
    //
    ofstream outfile("TypeFL_mamc_AlignMphi1000Msft100tb1.dat");
    ltini();
    //RealType mphi;
    ComplexType cb,cc,ctau;
    ComplexType cg,cga,cW,cZ;
    ComplexType SMkappab = SM::dKappahbb(Mh2,MB2,MB2);
    ComplexType SMkappac = SM::dKappahcc(Mh2,MC2,MC2);
    ComplexType SMkappatau = SM::dKappahtautau(Mh2,ML2,ML2);
    ComplexType SMkappaW = SM::dKappahWW(Mh2,30,30);
    ComplexType SMkappaZ = SM::dKappahZZ(Mh2,30,30);
    ComplexType SMkappag = SM::dghgg(Mh2,0.01,0.01);
    ComplexType SMkappaga = SM::dghgaga(Mh2,0.01,0.01);
    
    
    double mphi=1000.0;
    double M2=100*100;
    //
    double delZh=0.51/100;
    double delbb=0.28/100;
    double delcc=2.2/100;
    double delgg=1.6/100;
    double delWW=1.5/100;
    double deltautau=1.2/100;
    double delZZ=4.3/100;
    double delgaga=9.0/100;
    double delmumu=0.17;
    double delVBFbb=2.8/100;
    
    double brb=0.578;
    double brtau=0.0637;
    double brmu=0.000221;
    double brs=0.00044;
    double brc=0.0268;
    double brg=0.0856;
    double brga=0.0023;
    double brW=0.216;
    double brZ=0.0267;
    
    double cba=0;
    double tb=1.0;
    
    for ( double ma=-300; ma<=300; ma=ma+5 )
        for ( double mc=-300; mc<=300; mc=mc+5 )
        {
	  cb=dKappahbb(mphi*mphi,(mphi+ma)*(mphi+ma),(mphi+mc)*(mphi+mc),mphi*mphi-M2,atan(tb),atan(tb)-acos(cba),Mh2,MB2,MB2)/SMkappab;
          cc=dKappahcc(mphi*mphi,(mphi+ma)*(mphi+ma),(mphi+mc)*(mphi+mc),mphi*mphi-M2,atan(tb),atan(tb)-acos(cba),Mh2,MC2,MC2)/SMkappac;
          ctau=dKappahtautau(mphi*mphi,(mphi+ma)*(mphi+ma),(mphi+mc)*(mphi+mc),mphi*mphi-M2,atan(tb),atan(tb)-acos(cba),Mh2,ML2,ML2)/SMkappatau;
	  //
          cg=dghgg(mphi*mphi,(mphi+ma)*(mphi+ma),(mphi+mc)*(mphi+mc),mphi*mphi-M2,atan(tb),atan(tb)-acos(cba),Mh2,0.01,0.01)/SMkappag;
          cga=dghgaga(mphi*mphi,(mphi+ma)*(mphi+ma),(mphi+mc)*(mphi+mc),mphi*mphi-M2,atan(tb),atan(tb)-acos(cba),Mh2,0.01,0.01)/SMkappaga;
          cW=dKappahWW(mphi*mphi,(mphi+ma)*(mphi+ma),(mphi+mc)*(mphi+mc),mphi*mphi-M2,atan(tb),atan(tb)-acos(cba),Mh2,30,30)/SMkappaW;
          cZ=dKappahZZ(mphi*mphi,(mphi+ma)*(mphi+ma),(mphi+mc)*(mphi+mc),mphi*mphi-M2,atan(tb),atan(tb)-acos(cba),Mh2,30,30)/SMkappaZ;
	  //
            
            double ch=0;
            ch=Re(cb)*Re(cb)*brb+Re(cc)*Re(cc)*brc+Re(ctau)*Re(ctau)*brtau+Re(cg)*Re(cg)*brg+Re(cga)*Re(cga)*brga+Re(cW)*Re(cW)*brW+Re(cZ)*Re(cZ)*brZ;
            
            double chisq=0;
            chisq=(Re(cZ)*Re(cZ)-1)*(Re(cZ)*Re(cZ)-1)/delZh/delZh;
            chisq=chisq+pow(Re(cZ)*Re(cZ)*Re(cb)*Re(cb)/ch-1,2)/delbb/delbb;
            chisq=chisq+pow(Re(cZ)*Re(cZ)*Re(cc)*Re(cc)/ch-1,2)/delcc/delcc;
            chisq=chisq+pow(Re(cZ)*Re(cZ)*Re(cg)*Re(cg)/ch-1,2)/delgg/delgg;
            chisq=chisq+pow(Re(cZ)*Re(cZ)*Re(cW)*Re(cW)/ch-1,2)/delWW/delWW;
            chisq=chisq+pow(Re(cZ)*Re(cZ)*Re(ctau)*Re(ctau)/ch-1,2)/deltautau/deltautau;
            chisq=chisq+pow(Re(cZ)*Re(cZ)*Re(cZ)*Re(cZ)/ch-1,2)/delZZ/delZZ;
            chisq=chisq+pow(Re(cZ)*Re(cZ)*Re(cga)*Re(cga)/ch-1,2)/delgaga/delgaga;
            chisq=chisq+pow(Re(cZ)*Re(cZ)*Re(ctau)*Re(ctau)/ch-1,2)/delmumu/delmumu;
            chisq=chisq+pow(Re(cW)*Re(cW)*Re(cb)*Re(cb)/ch-1,2)/delVBFbb/delVBFbb;
            
            outfile<<ma<<" "<<mc<<" "<<Re(cb)<<" "<<Re(cc)<<" "<<Re(ctau)<<" "<<Re(cg)<<" "<<Re(cga)<<" "<<Re(cW)<<" "<<Re(cZ)<<endl;
        }
    
    
    /*
     double mphi=1000.0;
     double M2=300*300;
     for ( double cba = -1.0; cba <=1.0; cba = cba+0.01 )
     for ( double tb = 0.1; tb<=20.0; tb = tb+0.1 )
     {
     cb=TypeFL::dKappahbb(mphi*mphi,mphi*mphi,mphi*mphi,mphi*mphi-M2,atan(tb),atan(tb)-acos(cba),Mh2,MB2,MB2)/SMkappab;
     cc=TypeFL::dKappahcc(mphi*mphi,mphi*mphi,mphi*mphi,mphi*mphi-M2,atan(tb),atan(tb)-acos(cba),Mh2,MC2,MC2)/SMkappac;
     ctau=TypeFL::dKappahcc(mphi*mphi,mphi*mphi,mphi*mphi,mphi*mphi-M2,atan(tb),atan(tb)-acos(cba),Mh2,ML2,ML2)/SMkappatau;
     //
     cg=TypeFL::dghgg(mphi*mphi,mphi*mphi,mphi*mphi,mphi*mphi-M2,atan(tb),atan(tb)-acos(cba),Mh2,0.01,0.01)/SMkappag;
     cga=TypeFL::dghgaga(mphi*mphi,mphi*mphi,mphi*mphi,mphi*mphi-M2,atan(1),atan(1)-acos(cba),Mh2,0.01,0.01)/SMkappaga;
     cW=TypeFL::dKappahWW(mphi*mphi,mphi*mphi,mphi*mphi,mphi*mphi-M2,atan(tb),atan(tb)-acos(cba),Mh2,30,30)/SMkappaW;
     cZ=TypeFL::dKappahZZ(mphi*mphi,mphi*mphi,mphi*mphi,mphi*mphi-M2,atan(tb),atan(tb)-acos(cba),Mh2,30,30)/SMkappaZ;
     
     outfile<<cba<<" "<<tb<<" "<<Re(cb)<<" "<<Re(cc)<<" "<<Re(ctau)<<" "<<Re(cg)<<" "<<Re(cga)<<" "<<Re(cW)<<" "<<Re(cZ)<<endl;
     }
     */
    
    
    
    /*
     for (double i = log10(300); i < 4; i=i+0.01)
     {
     mphi = pow(10,i);
     // Parameters are MHH2, MA02, MHp2, M2, beta, alp, m12, m22, m32
     cb=TypeFL::dKappahbb(mphi*mphi,mphi*mphi+100*100,mphi*mphi+100*100,mphi*mphi-300*300,atan(1),atan(1)-PiHalf,Mh2,MB2,MB2)/SMkappab;
     cc=TypeFL::dKappahcc(mphi*mphi,mphi*mphi+100*100,mphi*mphi+100*100,mphi*mphi-300*300,atan(1),atan(1)-PiHalf,Mh2,MC2,MC2)/SMkappac;
     ctau=TypeFL::dKappahtautau(mphi*mphi,mphi*mphi+100*100,mphi*mphi+100*100,mphi*mphi-300*300,atan(1),atan(1)-PiHalf,Mh2,ML2,ML2)/SMkappatau;
     //
     cg=TypeFL::dghgg(mphi*mphi,mphi*mphi+100*100,mphi*mphi+100*100,mphi*mphi-300*300,atan(1),atan(1)-PiHalf,Mh2,0.01,0.01)/SMkappag;
     cga=TypeFL::dghgaga(mphi*mphi,mphi*mphi+100*100,mphi*mphi+100*100,mphi*mphi-300*300,atan(1),atan(1)-PiHalf,Mh2,0.01,0.01)/SMkappaga;
     cW=TypeFL::dKappahWW(mphi*mphi,mphi*mphi+100*100,mphi*mphi+100*100,mphi*mphi-300*300,atan(1),atan(1)-PiHalf,Mh2,30,30)/SMkappaW;
     cZ=TypeFL::dKappahZZ(mphi*mphi,mphi*mphi+100*100,mphi*mphi+100*100,mphi*mphi-300*300,atan(1),atan(1)-PiHalf,Mh2,30,30)/SMkappaZ;
     //
     outfile<<i<<" "<<Re(cb)<<" "<<Re(cc)<<" "<<Re(ctau)<<" "<<Re(cg)<<" "<<Re(cga)<<" "<<Re(cW)<<" "<<Re(cZ)<<endl;
     }
     */
    
    ltexi();
    
    outfile.close();
    return 0;
}

