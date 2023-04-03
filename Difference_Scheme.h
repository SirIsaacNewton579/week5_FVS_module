#include<math.h>
using namespace std;
void first_windward(double *f,double &fm,int pm) {fm = f[0];}
void VanLeerLim(double *ft,double &ftm,int pm){
    //pm = 1左值，pm=-1右值
    int j=1;
    double dL = ft[j]-ft[j-1*pm],dR=ft[j+1*pm]-ft[j];
    double eps = 1e-12;
    double phi = (abs(dL*dR) + dL*dR+eps)/(abs(dL*dR)+dR*dR+eps);
    ftm = ft[j]+0.5*phi*dR;
}
double minmod(double a,double b){
    if(!((a>0.0)^(b>0.0))) return (abs(a)>abs(b)?b:a);
    else return 0.0;
}
void NND(double *ft,double &ftm,int pm){
    //pm = 1左值，pm=-1右值
    int j=1;
    ftm = ft[j] + 0.5*minmod(ft[j]-ft[j-1*pm],ft[j+1*pm]-ft[j]);
}
void GVC2(double *ft,double &ftm,int pm){
    int j=1;
    if(abs(ft[j]-ft[j-1*pm]) < abs(ft[j+1*pm]-ft[j])){
        ftm = 0.5*(3*ft[j]-ft[j-1*pm]);
    }
    else{
        ftm = 0.5*(ft[j] + ft[j+1*pm]);
    }
}
void MUSCL(double *ft,double &ftm,int pm){
    int j=1;
    double dmf=ft[j]-ft[j-1*pm],dpf=ft[j+1*pm]-ft[j];
    double eps = 1e-6;
    double s = (2*dmf*dpf+eps)/(dmf*dmf+dpf*dpf+eps);
    ftm = ft[j]+s/4.0*((1-s/3.0)*dmf+(1+s/3.0)*dpf);
}