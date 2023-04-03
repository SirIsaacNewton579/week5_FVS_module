#include<iostream>
#include<fstream>
#include<math.h>
#include "gamma.h"   //声明全局变量gamma as g
#include "Matrix.h"  //矩阵运算
#include "FVS.h"  //通量矢量分裂相关函数 以及 差分格式
using namespace std;

double g = 1.4;
int main(){
    //开始计算
    int Nx = 101;
    double dx = 1./(Nx-1),dt = 0.001;
    double t_end = 0.14,Nt = round(t_end/dt);
    cout << "Nt=" << Nt << endl;
    double U[3][Nx];
    
    for(int i=0;i<Nx;i++){
        //cout << i*dx << endl;
        U[0][i] = (i>Nx/2 ? 0.125 : 1);  //rho
        U[1][i] = 0.;  // rho*u
        U[2][i] = (i>Nx/2 ? 0.1 : 1)/(g-1); //E = 1/2*rho*u^2 + p/(g-1)
    }
    updateU(U[0],Nx,dx,dt,Nt,GVC2,3,1,LFf);  
    //输出
    ofstream csvfile;
    csvfile.open("LF_GVC2_t=0.14.csv", ios::out | ios::trunc);
    csvfile <<"x" << "," << "rho"<<","<< "u" <<","<<"p"<< endl;
    for(int i=0;i<Nx;i++){
        csvfile <<i*dx <<","<< U[0][i]<<","<<U[1][i]/U[0][i]<<","<<(g-1)*(U[2][i]-0.5*U[1][i]*U[1][i]/U[0][i]) << endl;
    }
    csvfile.close();
    //system("pause");
    return 0;
}