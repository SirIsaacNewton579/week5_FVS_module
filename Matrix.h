#include<iostream>
using namespace std;
void mul_matrix(double *A,double *B,double *ret,int m=3,int o=3,int n=3){
    //矩阵乘法
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            ret[n*i+j] = 0.0;
            for(int k=0;k<o;k++){
                ret[n*i+j] += A[i*o+k]*B[k*n+j];
            }
        }
    }
}
void mul_matrix(double *A,double a,int m=1,int n=3){
    //矩阵乘标量
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            A[n*i+j] *= a;
        }
    }
}
void add_matrix(double *A,double *B,double *ret,int m=1,int n=3){
    //矩阵加法
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            ret[n*i+j] = A[n*i+j]+B[n*i+j];
        }
    }
}
void minus_matrix(double *A,double *B,double *ret,int m=1,int n=3){
    //矩阵减法
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            ret[n*i+j] = A[n*i+j]-B[n*i+j];
        }
    }
}
void print_matrix(double *p,int m=3,int n=3){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            cout << p[i*m+j] << "\t";
        }
        cout << endl;
    }
}