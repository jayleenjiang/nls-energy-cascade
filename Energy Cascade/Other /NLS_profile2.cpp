#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define NDEBUG
#define EIGEN_NO_DEBUG
#define LAPACK_COMPLEX_CUSTOM
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <cstdlib>
#include <float.h>
#include <random>
#include <stdint.h>
#include <climits>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>

#define _USE_MATH_DEFINES
using namespace std;
using namespace Eigen;


double dt = 0.0002;
double gam = 0.02;
double T1 = 1.0;
double T3 = 2.0;



VectorXd fun(const VectorXd& X)
//X = I1, I2, I3, theta1, theta3
{
    VectorXd tmp(5);
    double I1 = X(0);
    double I2 = X(1);
    double I3 = X(2);
    double theta1 = X(3);
    double theta3 = X(4);
    
    tmp(0) = 4*I2*I1*( sin(theta1) - gam*( 1 + cos(theta1) ) ) + 4*gam*T1 - 2*gam*(I1*I1 + 2*I3*I1);
    tmp(1) = -4*I2*( I1*sin(theta1) + I3*sin(theta3) );
    tmp(2) = 4*I2*I3*( sin(theta3) - gam*(1 + cos(theta3) ) ) + 4*gam*T3 - 2*gam*(I3*I3 + 2*I3*I1);
    tmp(3) = I2*( 1 + 2*cos(theta1) + 2*gam*sin(theta1) ) - I1 - 2*I1*cos(theta1) - 2*I3*cos(theta3);
    tmp(4) = I2*( 1 + 2*cos(theta3) + 2*gam*sin(theta3) ) - I3 - 2*I3*cos(theta3) - 2*I1*cos(theta1);
    return tmp;
}

VectorXd fun2(const VectorXd& X, const Vector4d& rnd)
//X = I1, I2, I3, theta1,  theta3
{
    VectorXd tmp(5);
    double I1 = X(0);
    double I2 = X(1);
    double I3 = X(2);
    tmp(0) = 2*sqrt(2*gam*T1*I1)*rnd(0) + 2*gam*T1*(pow(rnd(0), 2) - dt);
    tmp(1) = 0;
    tmp(2) = 2*sqrt(2*gam*T3*I3)*rnd(1) + 2*gam*T3*(pow(rnd(1), 2) - dt);
    tmp(3) = sqrt(2*gam*T1)*pow(I1, -0.5)*rnd(2);
    tmp(4) = sqrt(2*gam*T3)*pow(I3, -0.5)*rnd(3);
    return tmp;
}



int main()
{
    struct timeval t1, t2;
    gettimeofday(&t1,NULL);
    random_device rd;
    mt19937 mt(rd());
    VectorXd X0(5);
    uniform_real_distribution<double> u(0.0,1.0);
    X0 << 1, 0, 1, 0.0, 0.0;
    int N = 5000000;
    MatrixXd Traj(N + 1, 5);
    for(int j = 0; j < 5; j++)
        Traj(0, j) = X0(j);
    ofstream myfile;
    myfile.open("traj.txt");
    normal_distribution<double> nm(0.0, 1.0);
    for(int i = 0; i < N; i++)
    {
        Vector4d rnd;
        rnd << nm(mt), nm(mt), nm(mt), nm(mt);
        rnd *= sqrt(dt);
        //cout<<fun2(X0, rnd).transpose()<<endl;
        X0 += dt*fun(X0) + fun2(X0, rnd);
        if(X0(0) < 0)
        {
            X0(0)*= -1;
            cout<<"reflection!"<<endl;
        }
        if(X0(2) < 0)
        {
            X0(2)*= -1;
            cout<<"reflection!"<<endl;
        }
        for(int j = 0; j < 5; j++)
            Traj(i+1, j) = X0(j);
        if(X0(1) > 1.0)
        {
            N = i;
            break;
        }
        //cout<<X0.transpose()<<endl;
    }
    cout<<"stop at step "<<N<<endl;
    for(int i = 0; i < N + 1; i++)
    {
        for(int j = 0; j < 5; j++)
            myfile<<Traj(i,j)<<" ";
        myfile<<endl;
    }
    cout<<X0.transpose()<<endl;
    gettimeofday(&t2, NULL);
    double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u +
                    t2.tv_usec - t1.tv_usec) / 1.e6;
    
    cout << "total CPU time = " << delta <<endl;
}
