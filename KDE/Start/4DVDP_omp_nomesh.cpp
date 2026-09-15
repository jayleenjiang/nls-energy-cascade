#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define NDEBUG
#define EIGEN_NO_DEBUG
#define LAPACK_COMPLEX_CUSTOM
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <random>
#include <complex>
#include <omp.h>
double lowx = -3.0;
double Sp = 6.0;
int N = 300;
int N_box = 20000;
long long int N_sample = 3e9;
int Burn_in = 1000000;
double dt = 0.001;
double h = Sp/N;
double ratio = 1.0;
double mu = 1.5;
double delta = 0;//coupling strength
int gap = 1;//number of steps between samples


using namespace Eigen;
using namespace std;



Vector4d VDP4(Vector4d& X)
{
    Vector4d tmp;
    tmp(0) = mu*(X(0) - pow(X(0), 3)/3.0 - X(1)) + delta*(X(0) - X(2));
    tmp(1) = X(0)/mu;
    tmp(2) = mu*(X(2) - pow(X(2), 3)/3.0 - X(3)) + delta*(X(2) - X(0));
    tmp(3) = X(2)/mu;
    return tmp;
}

int which_box(vector<vector<int>>& list_of_box, Vector4d& xx)
{
    int ret_val = -1;
    //cout<<"xx = "<<xx.transpose()<<endl;
    if(xx(0) >= lowx && xx(0) < lowx + Sp && xx(1) >= lowx && xx(1) < lowx + Sp &&
       xx(2) >= lowx && xx(2) < lowx + Sp && xx(3) >= lowx && xx(3) < lowx + Sp)
    {
        int nx = int(floor((xx(0) - lowx)/h));
        int ny = int(floor((xx(1) - lowx)/h));
        int nz = int(floor((xx(2) - lowx)/h));
        int nw = int(floor((xx(3) - lowx)/h));
        int nxy = nx*N+ny;
        int nzw = nz*N+nw;
        //cout<<"nxy = "<<nxy<<" nzw = "<<nzw<<endl;
        if(list_of_box[nxy].size()!=0 && list_of_box[N*N + nzw].size()!=0)
        {
            vector<int> tmp;
            set_intersection(list_of_box[nxy].begin(), list_of_box[nxy].end(),list_of_box[N*N + nzw].begin(), list_of_box[N*N + nzw].end(), back_inserter(tmp));
            if(tmp.size() != 0)
            {
                ret_val = tmp[0];
            }
        }
    }
    return ret_val;
}

void create_boxes(vector<vector<int>>& list_of_box, MatrixXd& Boxes)
{
    Vector4d x_old, x_new;
    x_old<<1.0, 1.0, 1.0, 1.0;
    double eps0 = 0.15;
    random_device rd;
    mt19937 mt(1);
    uniform_real_distribution<double> u(0, 1.0);
    normal_distribution<double> normal(0.0, 1.0);
    Vector4d rnd;
    for(int i = 0; i < Burn_in; i++)
    {
        rnd << normal(mt), normal(mt), normal(mt), normal(mt);
        x_new = x_old + dt*VDP4(x_old) + eps0*sqrt(dt)*rnd;
        x_old = x_new;
    }
    int nx, ny, nz, nw, nxy, nzw;
    int count = 0;
    int flag = 0;
    while(count < N_box)
    {
        flag = 0;
        if(u(mt) < ratio)
        {
            while(flag == 0)
            {
                for(int i = 0; i < 5000; i++)
                {
                    rnd << normal(mt), normal(mt), normal(mt), normal(mt);
                    x_new = x_old + dt*VDP4(x_old) + eps0*sqrt(dt)*rnd;
                    x_old = x_new;
                }
                if(x_new(0) >= lowx && x_new(0) < lowx + Sp
                   && x_new(1) >= lowx && x_new(1) < lowx + Sp
                   && x_new(2) >= lowx && x_new(2) < lowx + Sp
                   && x_new(3) >= lowx && x_new(3) < lowx + Sp)
                {
                    flag = 1;
                }
            }
            nx = int(floor((x_new(0) - lowx)/h));
            ny = int(floor((x_new(1) - lowx)/h));
            nz = int(floor((x_new(2) - lowx)/h));
            nw = int(floor((x_new(3) - lowx)/h));
        }
        else
        {
            nx = int(N*u(mt));
            ny = int(N*u(mt));
            nz = int(N*u(mt));
            nw = int(N*u(mt));
        }
        nxy = nx*N+ny;
        nzw = nz*N+nw;
        vector<int> tmp;
        set_intersection(list_of_box[nxy].begin(), list_of_box[nxy].end(),list_of_box[N*N + nzw].begin(), list_of_box[N*N + nzw].end(), back_inserter(tmp));
        if(tmp.size() == 0)
        {
            list_of_box[nxy].push_back(count);
            sort(list_of_box[nxy].begin(), list_of_box[nxy].end());
            list_of_box[N*N + nzw].push_back(count);
            sort(list_of_box[N*N + nzw].begin(), list_of_box[N*N + nzw].end());
            Boxes(0, count) = lowx + nx*h + h/2.0;
            Boxes(1, count) = lowx + ny*h + h/2.0;
            Boxes(2, count) = lowx + nz*h + h/2.0;
            Boxes(3, count) = lowx + nw*h + h/2.0;
            count++;
        }
    }
    /*
    // for test only
    for(int i = 0; i < 2; i++)
    {
        cout<<"dimension "<<i+1<<" and "<<i+2<<endl;
        for(int j = 0; j < N*N; j++)
        {
            cout<<"index = "<<j<<" : ";
            if(list_of_box[i*N*N + j].size())
            {
                for(int k = 0; k < list_of_box[i*N*N + j].size(); k++)
                {
                    cout<<list_of_box[i*N*N+j][k]<<" ";
                }
            }
            cout<<endl;
        }
    }
    */
}


void MC(VectorXd& Box_count, vector<vector<int>>& list_of_box, const int N_thread, const double eps0)
{
    cout<<"begin parallel session"<<endl;
    VectorXi count(N_thread);
    int sum = 0;
#pragma omp parallel num_threads(N_thread)
    {
        int rank = omp_get_thread_num();
        int size = omp_get_num_threads();
        random_device rd;
        mt19937 mt(rd() + rank);
        normal_distribution<double> norm(0.0, 1.0);
        Vector4d x_old;
        x_old<<5*norm(mt), 5*norm(mt), 5*norm(mt), 5*norm(mt);
        Vector4d x_new;
        int count_f = 0;
        int index;
        for(int i = 0; i < Burn_in; i++)
        {
            Vector4d rnd;
            rnd << norm(mt), norm(mt), norm(mt), norm(mt);
            x_new = x_old + dt*VDP4(x_old) + eps0*sqrt(dt)*rnd;
            x_old = x_new;
        }
        for(int i = 0; i < N_sample; i++)
        {
            Vector4d rnd;
            rnd << norm(mt), norm(mt), norm(mt), norm(mt);
            x_new = x_old + dt*VDP4(x_old) + eps0*sqrt(dt)*rnd;
            index = which_box(list_of_box, x_new);
            if(index!= -1)
            {
//                cout<<"step = "<<i<<" x = "<< x_new.transpose()<<endl;
//                cout<<"index = "<<index<<endl;
                Box_count(rank*N_box + index) += 1.0;
                count_f++;
            }
            x_old = x_new;
        }
        count(rank) = count_f;
    }
    for(int i = 0; i < N_thread; i++)
    {
        sum += count(i);
    }
    cout<<"total effective count = "<<sum<<endl;
    //cout<<Box_count<<endl;
    for(int i = 1; i < N_thread; i++)
    {
        for(int j = 0; j < N_box; j++)
        {
            Box_count(j) += Box_count(i*N_box + j);
        }
    }
}

void create_boxes_slice(vector<vector<int>>& list_of_box, MatrixXd& Boxes)
{

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            int count = i*N + j;
            //cout<<"count = "<<count<<endl;
            int nx = 0;
            int ny = i;
            int nz = j;
            int nw = 0;
            int nxy = nx*N+ny;
            int nzw = nz*N+nw;
            vector<int> tmp;
            set_intersection(list_of_box[nxy].begin(), list_of_box[nxy].end(),list_of_box[N*N + nzw].begin(), list_of_box[N*N + nzw].end(), back_inserter(tmp));
            if(tmp.size() == 0)
            {
                list_of_box[nxy].push_back(count);
                sort(list_of_box[nxy].begin(), list_of_box[nxy].end());
                list_of_box[N*N + nzw].push_back(count);
                sort(list_of_box[N*N + nzw].begin(), list_of_box[N*N + nzw].end());
                Boxes(0, count) = lowx + nx*h + h/2.0;
                Boxes(1, count) = lowx + ny*h + h/2.0;
                Boxes(2, count) = lowx + nz*h + h/2.0;
                Boxes(3, count) = lowx + nw*h + h/2.0;
            }
        }
    }
    /*
     for(int i = 0; i < 2; i++)
     {
     cout<<"dimension "<<i+1<<" and "<<i+2<<endl;
     for(int j = 0; j < N*N; j++)
     {
     cout<<"index = "<<j<<" : ";
     if(list_of_box[i*N*N + j].size())
     {
     for(int k = 0; k < list_of_box[i*N*N + j].size(); k++)
     {
     cout<<list_of_box[i*N*N+j][k]<<" ";
     }
     }
     cout<<endl;
     }
     }
    */
    
}

/*
void training_data2(MatrixXd& X_train2)
//data without y
{
    Vector4d x_old;
    x_old<<1.0, 1.0, 1.0, 1.0;
    Vector4d x_new;
    random_device rd;
    mt19937 mt(rd());
    normal_distribution<double> normal(0.0, 1.0);
    uniform_real_distribution<double> u(0.0,1.0);
    Vector4d rnd;
    for(int i = 0; i < Burn_in; i++)
    {
        rnd << normal(mt), normal(mt), normal(mt), normal(mt);
        x_new = x_old + dt*Ring(x_old(0),x_old(1), x_old(2), x_old(3)) + eps0*sqrt(dt)*rnd;
        x_old = x_new;
    }
    int count = 0;
    while(count < N_train_2)
    {
        if(u(mt) < ratio)
        {
            for(int i = 0; i < 500; i++)
            {
                rnd << normal(mt), normal(mt), normal(mt), normal(mt);
                x_new = x_old + dt*Ring(x_old(0),x_old(1), x_old(2), x_old(3)) + eps0*sqrt(dt)*rnd;
                x_old = x_new;
            }
            if(x_new(0) >= lowx && x_new(0) < lowx + Sp
               && x_new(1) >= lowx && x_new(1) < lowx + Sp
               && x_new(2) >= lowx && x_new(2) < lowx + Sp
               && x_new(3) >= lowx && x_new(3) < lowx + Sp)
            {
                X_train2.col(count) = x_new;
                count++;
            }
        }
        else
        {
            X_train2.col(count) << lowx + Sp*u(mt), lowx + Sp*u(mt), lowx + Sp*u(mt), lowx + Sp*u(mt);
            count ++;
        }
        //cout<<"count = "<<count<<endl;
        //cout<<X_train2.col(count-1)<<endl;
    }
}

void slice_X_train(MatrixXd& X_test)
{
    for(int i = 0; i < N_X; i++)
    {
        for(int j = 0; j < N_X; j++)
        {
            double xx = lowx + i*h + h/2;
            double yy = lowx + j*h + h/2;
            X_test(0, i*N_X + j) = xx;
            X_test(1, i*N_X + j) = yy;
            X_test(2, i*N_X + j) = 0;
            X_test(3, i*N_X + j) = 0;
        }
    }
}

*/

int main ( )
{
    struct timeval t1, t2;
    random_device rd;
    mt19937 mt(rd());
    normal_distribution<double> norm(0.0, 1.0);
    uniform_real_distribution<double> u(0, 1.0);
    int N_thread = 10;
    
    gettimeofday(&t1,NULL);
    vector<vector<int>> list_of_box(2*N*N);//indices of boxes in each projection
    //first N_X * N_X is for x,y, second is for z, w
    // (x*N + y) and (z*N + w)
    //coordinates
    for(int i = 0; i < 2*N*N; i++)
    {
        list_of_box[i].reserve(10);
    }
    MatrixXd Boxes(4, N_box);//coordinate of boxes
    //Boxes <-> list_of_box
    VectorXd Box_count(N_thread*N_box);
    //MatrixXd X_train2(4, N_train_2);
    MatrixXd Y_train(9, N_box);
    create_boxes(list_of_box, Boxes);
    //cout<<Boxes<<endl;
    cout<<"boxes generated"<<endl;
    
    double eps_list[] = {0.05, 0.06, 0.075, 0.085, 0.1, 0.11, 0.125, 0.135, 0.15};
    for(int i = 0; i < 9; i++)
    {
        cout<<"noise = "<<eps_list[i]<<endl;
        Box_count.fill(0);
        MC(Box_count, list_of_box, N_thread, eps_list[i]);
        Y_train.row(i) = Box_count.head(N_box)/(N_sample*pow(h,4));
    }
    cout<<"Y train generated"<<endl;
    
    ofstream myfile1, myfile2;
    myfile1.open("VDP4_X_train.txt");
    myfile2.open("VDP4_Y_train.txt");
    for(int i = 0 ; i < N_box; i++)
    {
        myfile1<<Boxes(0, i)<<" "<<Boxes(1, i)<<" "<<Boxes(2, i)<<" "<<Boxes(3, i)<<endl;
    }
    for(int i = 0; i < N_box; i++)
    {
        for(int j = 0; j < 9; j++)
        {
            myfile2<<Y_train(j, i)<<" ";
        }
        myfile2<<endl;
    }
    
    

    
    gettimeofday(&t2, NULL);
    
    double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u +
                    t2.tv_usec - t1.tv_usec) / 1.e6;
    
    cout<<"total CPU time = " <<delta<<endl;

    
    return 0 ;
    
}
