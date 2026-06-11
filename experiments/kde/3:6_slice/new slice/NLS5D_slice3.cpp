#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>
#include <random>
#include <vector>
#include <algorithm>
#include <omp.h>

using namespace std;

double gamma_val = 0.1;
double T1        = 5.0;
double T3        = 5.0;
double dt_val    = 0.001;
double sqrt_dt_val;
const int Burn_in = 2000000;

double I_slice = 2.0;
double I_tol   = 0.5;

const double Th_lo = -M_PI, Th_hi = M_PI;
const int G_bin = 50;
const double d_bin = (Th_hi - Th_lo) / G_bin;
const int G_kde = 30;
const double d_kde = (Th_hi - Th_lo) / G_kde;

long long N_sample = 100000000;
int       N_thread = 8;

struct State6 {
    double v[6];
    double& operator()(int i) { return v[i]; }
    const double& operator()(int i) const { return v[i]; }
};

inline double wrap(double x) {
    return x - 2.0 * M_PI * floor((x + M_PI) / (2.0 * M_PI));
}


void step_EM(State6& X, mt19937& rng, normal_distribution<double>& dist)
{
    double I1 = X(0), I2 = X(1), I3 = X(2);
    double p1 = X(3), p2 = X(4), p3 = X(5);

    double I_prev[3], I_next[3];
    I_prev[0] = 0.0;   I_next[0] = I2;
    I_prev[1] = I1;    I_next[1] = I3;
    I_prev[2] = I2;    I_next[2] = 0.0;

    double d_phi_prev[3], d_phi_next[3];
    d_phi_prev[0] = 2.0*(p1 - 0.0);
    d_phi_next[0] = 2.0*(p1 - p2);
    d_phi_prev[1] = 2.0*(p2 - p1);
    d_phi_next[1] = 2.0*(p2 - p3);
    d_phi_prev[2] = 2.0*(p3 - p2);
    d_phi_next[2] = 2.0*(p3 - 0.0);

    double Ivec[3] = {I1, I2, I3};
    double M = I1 + I2 + I3;

    double dI[3], dphi[3];
    for(int j = 0; j < 3; j++) {
        dI[j] = 4.0 * Ivec[j] * (I_prev[j]*sin(d_phi_prev[j])
                                 + I_next[j]*sin(d_phi_next[j]));
        dphi[j] = 2.0*M - Ivec[j]
                  + 2.0*I_prev[j]*cos(d_phi_prev[j])
                  + 2.0*I_next[j]*cos(d_phi_next[j]);
    }

    dI[0] += 2.0*gamma_val*(2.0*T1 - (2.0*M*I1 - I1*I1 + 2.0*I_next[0]*I1*cos(d_phi_next[0])));

    // wrong: dphi[0] -= gamma_val*(2.0*I_next[0]*sin(d_phi_next[0]));
    // correct:
    dphi[0] += gamma_val*(2.0*I_next[0]*sin(d_phi_next[0]));

    dI[2] += 2.0*gamma_val*(2.0*T3 - (2.0*M*I3 - I3*I3 + 2.0*I_prev[2]*I3*cos(d_phi_prev[2])));
    
    // wrong: dphi[2] -= gamma_val*(2.0*I_prev[2]*sin(d_phi_prev[2]));
    // correct
    dphi[2] += gamma_val*(2.0*I_prev[2]*sin(d_phi_prev[2]));

    double sigma_I1  = 2.0*sqrt(max(2.0*gamma_val*T1*I1, 0.0));
    double sigma_I3  = 2.0*sqrt(max(2.0*gamma_val*T3*I3, 0.0));
    double sigma_p1  = sqrt(max(2.0*gamma_val*T1/max(I1,1e-14), 0.0));
    double sigma_p3  = sqrt(max(2.0*gamma_val*T3/max(I3,1e-14), 0.0));

    double dW[4];
    for(int j = 0; j < 4; j++) dW[j] = sqrt_dt_val * dist(rng);

    X(0) = I1  + dI[0]*dt_val + sigma_I1*dW[0];
    X(1) = I2  + dI[1]*dt_val;
    X(2) = I3  + dI[2]*dt_val + sigma_I3*dW[1];
    X(3) = p1  + dphi[0]*dt_val + sigma_p1*dW[2];
    X(4) = p2  + dphi[1]*dt_val;
    X(5) = p3  + dphi[2]*dt_val + sigma_p3*dW[3];

    X(0) = max(X(0), 1e-14);
    X(1) = max(X(1), 1e-14);
    X(2) = max(X(2), 1e-14);
    X(3) = wrap(X(3));
    X(4) = wrap(X(4));
    X(5) = wrap(X(5));
}

inline bool in_slice(double I1, double I2, double I3) {
    return fabs(I1 - I_slice) < I_tol
        && fabs(I2 - I_slice) < I_tol
        && fabs(I3 - I_slice) < I_tol;
}

inline double normpdf_periodic(double x, double mu, double h) {
    double d = x - mu;
    if(d >  M_PI) d -= 2.0 * M_PI;
    if(d < -M_PI) d += 2.0 * M_PI;
    double z = d / h;
    return exp(-0.5 * z * z) / (h * sqrt(2.0 * M_PI));
}

double bessel_I0(double x) {
    double sum = 1.0, term = 1.0;
    for(int k = 1; k < 200; k++) {
        term *= (x / (2.0 * k)) * (x / (2.0 * k));
        sum += term;
        if(term < 1e-16 * sum) break;
    }
    return sum;
}

// Theoretical density: exp(-a*(cos(th1)+cos(th3)))
// invariant measure exp(-H_code/(2T)): a = 2I²/(2T) = I²/T
void compute_theoretical_density(double I_val, double T, int G, double dg,
                                  vector<double>& out) {
    double a = I_val * I_val / T;  // I²/T 
    double I0a = bessel_I0(a);
    double Z = (2.0*M_PI*I0a) * (2.0*M_PI*I0a);

    cout << "  Theoretical: a = I^2/T = " << a
         << ",  I_0(a) = " << I0a << ",  Z = " << Z << endl;

    double rho_min = exp(-2.0*a) / Z;
    double rho_max = exp( 2.0*a) / Z;
    cout << "  Density range: [" << rho_min << ", " << rho_max << "]"
         << "  Dynamic range: " << rho_max/rho_min << "x" << endl;

    out.resize(G * G);
    for(int j1 = 0; j1 < G; j1++) {
        double th1 = Th_lo + (j1 + 0.5) * dg;
        for(int j3 = 0; j3 < G; j3++) {
            double th3 = Th_lo + (j3 + 0.5) * dg;
            out[j1*G+j3] = exp(-a*(cos(th1)+cos(th3))) / Z;
        }
    }
}

int main(int argc, char* argv[]) {
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL);

    if(argc > 1) N_thread  = atoi(argv[1]);
    if(argc > 2) gamma_val = atof(argv[2]);
    if(argc > 3) T1        = atof(argv[3]);
    if(argc > 4) T3        = atof(argv[4]);
    if(argc > 5) N_sample  = atoll(argv[5]);
    if(argc > 6) I_slice   = atof(argv[6]);
    if(argc > 7) I_tol     = atof(argv[7]);
    if(argc > 8) dt_val    = atof(argv[8]);

    sqrt_dt_val = sqrt(dt_val);
    bool equilibrium = (fabs(T1 - T3) < 1e-10);

    cout << "=== NLS 5D: 2D Slice===" << endl;
    cout << "Slice: I1=I2=I3=" << I_slice << " +/- " << I_tol << endl;
    cout << "gamma=" << gamma_val << " T1=" << T1 << " T3=" << T3;
    if(equilibrium) cout << "  [EQUILIBRIUM]";
    cout << endl;
    cout << "dt=" << dt_val << " Burn_in=" << Burn_in << endl;
    cout << "N_sample/thread=" << N_sample << " N_thread=" << N_thread << endl;
    if(equilibrium) {
        double a = I_slice*I_slice/T1;
        cout << "Theoretical: a = I^2/T = " << a << endl;
    }
    cout << endl;

    cout << "Collecting slice samples..." << endl;
    vector<double> all_th1, all_th3;
    long long total_steps = 0, total_in_slice = 0;

    #pragma omp parallel num_threads(N_thread)
    {
        int rank = omp_get_thread_num();
        mt19937 rng(random_device{}() + rank*137);
        normal_distribution<double> dist(0.0, 1.0);

        State6 X;
        X(0) = 1.0+0.1*rank; X(1) = 1.0; X(2) = 0.1;
        X(3) = 0.1*rank; X(4) = 0.0; X(5) = -0.1*rank;

        for(int i = 0; i < Burn_in; i++)
            step_EM(X, rng, dist);

        vector<double> loc_th1, loc_th3;
        long long loc_steps = 0, loc_sl = 0;

        for(long long step = 0; step < N_sample; step++) {
            step_EM(X, rng, dist);
            loc_steps++;
            if(in_slice(X(0), X(1), X(2))) {
                loc_th1.push_back(wrap(2.0*(X(3)-X(4))));
                loc_th3.push_back(wrap(2.0*(X(5)-X(4))));
                loc_sl++;
            }
            if(rank==0 && step%50000000LL==0 && step>0)
                cout << "  thread 0: " << step/1000000 << "M/" << N_sample/1000000 << "M" << endl;
        }

        #pragma omp critical
        {
            all_th1.insert(all_th1.end(), loc_th1.begin(), loc_th1.end());
            all_th3.insert(all_th3.end(), loc_th3.begin(), loc_th3.end());
            total_steps += loc_steps;
            total_in_slice += loc_sl;
        }
    }

    long long N_sl = all_th1.size();
    cout << "Total steps: " << total_steps << endl;
    cout << "In slice: " << N_sl << " (" << 100.0*N_sl/max(total_steps,1LL) << "%)" << endl;
    cout << endl;

    if(N_sl < 100) { cout << "Too few samples.\n"; return 1; }

    //  <cos(theta)>
    double mean_cos1 = 0, mean_cos3 = 0;
    for(long long s = 0; s < N_sl; s++) {
        mean_cos1 += cos(all_th1[s]);
        mean_cos3 += cos(all_th3[s]);
    }
    mean_cos1 /= N_sl; mean_cos3 /= N_sl;
    cout << "Measured <cos(theta1)> = " << mean_cos1 << endl;
    cout << "Measured <cos(theta3)> = " << mean_cos3 << endl;
    if(equilibrium) {
        double a = I_slice*I_slice/T1;
        // I_1(a)/I_0(a) via series
        double I0a = bessel_I0(a);
        double I1a = 0;
        for(int k = 0; k < 200; k++) {
            double t = 1.0;
            for(int j = 1; j <= k; j++) t *= (a/2.0)/(double)j;
            for(int j = 1; j <= k+1; j++) t *= (a/2.0)/(double)j;
            I1a += t;
        }
        cout << "Theory  <cos(theta)>  = " << (-I1a/I0a) << " (for a=" << a << ")" << endl;
    }
    cout << endl;

    // Binning
    cout << "=== Binning " << G_bin << "x" << G_bin << " ===" << endl;
    vector<double> hist(G_bin*G_bin, 0.0);
    for(long long s = 0; s < N_sl; s++) {
        int j1 = max(0, min((int)floor((all_th1[s]-Th_lo)/d_bin), G_bin-1));
        int j3 = max(0, min((int)floor((all_th3[s]-Th_lo)/d_bin), G_bin-1));
        hist[j1*G_bin+j3] += 1.0;
    }
    double ca = d_bin*d_bin;
    vector<double> dens_bin(G_bin*G_bin);
    for(int i = 0; i < G_bin*G_bin; i++) dens_bin[i] = hist[i]/((double)N_sl*ca);

    // KDE
    cout << "=== KDE " << G_kde << "x" << G_kde << " ===" << endl;
    double h_vals[2] = {0.05, 0.10};
    vector<vector<double>> kde_res(2);
    for(int ih=0; ih<2; ih++) {
        double h = h_vals[ih];
        cout << "  h=" << h << "..." << endl;
        kde_res[ih].resize(G_kde*G_kde, 0.0);
        #pragma omp parallel for num_threads(N_thread) schedule(dynamic,4)
        for(int j1=0; j1<G_kde; j1++) {
            double th1c = Th_lo+(j1+0.5)*d_kde;
            for(int j3=0; j3<G_kde; j3++) {
                double th3c = Th_lo+(j3+0.5)*d_kde;
                double sum=0;
                for(long long s=0; s<N_sl; s++)
                    sum += normpdf_periodic(th1c,all_th1[s],h)
                         * normpdf_periodic(th3c,all_th3[s],h);
                kde_res[ih][j1*G_kde+j3] = sum/(double)N_sl;
            }
        }
    }

    // Theory comparison
    if(equilibrium) {
        cout << endl << "=== Theory comparison ===" << endl;
        vector<double> theory30;
        compute_theoretical_density(I_slice, T1, G_kde, d_kde, theory30);

        cout << endl << "=== Error metrics ===" << endl;
        for(int ih=0; ih<2; ih++) {
            double l2=0, linf=0;
            for(int i=0;i<G_kde*G_kde;i++){
                double diff=fabs(kde_res[ih][i]-theory30[i]);
                l2+=diff*diff; if(diff>linf) linf=diff;
            }
            l2=sqrt(l2/(G_kde*G_kde));
            printf("  KDE h=%.2f: RMSE=%.6f  MaxErr=%.6f\n", h_vals[ih], l2, linf);
        }

        vector<double> theory50;
        compute_theoretical_density(I_slice, T1, G_bin, d_bin, theory50);
        {
            double l2=0, linf=0;
            for(int i=0;i<G_bin*G_bin;i++){
                double diff=fabs(dens_bin[i]-theory50[i]);
                l2+=diff*diff; if(diff>linf) linf=diff;
            }
            l2=sqrt(l2/(G_bin*G_bin));
            printf("  Binning 50x50: RMSE=%.6f  MaxErr=%.6f\n", l2, linf);
        }

        // Cross-section
        cout << endl << "  Cross-section theta3~0:" << endl;
        int j3m = G_kde/2;
        printf("  %8s %10s %10s %10s\n", "theta1","Theory","KDE_h0.05","KDE_h0.10");
        for(int j1=0; j1<G_kde; j1++) {
            double th1 = Th_lo+(j1+0.5)*d_kde;
            printf("  %8.3f %10.6f %10.6f %10.6f\n",
                   th1, theory30[j1*G_kde+j3m], kde_res[0][j1*G_kde+j3m], kde_res[1][j1*G_kde+j3m]);
        }
    }

    // === Write 2D density data to CSV files ===
    // Binning output
    {
        ofstream fout("density_bin.csv");
        fout << "theta1,theta3,density\n";
        for(int j1=0; j1<G_bin; j1++){
            double th1 = Th_lo + (j1+0.5)*d_bin;
            for(int j3=0; j3<G_bin; j3++){
                double th3 = Th_lo + (j3+0.5)*d_bin;
                fout << th1 << "," << th3 << "," << dens_bin[j1*G_bin+j3] << "\n";
            }
        }
        fout.close();
        cout << "Wrote density_bin.csv (" << G_bin << "x" << G_bin << ")" << endl;
    }

    // KDE outputs
    for(int ih=0; ih<2; ih++){
        char fname[64];
        snprintf(fname, sizeof(fname), "density_kde_h%.2f.csv", h_vals[ih]);
        ofstream fout(fname);
        fout << "theta1,theta3,density\n";
        for(int j1=0; j1<G_kde; j1++){
            double th1 = Th_lo + (j1+0.5)*d_kde;
            for(int j3=0; j3<G_kde; j3++){
                double th3 = Th_lo + (j3+0.5)*d_kde;
                fout << th1 << "," << th3 << "," << kde_res[ih][j1*G_kde+j3] << "\n";
            }
        }
        fout.close();
        cout << "Wrote " << fname << " (" << G_kde << "x" << G_kde << ")" << endl;
    }

    // Theory output (if equilibrium)
    if(equilibrium){
        vector<double> theory_out;
        compute_theoretical_density(I_slice, T1, G_kde, d_kde, theory_out);
        ofstream fout("density_theory.csv");
        fout << "theta1,theta3,density\n";
        for(int j1=0; j1<G_kde; j1++){
            double th1 = Th_lo + (j1+0.5)*d_kde;
            for(int j3=0; j3<G_kde; j3++){
                double th3 = Th_lo + (j3+0.5)*d_kde;
                fout << th1 << "," << th3 << "," << theory_out[j1*G_kde+j3] << "\n";
            }
        }
        fout.close();
        cout << "Wrote density_theory.csv (" << G_kde << "x" << G_kde << ")" << endl;
    }

    gettimeofday(&t_end, NULL);
    double sec = ((t_end.tv_sec-t_start.tv_sec)*1e6+t_end.tv_usec-t_start.tv_usec)/1e6;
    cout << endl << "Wall time = " << sec << "s" << endl;
    return 0;
}