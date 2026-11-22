// Four Hamiltonian integrators, compared pathwise (BC1, n = 25, T = (10,5), same initial state and Brownian path,
// fine grid hf = 5e-4/16), all in the FIXED frame (needed so that one noise realisation means the same thing for all):
//   1 RK4_C   : classic RK4 in Cartesian c = x + i y
//   2 RK4_A   : classic RK4 in action-angle (I, phi)
//   3 EPS_C   : explicit symplectic method in an extended phase space (EPS) for nonseparable H, 2nd order, in Cartesian
//               (method of M. Tao, Phys. Rev. E 94, 043303 (2016), arXiv:1609.02212)
//   4 EPS_A   : the same in action-angle
// plus NEW (Cartesian implicit midpoint = new code), OLD (old SIMD step) and REF (NEW at hf).
// Bath: by default every scheme uses the new code's Cartesian Euler-Maruyama bath step after its Hamiltonian step
// (order alternating).  Variants *_AB use the old action-angle bath step instead (for the action-angle schemes).
// EPS method: extended phase space (q,p,x,y), Hbar = H(q,y) + H(x,p) + omega/2 (|q-x|^2 + |p-y|^2),
//   step = A(h/2) B(h/2) C(h) B(h/2) A(h/2).  Canonical pairs: Cartesian (q,p) = (Im c, Re c) with H/2,
//   action-angle (q,p) = (phi, I) with H.  The bath step is applied to the first copy and the same increment to the
//   second copy (their difference is not changed by the bath).
// Usage: ./sch4 n T1 Tn tmax seed omega_C omega_A out.csv
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>
using cd = std::complex<double>;
const double PI = 3.14159265358979323846;
const cd I1(0, 1);
int n; double T1, Tn, g = 0.1; bool g_corot = false;
long old_small_dt = 0, old_steps = 0;

float fast_sin(float x) { const float B = 1.27323954f, C = -0.40528473f, P = 0.225f; float y = B * x + C * x * std::fabs(x); return P * (y * std::fabs(y) - y) + y; }
template <class R> R wrapr(R x) { return x - std::round(x / (R)(2 * PI)) * (R)(2 * PI); }

// ---------------- action-angle explicit Euler (copied logic of NLS_flux_SIMD_fixed.original.cpp, case 0 = BC1)
template <class R, bool FAST>
void aa_step(std::vector<R>& I, std::vector<R>& ph, double h, const double* dW, bool count) {
    std::vector<R> s(n, 0), c(n, 0), dI(n), dp(n);
    R M = 0; for (auto v : I) M += v;
    for (int j = 1; j < n; ++j) { R d = (R)2 * (ph[j] - ph[j - 1]);
        if (FAST) { s[j] = fast_sin(wrapr<float>(d)); c[j] = fast_sin(wrapr<float>(d + 1.5707963268f)); } else { s[j] = std::sin(d); c[j] = std::cos(d); } }
    for (int j = 0; j < n; ++j) {
        R tl = j > 0 ? I[j - 1] * s[j] : 0, tr = j < n - 1 ? -I[j + 1] * s[j + 1] : 0;
        R pl = j > 0 ? 2 * I[j - 1] * c[j] : 0, pr = j < n - 1 ? 2 * I[j + 1] * c[j + 1] : 0;
        dI[j] = 4 * I[j] * (tl + tr); dp[j] = (2 * M - I[j]) + pl + pr; }
    // boundary (BC1): exactly as in the old code
    int site[2] = {0, n - 1}; double T[2] = {T1, Tn};
    R nb_I[2], nb_s[2], nb_c[2];
    nb_I[0] = I[1]; nb_s[0] = -s[1]; nb_c[0] = c[1];
    nb_I[1] = I[n - 2]; nb_s[1] = s[n - 1]; nb_c[1] = c[n - 1];
    for (int b = 0; b < 2; ++b) { int r = site[b];
        dI[r] += (R)(2 * g) * ((R)(2 * T[b]) - ((R)2 * M * I[r] - I[r] * I[r] + (R)2 * nb_I[b] * I[r] * nb_c[b]));
        dp[r] += (R)g * ((R)2 * nb_I[b] * nb_s[b]); }
    if (count) { R mx = 1; for (int j = 0; j < n; ++j) mx = std::max(mx, std::max(std::fabs(dI[j]), std::fabs(dp[j])));
        ++old_steps; if (1.0 / (double)mx < h) ++old_small_dt; }
    R nI[2], np[2];
    for (int b = 0; b < 2; ++b) { int r = site[b]; R Ir = std::max(I[r], (R)1e-10);
        R cph = std::cos(ph[r]), sph = std::sin(ph[r]);
        R dBI = cph * (R)dW[2 * b] + sph * (R)dW[2 * b + 1], dBp = -sph * (R)dW[2 * b] + cph * (R)dW[2 * b + 1];
        nI[b] = (R)2 * std::sqrt((R)(2 * g * T[b]) * Ir) * dBI;
        np[b] = std::sqrt((R)(2 * g * T[b]) / Ir) * dBp; }
    for (int j = 0; j < n; ++j) { I[j] += dI[j] * (R)h; ph[j] += dp[j] * (R)h; }
    for (int b = 0; b < 2; ++b) { int r = site[b]; I[r] += nI[b]; ph[r] += np[b]; }
    for (int j = 0; j < n; ++j) { I[j] = std::max(I[j], (R)1e-10); ph[j] = wrapr<R>(ph[j]); }
}
// ---------------- Cartesian
void force(const std::vector<cd>& c, std::vector<cd>& F) {
    double M = 0; for (auto z : c) M += std::norm(z);
    for (int j = 0; j < n; ++j) { cd N = 0; if (j > 0) N += c[j-1] * c[j-1]; if (j < n - 1) N += c[j+1] * c[j+1];
        F[j] = ((g_corot ? 0.0 : 2 * M) - std::norm(c[j])) * c[j] + 2.0 * N * std::conj(c[j]); }
}
void ce_step(std::vector<cd>& c, double h, const double* dW) {
    std::vector<cd> F(n); force(c, F); std::vector<cd> nc(n);
    for (int j = 0; j < n; ++j) nc[j] = c[j] + h * I1 * F[j];
    int site[2] = {0, n - 1}; double T[2] = {T1, Tn};
    for (int b = 0; b < 2; ++b) { int r = site[b]; nc[r] += -g * F[r] * h + std::sqrt(2 * g * T[b]) * cd(dW[2*b], dW[2*b+1]); }
    c = nc;
}
struct Mid { std::vector<cd> g0, gs, m, F; Mid() : g0(n), gs(n), m(n), F(n) {} };
void midpoint(std::vector<cd>& c, double h, Mid& w) {
    w.g0 = c; force(c, w.F); for (int j = 0; j < n; ++j) w.gs[j] = c[j] + h * I1 * w.F[j];
    for (int it = 0; it < 30; ++it) { for (int j = 0; j < n; ++j) w.m[j] = 0.5 * (w.g0[j] + w.gs[j]);
        force(w.m, w.F); double ch = 0;
        for (int j = 0; j < n; ++j) { cd nw = w.g0[j] + h * I1 * w.F[j]; ch = std::max(ch, std::abs(nw - w.gs[j])); w.gs[j] = nw; }
        if (ch < 1e-14) break; }
    c = w.gs;
}
void bath(std::vector<cd>& c, int r, double T, double h, double w1, double w2, Mid& w) {
    force(c, w.F); c[r] += -g * w.F[r] * h + std::sqrt(2 * g * T) * cd(w1, w2);
}
void cm_step(std::vector<cd>& c, double h, const double* dW, long k, Mid& w) {
    midpoint(c, h, w);
    if (k % 2 == 0) { bath(c, 0, T1, h, dW[0], dW[1], w); bath(c, n - 1, Tn, h, dW[2], dW[3], w); }
    else            { bath(c, n - 1, Tn, h, dW[2], dW[3], w); bath(c, 0, T1, h, dW[0], dW[1], w); }
}

// ---- action-angle Hamiltonian drift (fixed frame, includes 2M)
void aa_rhs(const std::vector<double>& I, const std::vector<double>& p, std::vector<double>& dI, std::vector<double>& dp) {
    double M = 0; for (auto v : I) M += v;
    for (int j = 0; j < n; ++j) { double tI = 0, tp = 0;
        if (j > 0)   { double d = 2*(p[j]-p[j-1]); tI += I[j-1]*std::sin(d); tp += 2*I[j-1]*std::cos(d); }
        if (j < n-1) { double d = 2*(p[j+1]-p[j]); tI -= I[j+1]*std::sin(d); tp += 2*I[j+1]*std::cos(d); }
        dI[j] = 4*I[j]*tI; dp[j] = (2*M - I[j]) + tp; }
}
void aa_rk4(std::vector<double>& I, std::vector<double>& p, double h) {
    std::vector<double> a1(n),b1(n),a2(n),b2(n),a3(n),b3(n),a4(n),b4(n),It(n),pt(n);
    aa_rhs(I,p,a1,b1); for(int j=0;j<n;++j){It[j]=I[j]+0.5*h*a1[j]; pt[j]=p[j]+0.5*h*b1[j];}
    aa_rhs(It,pt,a2,b2); for(int j=0;j<n;++j){It[j]=I[j]+0.5*h*a2[j]; pt[j]=p[j]+0.5*h*b2[j];}
    aa_rhs(It,pt,a3,b3); for(int j=0;j<n;++j){It[j]=I[j]+h*a3[j]; pt[j]=p[j]+h*b3[j];}
    aa_rhs(It,pt,a4,b4);
    for(int j=0;j<n;++j){ I[j]+=h/6*(a1[j]+2*a2[j]+2*a3[j]+a4[j]); p[j]+=h/6*(b1[j]+2*b2[j]+2*b3[j]+b4[j]); }
}
void aa_euler(std::vector<double>& I, std::vector<double>& p, double h) {
    std::vector<double> a(n), b(n); aa_rhs(I,p,a,b); for(int j=0;j<n;++j){ I[j]+=h*a[j]; p[j]+=h*b[j]; }
}
// old action-angle BC1 bath step at site r (drift of the old code + rotated noise), Euler-Maruyama, clamp
void aa_bath(std::vector<double>& I, std::vector<double>& p, int r, double T, double h, double w1, double w2) {
    double M = 0; for (auto v : I) M += v;
    int nb = (r == 0) ? 1 : n - 2; double d = 2*(p[r]-p[nb]);   // theta seen from site r
    double dI = 2*g*(2*T - (2*M*I[r] - I[r]*I[r] + 2*I[nb]*I[r]*std::cos(d)));
    double dp = g*(2*I[nb]*std::sin(d));
    double Ir = std::max(I[r], 1e-10), cph = std::cos(p[r]), sph = std::sin(p[r]);
    double dBI = cph*w1 + sph*w2, dBp = -sph*w1 + cph*w2;
    I[r] += dI*h + 2*std::sqrt(2*g*T*Ir)*dBI; p[r] += dp*h + std::sqrt(2*g*T/Ir)*dBp;
    I[r] = std::max(I[r], 1e-10);
}
void to_c(const std::vector<double>& I, const std::vector<double>& p, std::vector<cd>& c) { for (int j=0;j<n;++j) c[j] = std::polar(std::sqrt(I[j]), p[j]); }
void to_a(const std::vector<cd>& c, std::vector<double>& I, std::vector<double>& p) { for (int j=0;j<n;++j){ I[j]=std::norm(c[j]); p[j]=std::arg(c[j]); } }
void cart_rk4(std::vector<cd>& c, double h) {
    std::vector<cd> k1(n),k2(n),k3(n),k4(n),t(n),F(n);
    force(c,F); for(int j=0;j<n;++j){k1[j]=I1*F[j]; t[j]=c[j]+0.5*h*k1[j];}
    force(t,F); for(int j=0;j<n;++j){k2[j]=I1*F[j]; t[j]=c[j]+0.5*h*k2[j];}
    force(t,F); for(int j=0;j<n;++j){k3[j]=I1*F[j]; t[j]=c[j]+h*k3[j];}
    force(t,F); for(int j=0;j<n;++j){k4[j]=I1*F[j]; c[j]+=h/6.0*(k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j]);}
}

// ---- EPS, Cartesian.  First copy c = X + iY (q = Y, p = X), second copy d = U + iV (q = V, p = U).
//   A (H(q of c, p of d)): z = U + iY;  X -= d*Im F(z);  V += d*Re F(z)
//   B (H(q of d, p of c)): z = X + iV;  Y += d*Re F(z);  U -= d*Im F(z)
//   C: rotate (q_c - q_d, p_c - p_d) by angle 2 omega d, keep sums.
void epsC_A(std::vector<cd>& c, std::vector<cd>& d, double dl) { std::vector<cd> z(n), F(n);
    for (int j=0;j<n;++j) z[j] = cd(d[j].real(), c[j].imag()); force(z, F);
    for (int j=0;j<n;++j) { c[j] = cd(c[j].real() - dl*F[j].imag(), c[j].imag()); d[j] = cd(d[j].real(), d[j].imag() + dl*F[j].real()); } }
void epsC_B(std::vector<cd>& c, std::vector<cd>& d, double dl) { std::vector<cd> z(n), F(n);
    for (int j=0;j<n;++j) z[j] = cd(c[j].real(), d[j].imag()); force(z, F);
    for (int j=0;j<n;++j) { c[j] = cd(c[j].real(), c[j].imag() + dl*F[j].real()); d[j] = cd(d[j].real() - dl*F[j].imag(), d[j].imag()); } }
void rotC(double& q1, double& p1, double& q2, double& p2, double ang) {
    double sq = q1+q2, sp = p1+p2, dq = q1-q2, dp = p1-p2, cs = std::cos(ang), sn = std::sin(ang);
    double dq2 = cs*dq + sn*dp, dp2 = -sn*dq + cs*dp;
    q1 = 0.5*(sq+dq2); q2 = 0.5*(sq-dq2); p1 = 0.5*(sp+dp2); p2 = 0.5*(sp-dp2); }
void epsC_C(std::vector<cd>& c, std::vector<cd>& d, double dl, double om) {
    for (int j=0;j<n;++j) { double q1=c[j].imag(), p1=c[j].real(), q2=d[j].imag(), p2=d[j].real(); rotC(q1,p1,q2,p2,2*om*dl); c[j]=cd(p1,q1); d[j]=cd(p2,q2); } }
void eps_cart(std::vector<cd>& c, std::vector<cd>& d, double h, double om) {
    epsC_A(c,d,h/2); epsC_B(c,d,h/2); epsC_C(c,d,h,om); epsC_B(c,d,h/2); epsC_A(c,d,h/2); }
// ---- action-angle gradients (fixed frame, H includes M^2): HI = dH/dI, Hp = dH/dphi
void aa_grad(const std::vector<double>& I, const std::vector<double>& p, std::vector<double>& HI, std::vector<double>& Hp) {
    double M = 0; for (auto v : I) M += v;
    for (int j=0;j<n;++j) { double a = 2*M - I[j], b = 0;
        if (j>0)   { double t = 2*(p[j]-p[j-1]); a += 2*I[j-1]*std::cos(t); b += -4*I[j]*I[j-1]*std::sin(t); }
        if (j<n-1) { double t = 2*(p[j+1]-p[j]); a += 2*I[j+1]*std::cos(t); b +=  4*I[j]*I[j+1]*std::sin(t); }
        HI[j] = a; Hp[j] = b; } }
// EPS in action-angle: copy 1 (phi1, I1), copy 2 (phi2, I2); q = phi, p = I.
//   A: H(q1, p2): grads at (I2, phi1): I1 -= d*Hp;  phi2 += d*HI
//   B: H(q2, p1): grads at (I1, phi2): phi1 += d*HI; I2 -= d*Hp
void eps_aa(std::vector<double>& I1v, std::vector<double>& P1, std::vector<double>& I2v, std::vector<double>& P2, double h, double om) {
    std::vector<double> HI(n), Hp(n);
    auto A = [&](double dl){ aa_grad(I2v, P1, HI, Hp); for (int j=0;j<n;++j){ I1v[j] -= dl*Hp[j]; P2[j] += dl*HI[j]; } };
    auto B = [&](double dl){ aa_grad(I1v, P2, HI, Hp); for (int j=0;j<n;++j){ P1[j] += dl*HI[j]; I2v[j] -= dl*Hp[j]; } };
    A(h/2); B(h/2);
    for (int j=0;j<n;++j) { double q1=P1[j], p1=I1v[j], q2=P2[j], p2=I2v[j]; rotC(q1,p1,q2,p2,2*om*h); P1[j]=q1; I1v[j]=p1; P2[j]=q2; I2v[j]=p2; }
    B(h/2); A(h/2); }
struct S { int kind; bool aabath; double h; int sub; std::vector<double> I, p, I2, p2; std::vector<cd> c, d; long k = 0; double acc[4] = {0,0,0,0}; };
double omC = 100, omA = 100;
void bath_pair_cart(std::vector<cd>& c, double h, const double* dW, long k, Mid& w) {
    if (k % 2 == 0) { bath(c, 0, T1, h, dW[0], dW[1], w); bath(c, n-1, Tn, h, dW[2], dW[3], w); }
    else            { bath(c, n-1, Tn, h, dW[2], dW[3], w); bath(c, 0, T1, h, dW[0], dW[1], w); } }
void bath_pair_aa(std::vector<double>& I, std::vector<double>& p, double h, const double* dW, long k) {
    if (k % 2 == 0) { aa_bath(I,p,0,T1,h,dW[0],dW[1]); aa_bath(I,p,n-1,Tn,h,dW[2],dW[3]); }
    else            { aa_bath(I,p,n-1,Tn,h,dW[2],dW[3]); aa_bath(I,p,0,T1,h,dW[0],dW[1]); } }
void step(S& s, Mid& w) {
    const double* dW = s.acc;
    if (s.kind == 1) { cart_rk4(s.c, s.h); bath_pair_cart(s.c, s.h, dW, s.k, w); }
    else if (s.kind == 3) { eps_cart(s.c, s.d, s.h, omC); std::vector<cd> before = s.c; bath_pair_cart(s.c, s.h, dW, s.k, w);
        for (int j=0;j<n;++j) s.d[j] += s.c[j] - before[j]; }
    else if (s.kind == 6) cm_step(s.c, s.h, dW, s.k, w);
    else {   // 2 RK4_A, 4 EPS_A (action-angle state)
        if (s.kind == 2) aa_rk4(s.I, s.p, s.h); else eps_aa(s.I, s.p, s.I2, s.p2, s.h, omA);
        std::vector<double> Ib = s.I, pb = s.p;
        if (s.aabath) bath_pair_aa(s.I, s.p, s.h, dW, s.k);
        else { std::vector<cd> c(n); to_c(s.I, s.p, c); bath_pair_cart(c, s.h, dW, s.k, w);
               for (int j=0;j<n;++j) { s.I[j] = std::norm(c[j]); s.p[j] = pb[j] + std::remainder(std::arg(c[j]) - pb[j], 2*PI); } }
        if (s.kind == 4) for (int j=0;j<n;++j) { s.I2[j] += s.I[j]-Ib[j]; s.p2[j] += s.p[j]-pb[j]; }
    }
    ++s.k; }
std::vector<double> Iget(const S& s) { if (s.kind==1||s.kind==3||s.kind==6) { std::vector<double> I(n); for (int j=0;j<n;++j) I[j]=std::norm(s.c[j]); return I; } return s.I; }
int main(int argc, char** argv) {
    n = atoi(argv[1]); T1 = atof(argv[2]); Tn = atof(argv[3]); double tmax = atof(argv[4]); unsigned seed = atoi(argv[5]);
    omC = atof(argv[6]); omA = atof(argv[7]); FILE* out = fopen(argv[8], "w"); if (argc > 9) g = atof(argv[9]);   // g = 0 and T = 0: no baths
    const double h5 = 5e-4; const int SUB = 16; const double hf = h5 / SUB;
    std::mt19937_64 rng(seed); std::normal_distribution<double> N01(0, 1);
    std::exponential_distribution<double> ex(1.0 / std::sqrt(0.5 * (T1 + Tn) / n)); std::uniform_real_distribution<double> un(-PI, PI);
    std::vector<double> I0(n), p0(n); for (int j = 0; j < n; ++j) { I0[j] = ex(rng); p0[j] = un(rng); }
    std::vector<float> Iof_(I0.begin(), I0.end()), pof(p0.begin(), p0.end()); double aold[4] = {0,0,0,0};
    auto mk = [&](int kind, bool ab, double h, int sub) { S s; s.kind=kind; s.aabath=ab; s.h=h; s.sub=sub; s.I=I0; s.p=p0; s.I2=I0; s.p2=p0; s.c.resize(n); to_c(I0,p0,s.c); s.d=s.c; return s; };
    std::vector<S> sc = { mk(6,false,hf,1), mk(1,false,h5,SUB), mk(2,false,h5,SUB), mk(3,false,h5,SUB), mk(4,false,h5,SUB), mk(6,false,h5,SUB), mk(2,true,h5,SUB), mk(4,true,h5,SUB) };
    const char* nm[9] = {"REF","RK4_C","RK4_A","EPS_C","EPS_A","NEW","RK4_A_AB","EPS_A_AB","OLD"};
    fprintf(out, "t"); for (int s=0;s<9;++s) for (int j=1;j<=n;++j) fprintf(out, ",%s_I%d", nm[s], j);
    fprintf(out, ",EPS_C_sep,EPS_A_sep\n");
    Mid w; long nfine = std::llround(tmax / hf);
    for (long k = 0; k < nfine; ++k) {
        double f[4]; for (auto& v : f) v = N01(rng) * std::sqrt(hf);
        for (int q = 0; q < 4; ++q) aold[q] += f[q];
        for (auto& s : sc) { for (int q = 0; q < 4; ++q) s.acc[q] += f[q]; if ((k+1) % s.sub == 0) { step(s, w); for (auto& v : s.acc) v = 0; } }
        if ((k+1) % SUB == 0) { aa_step<float,true>(Iof_, pof, h5, aold, false); for (auto& v : aold) v = 0; }
        if ((k+1) % (2*SUB) == 0) {
            fprintf(out, "%.4f", (k+1)*hf);
            for (auto& s : sc) { auto I = Iget(s); for (double v : I) fprintf(out, ",%.10g", v); }
            for (float v : Iof_) fprintf(out, ",%.8g", v);
            double sC = 0, sA = 0; for (int j=0;j<n;++j) { sC = std::max(sC, std::abs(sc[3].c[j]-sc[3].d[j])); sA = std::max(sA, std::max(std::fabs(sc[4].I[j]-sc[4].I2[j]), std::fabs(sc[4].p[j]-sc[4].p2[j]))); }
            fprintf(out, ",%.4g,%.4g\n", sC, sA); }
    }
    fclose(out); }
