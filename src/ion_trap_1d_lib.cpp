#define _USE_MATH_DEFINES

#include "ion_trap_1d_lib.h"
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

// physical constants
static const double Z = 1;
static const double e = 1.60217883e-19;
static const double eps0 = 8.854187817e-12;
static const double M_Yb = 2.8733965e-25;
static const double nu = 0.25 * 2*M_PI*1e6;

// dimensionless conversion factors
static const double l0 = pow(((Z*Z * e*e) / (4 * M_PI * eps0 * M_Yb * nu*nu)), 1.0/3.0);
static const double m0 = M_Yb;
static const double t0 = 1.0 / nu;

using namespace std;

double total_energy(
    const int n,
    vector<double>& x
    )
{
    double energy = 0;
    for (int i = 0; i < n; i++) energy += 0.5 * M_Yb * (nu * nu) * (x[i] * x[i]);
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            energy += ((Z * Z * e * e) / (4 * M_PI * eps0)) * (1 / abs(x[i] - x[j]));
    return energy;
} 

double neg_grad_energy(
    const int n,
    vector<double>& x,
    const int k,
    const double dx
    )
{
    x[k] = x[k] - dx;
    double energy1 = total_energy(n, x);
    x[k] = x[k] + 2*dx;
    double energy2 = total_energy(n, x);
    x[k] = x[k] - dx;
    return -(energy2 - energy1) / (2*dx);
}

double force(
    const int n,
    vector<double>& x,
    const int k
    )
{
    double force = -M_Yb * (nu * nu) * x[k];
    for (int i = 0; i < n; i++)
        if (i != k)
            force += (x[i] < x[k] ? 1 : -1) * ((Z * Z * e * e) / (4*M_PI*eps0)) * (1 / ((x[k] - x[i]) * (x[k] - x[i])));
    return force;
}

double force(
    const int n,
    vector<double>& x,
    vector<double>& v,
    const int k,
    const double b
    )
{
    double force = -M_Yb * (nu * nu) * x[k] - b*v[k];
    for (int i = 0; i < n; i++)
        if (i != k)
            force += (x[i] < x[k] ? 1 : -1) * ((Z * Z * e * e) / (4*M_PI*eps0)) * (1 / ((x[k] - x[i]) * (x[k] - x[i])));
    return force;
}


pair<vector<vector<double> >, vector<vector<double> > > sim_leapfrog(
    const int n,
    const double T,
    const double dt,
    const double M,
    vector<double>& x_0,
    vector<double>& v_0,
    const double dx
    )
{
    const int n_tsteps = T / dt;
    vector<vector<double> > x(n_tsteps + 1, vector<double>(n));
    vector<vector<double> > v(n_tsteps + 1, vector<double>(n));
    vector<vector<double> > a(n_tsteps + 1, vector<double>(n));

    for (int k = 0; k < n; k++) {
        x[0][k] = x_0[k];
        v[0][k] = v_0[k];
        a[0][k] = force(n, x_0, k) / M;
    }

    for (int i = 0; i < n_tsteps; i++) {
        for (int k = 0; k < n; k++)
            x[i+1][k] = x[i][k] + dt*v[i][k] + 0.5*(dt*dt)*a[i][k];
        for (int k = 0; k < n; k++) {
            a[i+1][k] = force(n, x[i+1], k) / M;
            v[i+1][k] = v[i][k] + 0.5*dt*(a[i][k] + a[i+1][k]);
        }
    }

    pair<vector<vector<double> >, vector<vector<double> > > ret = make_pair(x, v);
    return ret;
}

vector<vector<vector<double > > > sim_er(
    const int n,
    const int n_tsteps,
    double dt,
    double etol,
    const double M,
    const double b,
    vector<double>& x_0,
    vector<double>& v_0
    )
{
    vector<vector<double> > x(n_tsteps + 1, vector<double>(n, 0));
    vector<vector<double> > v(n_tsteps + 1, vector<double>(n, 0));
    vector<vector<double> > a(n_tsteps + 1, vector<double>(n, 0));
    vector<vector<double> > t(n_tsteps + 1, vector<double>(n ,0));
    vector<vector<double> > err(n_tsteps + 1, vector<double>(n, 0));
    vector<double> xhalf(n, 0);
    vector<double> vhalf(n, 0);
    double ahalf;
    double xerr;
    double verr;
    double maxerr;

    for (int k = 0; k < n; k++) {
        x[0][k] = x_0[k];
        v[0][k] = v_0[k];
        a[0][k] = force(n, x_0, v_0, k, b) / M;
    }

    for (int i = 0; i < n_tsteps; i++) {
        for (int k = 0; k < n; k++) {
            xhalf[k] = x[i][k] + v[i][k] * (dt/2);
            vhalf[k] = v[i][k] + a[i][k] * (dt/2);
        }
        for (int k = 0; k < n; k++) {
            ahalf = force(n, xhalf, vhalf, k, b) / M;
            x[i+1][k] = x[i][k] + vhalf[k]*dt;
            v[i+1][k] = v[i][k] + ahalf*dt;
            xerr = abs( ( (v[i][k] - vhalf[k])*dt ) / 2);
            verr = abs( ( (a[i][k] - ahalf)*dt ) / 2);
            err[i][k] = max(xerr, verr);
        }
        for (int k = 0; k < n; k++)
            a[i+1][k] = force(n, x[i+1], v[i+1], k, b) / M;
        
        t[i+1][0] = t[i][0] + dt;
        maxerr = 0;
        for (int k = 0; k < n; k++) {
            maxerr = max(maxerr, err[i][k]);
        }
        // dt = 0.9 * sqrt(etol / maxerr) * dt;
    }

    vector<vector<vector<double> > > ret(5);
    ret[0] = x;
    ret[1] = v;
    ret[2] = a;
    ret[3] = t;
    ret[4] = err;
    return ret;
}

double a_dless(
    const int n,
    vector<double>& x, 
    const int k
    )
{
    double a = -x[k];
    for (int i = 0; i < n; i++)
        if (i != k)
            a += (x[i] < x[k] ? 1 : -1) * (1.0 / ((x[k]-x[i]) * (x[k]-x[i])));
    return a;
}

pair<vector<vector<double> >, vector<vector<double> > > sim_leapfrog_dless(
        const int n,
        double T,
        double dt, 
        vector<double>& x_0,
        vector<double>& v_0
    )
{
    T = T / t0;
    dt = dt / t0;
    for (int k = 0; k < n; k++) {
        x_0[k] = x_0[k] / l0;
        v_0[k] = v_0[k] / (l0 / t0);
    }

    const int n_tsteps = T / dt;
    vector<vector<double> > x(n_tsteps + 1, vector<double>(n));
    vector<vector<double> > v(n_tsteps + 1, vector<double>(n));
    vector<vector<double> > a(n_tsteps + 1, vector<double>(n));

    for (int k = 0; k < n; k++) {
        x[0][k] = x_0[k];
        v[0][k] = v_0[k];
        a[0][k] = a_dless(n, x_0, k);
    }

    for (int i = 0; i < n_tsteps; i++) {
        for (int k = 0; k < n; k++)
            x[i+1][k] = x[i][k] + dt*v[i][k] + 0.5*(dt*dt)*a[i][k];
        for (int k = 0; k < n; k++) {
            a[i+1][k] = a_dless(n, x[i+1], k);
            v[i+1][k] = v[i][k] + 0.5*dt*(a[i][k] + a[i+1][k]);
        }
    }

    for (int i = 0; i <= n_tsteps; i++) {
        for (int k = 0; k < n; k++) {
            x[i][k] = x[i][k] * l0;
            v[i][k] = v[i][k] * (l0 / t0);
        }
    }

    pair<vector<vector<double> >, vector<vector<double> > > ret = make_pair(x, v);
    return ret;
}
