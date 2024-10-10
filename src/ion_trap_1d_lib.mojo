from collections.list import List
from mojo_utils import init_vector

alias PI = 3.14159265358979323846
alias Z = 1.0
alias e = 1.60217883e-19
alias eps0 = 8.854187817e-12
alias M_Yb = 2.8733965e-25
alias nu = 0.25 * 2*PI*1e6

alias l0 = pow(((Z*Z * e*e) / (4.0 * PI * eps0 * M_Yb * nu*nu)), 1.0/3.0)
alias m0 = M_Yb
alias t0 = 1.0 / nu

fn total_energy(x: List[Float64]) -> Float64:
    var energy: Float64 = 0
    var n = len(x)

    for i in range(n):
        energy += 0.5 * M_Yb * (nu*nu) * (x[i]*x[i])

    for i in range(n):
        for j in range(i+1, n):
            energy += ((Z*Z * e*e) / (4 * PI * eps0)) * (1 / abs(x[i] - x[j]))

    return energy

fn neg_grad_energy(inout x: List[Float64], k: Int, dx: Float64) -> Float64:
    x[k] -= dx
    var energy1 = total_energy(x)
    x[k] += 2*dx
    energy2 = total_energy(x)
    x[k] -= dx
    return -(energy2 - energy1) / (2*dx)

fn force(x: List[Float64], k: Int) -> Float64:
    var force = -M_Yb * (nu*nu) * x[k]

    for i in range(len(x)):
        if i != k:
            force += (1 if x[i] < x[k] else -1) * ((Z * Z * e * e) / (4*PI*eps0)) * (1 / ((x[k] - x[i]) * (x[k] - x[i])));

    return force

fn force(
    x: List[Float64],
    v: List[Float64],
    k: Int,
    b: Float64
) -> Float64:
    var force = -M_Yb * (nu * nu) * x[k] - b*v[k]

    for i in range(len(x)):
        if i != k:
            force += (1 if x[i] < x[k] else -1) * ((Z * Z * e * e) / (4*PI*eps0)) * (1 / ((x[k] - x[i]) * (x[k] - x[i])))

    return force

fn sim_leapfrog(
    T: Float64,
    dt: Float64,
    M: Float64,
    x_0: List[Float64],
    v_0: List[Float64],
    dx: Float64 = 1e-9
) -> (List[List[Float64]], List[List[Float64]]):
    var n_tsteps: Int = int(T / dt);
    var n = len(x_0)
    var x: List[List[Float64]] = init_vector[List[Float64]](n_tsteps + 1, init_vector[Float64](n, 0.0))
    var v: List[List[Float64]] = init_vector[List[Float64]](n_tsteps + 1, init_vector[Float64](n, 0.0))
    var a: List[List[Float64]] = init_vector[List[Float64]](n_tsteps + 1, init_vector[Float64](n, 0.0))

    for k in range(n):
        x[0][k] = x_0[k]
        v[0][k] = v_0[k]
        a[0][k] = force(x_0, k) / M

    for i in range(n_tsteps):
        for k in range(n):
            x[i+1][k] = x[i][k] + dt*v[i][k] + 0.5 * (dt*dt) * a[i][k]

        for k in range(n):
            a[i+1][k] = force(x[i+1], k) / M
            v[i+1][k] = v[i][k] + 0.5*dt*(a[i][k] + a[i+1][k])

    return (x, v)

fn sim_er(
    n_tsteps: Int,
    dt: Float64,
    etol: Float64,
    M: Float64,
    b: Float64,
    x_0: List[Float64],
    v_0: List[Float64]
) -> (List[List[Float64]], List[List[Float64]], List[List[Float64]], List[List[Float64]], List[List[Float64]]):
    var n = len(x_0)
    var x: List[List[Float64]] = init_vector[List[Float64]](n_tsteps + 1, init_vector[Float64](n, 0.0))
    var v: List[List[Float64]] = init_vector[List[Float64]](n_tsteps + 1, init_vector[Float64](n, 0.0))
    var a: List[List[Float64]] = init_vector[List[Float64]](n_tsteps + 1, init_vector[Float64](n, 0.0))
    var t: List[List[Float64]] = init_vector[List[Float64]](n_tsteps + 1, init_vector[Float64](n, 0.0))
    var err: List[List[Float64]] = init_vector[List[Float64]](n_tsteps + 1, init_vector[Float64](n, 0.0))
    var xhalf: List[Float64] = init_vector[Float64](n, 0.0)
    var vhalf: List[Float64] = init_vector[Float64](n, 0.0)
    var ahalf: Float64
    var xerr: Float64
    var verr: Float64
    var maxerr: Float64

    for k in range(n):
        x[0][k] = x_0[k]
        v[0][k] = v_0[k]
        a[0][k] = force(x_0, v_0, k, b) / M

    for i in range(n_tsteps):
        for k in range(n):
            xhalf[k] = x[i][k] + v[i][k] * (dt/2)
            vhalf[k] = v[i][k] + a[i][k] * (dt/2)

        for k in range(n):
            ahalf = force(xhalf, vhalf, k, b) / M
            x[i+1][k] = x[i][k] + vhalf[k] * dt
            v[i+1][k] = v[i][k] + ahalf * dt
            xerr = abs( ( (v[i][k] - vhalf[k])*dt ) / 2)
            verr = abs( ( (a[i][k] - ahalf)*dt ) / 2)
            err[i][k] = max(xerr, verr)

        for k in range(n):
            a[i+1][k] = force(x[i+1], v[i+1], k, b) / M

        t[i+1][0] = t[i][0] + dt

        maxerr = 0
        for k in range(n):
            maxerr = max(maxerr, err[i][k])

        # dt = 0.9 * sqrt(etol / maxerr) * dt;

    return (x, v, a, t, err)

fn a_dless(x: List[Float64], k: Int) -> Float64:
    var a: Float64 = -x[k]
    var n = len(x)

    for i in range(n):
        if i != k:
            a += (1 if x[i] < x[k] else -1) * (1.0 / ((x[k]-x[i]) * (x[k]-x[i])))

    return a

fn sim_leapfrog_dless(
    owned T: Float64,
    owned dt: Float64,
    owned x_0: List[Float64],
    owned v_0: List[Float64]
) -> (List[List[Float64]], List[List[Float64]]):
    var l0 = pow(((Z*Z * e*e) / (4.0 * PI * eps0 * M_Yb * nu*nu)), 1.0/3.0)

    T = T / t0
    dt = dt / t0
    var n = len(x_0)
    for k in range(n):
        x_0[k] = x_0[k] / l0
        v_0[k] = v_0[k] / (l0 / t0)

    var n_tsteps: Int = int(T / dt)
    var x: List[List[Float64]] = init_vector[List[Float64]](n_tsteps + 1, init_vector[Float64](n, 0.0))
    var v: List[List[Float64]] = init_vector[List[Float64]](n_tsteps + 1, init_vector[Float64](n, 0.0))
    var a: List[List[Float64]] = init_vector[List[Float64]](n_tsteps + 1, init_vector[Float64](n, 0.0))

    for k in range(n):
        x[0][k] = x_0[k]
        v[0][k] = v_0[k]
        a[0][k] = a_dless(x_0, k)

    for i in range(n_tsteps):
        for k in range(n):
            x[i+1][k] = x[i][k] + dt*v[i][k] + 0.5*(dt*dt)*a[i][k]
        for k in range(n):
            a[i+1][k] = a_dless(x[i+1], k)
            v[i+1][k] = v[i][k] + 0.5*dt*(a[i][k] + a[i+1][k])

    for i in range(n_tsteps+1):
        for k in range(n):
            x[i][k] = x[i][k] * l0
            v[i][k] = v[i][k] * (l0 / t0)

    return (x, v)