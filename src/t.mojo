from mojo_utils import *
from math import pi


fn psuedopotential_acc_laser[dims: Int](
    x: List[SIMD[DType.float64, dims]],
    v: List[SIMD[DType.float64, dims]],
    k: Int,
    laserWidth: Float64
):
    # Harmonic Potential force
    var a = -w*w * r[k]

    # Laser cooling
    var kdotv = dot(kvector, v[k])
    var satParam = s * exp( -2 * pow(dist_to_laser(r, k), 2) / (laserWidth * laserWidth))
    var scattering_rate = (decayRate * satParam) / (1 + 2*satParam + pow( (2*(detuning - kdotv)) / decayRate, 2) )
    a += kvector * (hbar * scattering_rate) / M_Yb

    # Coulomb Force
    for i in range(len(r)):
        if i != k:
            var dri_mag = length(r[k] - r[i])
            var dri = (r[k] - r[i]) / dri_mag
            var ai_mag = ((Z*Z * e*e) / (4 * pi * eps0 * M_Yb)) * (1 / pow(dri_mag, 2));

            a += ai_mag * dri;

    return a;

fn sim_er[
    dims: Int,
    acceleration: fn(
        List[SIMD[DType.float64, dims]],
        List[SIMD[DType.float64, dims]],
        Int,
        Float64
    ) -> SIMD[DType.float64, dims]
](
    n_tsteps: Int,
    owned dt: Float64,
    etol: Float64,
    r_0: List[SIMD[DType.float64, dims]],
    v_0: List[SIMD[DType.float64, dims]],
) -> (
    List[SIMD[DType.float64, dims]],
    List[SIMD[DType.float64, dims]],
    List[SIMD[DType.float64, dims]],
    List[Float64], List[Float64]
):
    alias VecN = SIMD[DType.float64, dims]

    var n = len(xs_0)
    var r = init_vector[List[VecN]](n_tsteps + 1, init_vector[VecN](n, vec3(0, 0, 0)))
    var v = init_vector[List[VecN]](n_tsteps + 1, init_vector[VecN](n, vec3(0, 0, 0)))
    var a = init_vector[List[VecN]](n_tsteps + 1, init_vector[VecN](n, vec3(0, 0, 0)))
    var t = init_vector[Float64](n_tsteps + 1, 0.0)
    var err = init_vector[Float64](n_tsteps + 1, 0.0)
    var ahalf: VecN
    var rerr: Float64
    var verr: Float64
    var maxerr: Float64

    for k in range(n):
        r[0][k] = r_0[k]
        v[0][k] = v_0[k]
        a[0][k] = acceleration(r_0, v_0, k, laserWidth)

    var rhalf = init_vector[Vec3](n, vec3(0, 0, 0))
    var vhalf = init_vector[Vec3](n, vec3(0, 0, 0))

    for i in range(n_tsteps):
        err[i+1] = etol + 1;
        while err[i+1] >= etol:
            err[i+1] = 0 # reset error
            for k in range(n):
                rhalf[k] = r[i][k] + v[i][k]*(dt/2)
                vhalf[k] = v[i][k] + a[i][k]*(dt/2)

            for k in range(n):
                for j in range(3):
                    rerr = abs( ( (v[i][k][j] - vhalf[k][j])*dt ) / 2);
                    err[i+1] = max(err[i+1], rerr);

            if err[i+1] < etol:
                for k in range(n):
                    ahalf = acceleration(rhalf, vhalf, k, laserWidth);
                    r[i+1][k] = r[i][k] + vhalf[k]*dt;
                    v[i+1][k] = v[i][k] + ahalf*dt;

                for k in range(n):
                    a[i+1][k] = acceleration(r[i+1], v[i+1], k, laserWidth);

                t[i+1] = t[i] + dt;
            else:
                dt = 0.9 * sqrt(etol / err[i+1]) * dt;

    return (r, v, a, t, err)