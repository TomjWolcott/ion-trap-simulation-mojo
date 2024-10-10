from math import pi, sqrt, exp
from mojo_utils import *

# physical constants
alias Z: Float64 = 1
alias e: Float64 = 1.60217883e-19
alias eps0: Float64 = 8.854187817e-12
alias M_Yb: Float64 = 2.8733965e-25
alias hbar: Float64 = 1.05457182e-34

# NOTE: simd width must be a power of 2, and runtime functions
# are not allowed in aliases, so just ignore the extra zero on all
# of the Vec3's they're never relevant

# trap parameters
alias w: Vec3 = Vec3(5.7 * 2*pi*1e6, 5.7 * 2*pi*1e6, 1.5 * 2*pi*1e6, 0)
# alias wx: Float64 = 5.7 * 2*pi*1e6
# alias wy: Float64 = 5.7 * 2*pi*1e6
# alias wz: Float64 = 1.5 * 2*pi*1e6

# laser cooling parameters
alias laserWavelength: Float64 = 369 * 1e-9
alias k: Float64 = (2*pi) / laserWavelength
alias laserOrigin: Vec3 = Vec3(0, 0, 0, 0)
alias sqrt_3: Float64 = 1.7320508075688772935274463415058723669428
alias kvector: Vec3 = Vec3( -k/sqrt_3, -k/sqrt_3, -k/sqrt_3, 0 )
alias kvector_length: Float64 = k # Update if kvector is going to regularly change
alias laserWidth: Float64 = 1e-6

alias decayRate: Float64 = 19.6 * 1e6 * 2*pi
alias detuning: Float64 = -decayRate / 2
alias s: Float64 = 1

fn acceleration(
    r: List[Vec3],
    k: Int
) -> Vec3:
    var a = -w*w * r[k]

    # Coulomb Force
    for i in range(len(r)):
        if i != k:
            var dri_mag = length(r[k] - r[i])
            var dri = (r[k] - r[i]) / dri_mag
            var ai_mag = ((Z*Z * e*e) / (4 * pi * eps0 * M_Yb)) * (1 / (dri_mag*dri_mag));

            a += ai_mag * dri;

    return a;

fn acceleration(
    r: List[Vec3],
    v: List[Vec3],
    k: Int
) -> Vec3:
    var a = -w*w * r[k]

    # Laser cooling
    var kdotv = dot(kvector, v[k])
    var rhoee = (s/2) / (1 + s + pow( (2*(detuning - kdotv)) / decayRate, 2) )
    a += kvector * (hbar * decayRate * rhoee) / M_Yb

    # Coulomb Force
    for i in range(len(r)):
        if i != k:
            var dri_mag = length(r[k] - r[i])
            var dri = (r[k] - r[i]) / dri_mag
            var ai_mag = ((Z*Z * e*e) / (4 * pi * eps0 * M_Yb)) * (1 / (dri_mag*dri_mag));

            a += ai_mag * dri;

    return a;

fn dist_to_laser(r: List[Vec3], k: Int) -> Float64:
    var rminuso = r[k] - laserOrigin
    var cross_product_length = length(cross(rminuso, kvector))

    return cross_product_length / kvector_length

fn acceleration(
    r: List[Vec3],
    v: List[Vec3],
    k: Int,
    laserWidth: Float64
) -> Vec3:
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
            var ai_mag = ((Z*Z * e*e) / (4 * pi * eps0 * M_Yb)) * (1 / (dri_mag*dri_mag));

            a += ai_mag * dri;

    return a;


fn sim_leapfrog(
    T: Float64,
    dt: Float64,
    r_0: List[Vec3],
    v_0: List[Vec3]
) -> (List[List[Vec3]], List[List[Vec3]]):
    var n_tsteps: Int = int(T / dt);
    var n = len(r_0)
    var r = init_vector[List[Vec3]](n_tsteps+1, init_vector[Vec3](n, vec3(0, 0, 0)))
    var v = init_vector[List[Vec3]](n_tsteps+1, init_vector[Vec3](n, vec3(0, 0, 0)))
    var vhalf = init_vector[List[Vec3]](n_tsteps+1, init_vector[Vec3](n, vec3(0, 0, 0)))
    var a = init_vector[List[Vec3]](n_tsteps+1, init_vector[Vec3](n, vec3(0, 0, 0)))

    for k in range(n):
        r[0][k] = r_0[k]
        v[0][k] = v_0[k]
        a[0][k] = acceleration(r_0, k)

    for i in range(n_tsteps):
        for k in range(n):
            vhalf[i][k] = v[i][k] + 0.5*dt*a[i][k]
            r[i+1][k] = r[i][k] + dt*vhalf[i][k]

        for k in range(n):
            a[i+1][k] = acceleration(r[i+1], k)
            v[i+1][k] = vhalf[i][k] + 0.5*dt*a[i+1][k]

    return (r, v)

fn sim_er(
    n_tsteps: Int,
    owned dt: Float64,
    etol: Float64,
    r_0: List[Vec3],
    v_0: List[Vec3]
) -> (List[List[Vec3]], List[List[Vec3]], List[List[Vec3]], List[Float64], List[Float64]):
    var n = len(r_0)
    var r = init_vector[List[Vec3]](n_tsteps + 1, init_vector[Vec3](n, vec3(0, 0, 0)))
    var v = init_vector[List[Vec3]](n_tsteps + 1, init_vector[Vec3](n, vec3(0, 0, 0)))
    var a = init_vector[List[Vec3]](n_tsteps + 1, init_vector[Vec3](n, vec3(0, 0, 0)))
    var t = init_vector[Float64](n_tsteps + 1, 0.0)
    var err = init_vector[Float64](n_tsteps + 1, 0.0)
    var ahalf: Vec3
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