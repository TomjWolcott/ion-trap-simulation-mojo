from mojo_utils import *
from math import pi, sin, cos, sqrt, exp
from python import Python

alias Z = 1.0
alias e = 1.60217883e-19
alias eps0 = 8.854187817e-12
alias M_Yb = 2.8733965e-25
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

# run-time functions aren't allowed in alias definitions, so sqrt_3 is here pre-calculated
alias sqrt_3: Float64 = 1.7320508075688772935274463415058723669428
alias kvector: Vec3 = Vec3( -k/sqrt_3, -k/sqrt_3, -k/sqrt_3, 0 )
alias laserWidth: Float64 = 1e-6

alias decayRate: Float64 = 19.6 * 1e6 * 2*pi
alias detuning: Float64 = -decayRate / 2
alias s: Float64 = 1

fn main():
    var trap_psuedo = IonTrapSim3D(PsuedoPotential3D(), 5, n_tsteps = 10000, dt=1e-9)
    var r_psuedo: List[List[Vec3]]
    var r_motion: List[List[Vec3]]
    var t: List[Float64]

    r_psuedo,_,_,t,_ = trap_psuedo.sim_er()
    var trap_motion = IonTrapSim3D(MicroMotionPotential3D(), trap_psuedo^)
    r_motion,_,_,_,_ = trap_motion.sim_er()

    print(r_motion[100].__str__())
    plot_simulation(t, "Psuedopotential", r_psuedo, "Micromotion", r_motion)

fn plot_simulation(t: List[Float64], name1: String, x1: List[List[Vec3]], name2: String, x2: List[List[Vec3]]):
    try:
        var np = Python.import_module("numpy")
        var plt = Python.import_module("matplotlib.pyplot")
        var n = 1

        for i in range(len(x1[0])):
            var xs1 = Python.evaluate("[]")
            var xs2 = Python.evaluate("[]")
            var ts = Python.evaluate("[]")

            for j in range(n/2, len(t)-n/2):
                var x_2: Float64 = 0
                for k in range(n):
                    x_2 += x2[int(j+k-(n-1)/2)][i][2] / n

                xs1.append(x1[j][i][2])
                xs2.append(x_2)
                ts.append(t[j])

            plt.plot(ts, np.array(xs1), label=str(name1).__add__(str(" - q_")).__add__(str(i)))
            plt.plot(ts, np.array(xs2), label=str(name2).__add__(str(" - q_")).__add__(str(i)))

        plt.legend()
        plt.show()
    except:
        pass

struct PsuedoPotential3D(Potential3D):
    fn __init__(inout self):
        pass

    fn __copyinit__(inout self, other: Self):
        pass

    fn acceleration(self, r: Vec3, t: Float64) -> Vec3:
        return -w*w * r

struct MicroMotionPotential3D(Potential3D):
    fn __init__(inout self):
        pass

    fn __copyinit__(inout self, other: Self):
        pass

    fn acceleration(self, x: Vec3, t: Float64) -> Vec3:
        return -4*40*w * e * cos(40*w*t) / M_Yb * x


trait Potential3D(Copyable):
    fn acceleration(self, x: Vec3, t: Float64) -> Vec3: ...

fn dist_to_laser(r: List[Vec3], k: Int) -> Float64:
    var rminuso = r[k] - laserOrigin
    var kvector_length = sqrt(dot(kvector, kvector));
    var cross_product_length = length(cross(rminuso, kvector))

    return cross_product_length / kvector_length

struct IonTrapSim3D[P: Potential3D]:
    var r_0: List[Vec3]
    var v_0: List[Vec3]
    var n: Int
    var dt: Float64
    var n_tsteps: Int
    var etol: Float64
    var mass: Float64
    var potential: P

    fn __init__(
        inout self,
        potential: P,
        num_qubits: Int,
        n_tsteps: Int = 10000,
        dt: Float64 = 1e-8,
        mass: Float64 = M_Yb
    ):
        self.potential = potential

        self.r_0 = List[Vec3](capacity=num_qubits)
        self.v_0 = List[Vec3](capacity=num_qubits)
        self.n = num_qubits

        for i in range(num_qubits):
            self.r_0.append(vec3(0,0,6e-6 * (i - (num_qubits - 1.0) / 2.0) + 6e-7))
            self.v_0.append(vec3(0,0,0))

        self.n_tsteps = n_tsteps
        self.dt = dt
        self.etol = 1e-9
        self.mass = mass

    fn __init__[P2: Potential3D](
        inout self,
        potential: P,
        owned other: IonTrapSim3D[P2]
    ):
        self.potential = potential
        self.r_0 = other.r_0
        self.v_0 = other.v_0
        self.n = other.n
        self.dt = other.dt
        self.n_tsteps = other.n_tsteps
        self.etol = other.etol
        self.mass = other.mass

    fn force(self, t: Float64, r: List[Vec3], v: List[Vec3], k: Int, laserWidth: Float64) -> Vec3:
        # Harmonic Potential force
        var a = self.mass * self.potential.acceleration(r[k], t)

        # Laser cooling
        var kdotv = dot(kvector, v[k])
        var satParam = s * exp( -2 * pow(dist_to_laser(r, k), 2) / (laserWidth * laserWidth))
        var scattering_rate = (decayRate * satParam) / (1 + 2*satParam + pow( (2*(detuning - kdotv)) / decayRate, 2) )
        a += kvector * (hbar * scattering_rate)

        # Coulomb Force
        for i in range(len(r)):
            if i != k:
                var dri_mag = length(r[k] - r[i])
                var dri = (r[k] - r[i]) / dri_mag
                var ai_mag = ((Z*Z * e*e) / (4 * pi * eps0)) * (1 / pow(dri_mag, 2));

                a += ai_mag * dri;

        return a;

    fn sim_er(
        self
    ) -> (List[List[Vec3]], List[List[Vec3]], List[List[Vec3]], List[Float64], List[List[Float64]]):
        var r = init_vector[List[Vec3]](self.n_tsteps + 1, init_vector[Vec3](self.n, 0.0))
        var v = init_vector[List[Vec3]](self.n_tsteps + 1, init_vector[Vec3](self.n, 0.0))
        var a = init_vector[List[Vec3]](self.n_tsteps + 1, init_vector[Vec3](self.n, 0.0))
        var t = init_vector[Float64](self.n_tsteps + 1, 0.0)
        var err = init_vector[List[Float64]](self.n_tsteps + 1, init_vector[Float64](self.n, 0.0))
        var rhalf = init_vector[Vec3](self.n, 0.0)
        var vhalf = init_vector[Vec3](self.n, 0.0)
        var ahalf: Vec3
        var rerr: Float64
        var verr: Float64
        var maxerr: Float64

        for k in range(self.n):
            r[0][k] = self.r_0[k]
            v[0][k] = self.v_0[k]
            a[0][k] = self.force(0, self.r_0, self.v_0, k, laserWidth) / self.mass

        for i in range(self.n_tsteps):
            thalf = t[i] + self.dt/2

            for k in range(self.n):
                rhalf[k] = r[i][k] + v[i][k] * (self.dt/2)
                vhalf[k] = v[i][k] + a[i][k] * (self.dt/2)

            for k in range(self.n):
                ahalf = self.force(thalf, rhalf, vhalf, k, laserWidth) / self.mass
                r[i+1][k] = r[i][k] + vhalf[k] * self.dt
                v[i+1][k] = v[i][k] + ahalf * self.dt
#                 rerr = abs( ( (v[i][k] - vhalf[k])*self.dt ) / 2)
#                 verr = abs( ( (a[i][k] - ahalf)*self.dt ) / 2)
#                 err[i][k] = max(rerr, verr)

            for k in range(self.n):
                a[i+1][k] = self.force(t[i], r[i+1], v[i+1], k, laserWidth) / self.mass

            t[i+1] = t[i] + self.dt

#             maxerr = 0
#             for k in range(self.n):
#                 maxerr = max(maxerr, err[i][k])

            # dt = 0.9 * sqrt(etol / maxerr) * dt;

        return (r, v, a, t, err)
