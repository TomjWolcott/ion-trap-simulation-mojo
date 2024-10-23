from mojo_utils import *
from math import pi, sin, cos, floor, sqrt
from python import Python

alias Z = 1.0
alias e = 1.60217883e-19
alias eps0 = 8.854187817e-12
alias M_Yb = 2.8733965e-25
alias nu = 0.25 * 2*pi*1e6
alias w = 30*nu

fn main():
    var trap_psuedo = IonTrapSim1D(PsuedoPotential1D(), 3, n_tsteps = 30000, dt=1e-9)
    var x_psuedo: List[List[Float64]]
    var x_motion: List[List[Float64]]
    var t: List[Float64]
    x_psuedo,_,_,t,_ = trap_psuedo.sim_er()
    var trap_motion = IonTrapSim1D(MicroMotionPotential1D(), trap_psuedo^)
    x_motion,_,_,_,_ = trap_motion.sim_er()

    plot_simulation(t, "Psuedopotential", x_psuedo, "Micromotion", x_motion)

fn plot_simulation(t: List[Float64], name1: String, x1: List[List[Float64]], name2: String, x2: List[List[Float64]]):
    try:
        var np = Python.import_module("numpy")
        var plt = Python.import_module("matplotlib.pyplot")

        var dt = t[1]-t[0]
        var n = int((pi / w) / dt)

        print("n:", n, "len(t)", len(t))

        for i in range(len(x1[0])):
            var xs1 = Python.evaluate("[]")
            var xs2 = Python.evaluate("[]")
            var xs3 = Python.evaluate("[]")
            var ts = Python.evaluate("[]")

            for j in range(n, len(t)-n-1):
                var x_2: Float64 = 0
                var avr_size = n#min(j, min(int(floor(n/2)), len(t) - j))
                for k in range(-avr_size, avr_size+1):
                    x_2 += x2[int(j+k-(n-1)/2)][i] / (2*avr_size+1)

                xs1.append(x1[j][i])
                xs2.append(x_2)
                xs3.append(x2[j][i])
                ts.append(t[j])

            plt.plot(ts, np.array(xs1), label=str(name1).__add__(str(" - q_")).__add__(str(i+1)))
            plt.plot(ts, np.array(xs2), label=str(name2).__add__(str(" - q_")).__add__(str(i+1)))
            plt.plot(ts, np.array(xs3), label=str("unaveraged").__add__(str(" - q_")).__add__(str(i+1)))

        plt.legend()
        plt.show()
    except e:
        print("Error", e)

struct PsuedoPotential1D(Potential1D):
    fn __init__(inout self):
        pass

    fn __copyinit__(inout self, other: Self):
        pass

    fn acceleration(self, x: Float64, t: Float64) -> Float64:
        return -(nu*nu) * x

struct MicroMotionPotential1D(Potential1D):
    fn __init__(inout self):
        pass

    fn __copyinit__(inout self, other: Self):
        pass

    fn acceleration(self, x: Float64, t: Float64) -> Float64:
        return -nu * w * -cos(w*t) * x * sqrt(2.0)


trait Potential1D(Copyable):
    fn acceleration(self, x: Float64, t: Float64) -> Float64: ...

struct IonTrapSim1D[P: Potential1D]:
    var x_0: List[Float64]
    var v_0: List[Float64]
    var n: Int
    var dt: Float64
    var b: Float64
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
        mass: Float64 = M_Yb,
        b: Float64 = 3e-20
    ):
        self.potential = potential

        self.x_0 = List[Float64](capacity=num_qubits)
        self.v_0 = List[Float64](capacity=num_qubits)
        self.n = num_qubits

        for i in range(num_qubits):
            self.x_0.append(6e-6 * (i - (num_qubits - 1.0) / 2.0) + 6e-7)
            self.v_0.append(0)

        self.n_tsteps = n_tsteps
        self.dt = dt
        self.etol = 1e-9
        self.mass = mass
        self.b = b

    fn __init__[P2: Potential1D](
        inout self,
        potential: P,
        owned other: IonTrapSim1D[P2]
    ):
        self.potential = potential
        self.x_0 = other.x_0
        self.v_0 = other.v_0
        self.n = other.n
        self.dt = other.dt
        self.b = other.b
        self.n_tsteps = other.n_tsteps
        self.etol = other.etol
        self.mass = other.mass

    fn force(self, t: Float64, x: List[Float64], v: List[Float64], k: Int, b: Float64) -> Float64:
        var force = M_Yb * self.potential.acceleration(x[k], t) - b*v[k]

        for i in range(self.n):
            if i != k:
                force += (1 if x[i] < x[k] else -1) * ((Z * Z * e * e) / (4*pi*eps0)) * (1 / ((x[k] - x[i]) * (x[k] - x[i])))

        return force

    fn sim_er(
        self
    ) -> (List[List[Float64]], List[List[Float64]], List[List[Float64]], List[Float64], List[List[Float64]]):
        var x = init_vector[List[Float64]](self.n_tsteps + 1, init_vector[Float64](self.n, 0.0))
        var v = init_vector[List[Float64]](self.n_tsteps + 1, init_vector[Float64](self.n, 0.0))
        var a = init_vector[List[Float64]](self.n_tsteps + 1, init_vector[Float64](self.n, 0.0))
        var t = init_vector[Float64](self.n_tsteps + 1, 0.0)
        var err = init_vector[List[Float64]](self.n_tsteps + 1, init_vector[Float64](self.n, 0.0))
        var xhalf = init_vector[Float64](self.n, 0.0)
        var vhalf = init_vector[Float64](self.n, 0.0)
        var ahalf: Float64
        var xerr: Float64
        var verr: Float64
        var maxerr: Float64

        for k in range(self.n):
            x[0][k] = self.x_0[k]
            v[0][k] = self.v_0[k]
            a[0][k] = self.force(0, self.x_0, self.v_0, k, self.b) / self.mass

        for i in range(self.n_tsteps):
            thalf = t[i] + self.dt/2

            for k in range(self.n):
                xhalf[k] = x[i][k] + v[i][k] * (self.dt/2)
                vhalf[k] = v[i][k] + a[i][k] * (self.dt/2)

            for k in range(self.n):
                ahalf = self.force(thalf, xhalf, vhalf, k, self.b) / self.mass
                x[i+1][k] = x[i][k] + vhalf[k] * self.dt
                v[i+1][k] = v[i][k] + ahalf * self.dt
                xerr = abs( ( (v[i][k] - vhalf[k])*self.dt ) / 2)
                verr = abs( ( (a[i][k] - ahalf)*self.dt ) / 2)
                err[i][k] = max(xerr, verr)

            for k in range(self.n):
                a[i+1][k] = self.force(t[i], x[i+1], v[i+1], k, self.b) / self.mass

            t[i+1] = t[i] + self.dt

            maxerr = 0
            for k in range(self.n):
                maxerr = max(maxerr, err[i][k])

            # dt = 0.9 * sqrt(etol / maxerr) * dt;

        return (x, v, a, t, err)
