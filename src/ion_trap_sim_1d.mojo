from mojo_utils import *
from math import pi, sin, cos, floor, sqrt, ceil
from python import Python, PythonObject
import time

alias Z = 1.0
alias e = 1.60217883e-19
alias eps0 = 8.854187817e-12
alias M_Yb = 2.8733965e-25
alias nu = 0.25 * 2*pi*1e6
alias w = 10.0 * 2*pi*1e6

fn main() raises:
    fn move_particle(k: Int, x: Float64) -> Float64:
        return x#x*1.035

    var trap_psuedo = IonTrapSim1D(10, PsuedoPotential1D(), MagneticPotential(3e-20), CoulombPotential(), T=1e-4, dt=1e-9)
    trap_psuedo.set_to_equilibrium(1e-16)
    var x_psuedo: List[List[Float64]]
    var x_motion: List[List[Float64]]
    var t: List[Float64]


    var t_start = time.perf_counter_ns()
    x_psuedo,_,_,t,_ = trap_psuedo.sim_er()
    var t_end = time.perf_counter_ns()
    print("sim_er took", (t_end-t_start)/1e9, "seconds")
    var trap_motion = IonTrapSim1D(trap_psuedo^, MicroMotionPotential1D(), MagneticPotential(3e-20), CoulombPotential())
    trap_motion.map_x_0[move_particle]()
    x_motion,_,_,_,_ = trap_motion.sim_er()

    plot_simulation(t, "Psuedopotential", x_psuedo, "Micromotion", x_motion, 20)

#     var errs = Python.evaluate("[]")
#     var dts = Python.evaluate("[]")
#
#     var T = 1e-4
#     var dt_0 = 1e-10
#     var trap_psuedo = IonTrapSim1D(10, PsuedoPotential1D(), MagneticPotential(), CoulombPotential(), T=T, dt=dt_0)
#     var x1: List[List[Float64]]
#     var t1: List[Float64]
#     x1,_,_,t1,_ = trap_psuedo.sim_er()
#
# #     plot_sim_er(List(t1), List(x1), 10, 0)
#
#     for i in range(300):
#         var dt = pow(10.0, (i * 1e-2)) * dt_0
#         var n_tsteps = int(T / dt)
#         var x2: List[List[Float64]]
#         var t2: List[Float64]
#
#         trap_psuedo.dt = dt
#         trap_psuedo.n_tsteps = n_tsteps
#
#         print("sim_err on dt,n_tsteps", dt, n_tsteps)
#         x2,_,_,t2,_ = trap_psuedo.sim_er()
# #         plot_sim_er(List(t1, t2), List(x1, x2), 10, 0)
#
#         print("get_err")
#         errs.append(get_err(t1, x1, t2, x2, 0))
#         dts.append(dt)
#
#     var np = Python.import_module("numpy")
#     var plt = Python.import_module("matplotlib.pyplot")
#
#     plt.plot(dts, np.array(errs))
#     plt.xscale("log")
#     plt.yscale("log")
#
#     plt.title("Average %Error in psuedo-potential vs dt")
#     plt.show()

fn plot_sim_er(t: List[List[Float64]], x: List[List[List[Float64]]], scale: Int, k: Int) raises:
    var np = Python.import_module("numpy")
    var plt = Python.import_module("matplotlib.pyplot")

    for i in range(len(t)):
        ts = Python.evaluate("[]")
        xs = Python.evaluate("[]")

        for j in range(len(t[i])/scale):
            ts.append(t[i][j*scale])
            xs.append(x[i][j*scale][k])

        plt.plot(ts, np.array(xs), label=str(i))

    plt.legend()
    plt.show()

fn get_err(t1: List[Float64], x1: List[List[Float64]], t2: List[Float64], x2: List[List[Float64]], k: Int) -> Float64:
    var dt1 = t1[1] - t1[0]
    var dt2 = t2[1] - t2[0]
    var err = 0.0

    for i in range(len(t2)):
        var j = t2[i] / dt1
        var j_r = j % 1.0

        if int(ceil(j)) >= len(t1):
            print("break", err)
            break

        var other_x = x1[int(floor(j))][k] * (1.0 - j_r) + j_r * x1[int(ceil(j))][k]
        var percent_err = abs(other_x - x2[i][k]) / abs(x2[i][k])

#         err = max(err, percent_err)
        err += percent_err / len(t2)

    print("return", err)
    return err

fn plot_simulation(t: List[Float64], name1: String, x1: List[List[Float64]], name2: String, x2: List[List[Float64]], skip: Int):
    try:
        var np = Python.import_module("numpy")
        var plt = Python.import_module("matplotlib.pyplot")

        var dt = t[1]-t[0]
        var n = int((pi / w) / dt)

        print("n:", n, "len(t)", len(t))

        for i_ in range(len(x1[0])/skip):
            var i = i_ * skip
            var xs1 = Python.evaluate("[]")
            var xs2 = Python.evaluate("[]")
            var xs3 = Python.evaluate("[]")
            var ts = Python.evaluate("[]")

            for j in range(len(t)-n+1):
                var x_2: Float64 = 0
                var avr_size = min(n, j)
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

    fn force(self, t: Float64, xs: List[Float64], vs: List[Float64], k: Int, mass: Float64) -> Float64:
        return -mass * (nu*nu) * xs[k]

struct MicroMotionPotential1D(Potential1D):
    fn __init__(inout self):
        pass

    fn __copyinit__(inout self, other: Self):
        pass

    fn force(self, t: Float64, xs: List[Float64], vs: List[Float64], k: Int, mass: Float64) -> Float64:
        return -mass * nu * w * cos(w*t) * xs[k] * sqrt(2.0)

struct MagneticPotential(Potential1D):
    var b: Float64

    fn __init__(inout self):
        self.b = 3e-20
        pass

    fn __init__(inout self, b: Float64):
        self.b = b
        pass

    fn __copyinit__(inout self, other: Self):
        self.b = other.b
        pass

    fn force(self, t: Float64, xs: List[Float64], vs: List[Float64], k: Int, mass: Float64) -> Float64:
        return -self.b * vs[k]

struct CoulombPotential(Potential1D):
    fn __init__(inout self):
        pass

    fn __copyinit__(inout self, other: Self):
        pass

    fn force(self, t: Float64, xs: List[Float64], vs: List[Float64], k: Int, mass: Float64) -> Float64:
        var force = 0.0

        for i in range(len(xs)):
            if i != k:
                force += (1 if xs[i] < xs[k] else -1) * ((Z * Z * e * e) / (4*pi*eps0)) * (1 / ((xs[k] - xs[i]) * (xs[k] - xs[i])))

        return force

struct PartialCoulombPotential(Potential1D):
    # Number of ions to compare against on either side
    var ions: Int

    fn __init__(inout self):
        self.ions = 1000000
        pass

    fn __init__(inout self, num_ions: Int):
        self.ions = num_ions
        pass

    fn __copyinit__(inout self, other: Self):
        self.ions = other.ions
        pass

    fn force(self, t: Float64, xs: List[Float64], vs: List[Float64], k: Int, mass: Float64) -> Float64:
        var force = 0.0
        var min_i = max(k - self.ions, 0)
        var max_i = min(k + self.ions, len(xs) - 1)

        for i in range(min_i, max_i+1):
            if i != k:
                force += (1 if xs[i] < xs[k] else -1) * ((Z * Z * e * e) / (4*pi*eps0)) * (1 / ((xs[k] - xs[i]) * (xs[k] - xs[i])))

        return force

struct NoPotential(Potential1D):
    fn __init__(inout self):
        pass

    fn __copyinit__(inout self, other: Self):
        pass

    fn force(self, t: Float64, xs: List[Float64], vs: List[Float64], k: Int, mass: Float64) -> Float64:
        return 0.0

trait Potential1D(Copyable, Defaultable):
    fn force(self, t: Float64, xs: List[Float64], vs: List[Float64], k: Int, mass: Float64) -> Float64: ...

struct IonTrapSim1D[
    P1: Potential1D = NoPotential,
    P2: Potential1D = NoPotential,
    P3: Potential1D = NoPotential,
    P4: Potential1D = NoPotential,
    P5: Potential1D = NoPotential
]:
    var x_0: List[Float64]
    var v_0: List[Float64]
    var n: Int
    var dt: Float64
    var n_tsteps: Int
    var etol: Float64
    var mass: Float64
    var potential1: P1
    var potential2: P2
    var potential3: P3
    var potential4: P4
    var potential5: P5

    fn __init__(
        inout self,
        num_qubits: Int,
        potential1: P1 = P1(),
        potential2: P2 = P2(),
        potential3: P3 = P3(),
        potential4: P4 = P4(),
        potential5: P5 = P5(),
        n_tsteps: Int = 10000,
        dt: Float64 = 1e-8,
        mass: Float64 = M_Yb,
    ):
        self.potential1 = potential1
        self.potential2 = potential2
        self.potential3 = potential3
        self.potential4 = potential4
        self.potential5 = potential5

        self.x_0 = List[Float64](capacity=num_qubits)
        self.v_0 = List[Float64](capacity=num_qubits)
        self.n = num_qubits

        for i in range(num_qubits):
            self.x_0.append(6e-5 * (i - (num_qubits - 1.0) / 2.0) / sqrt(num_qubits * 1.0))
            self.v_0.append(0)

        self.n_tsteps = n_tsteps
        self.dt = dt
        self.etol = 1e-9
        self.mass = mass

    fn __init__(
        inout self,
        num_qubits: Int,
        potential1: P1 = P1(),
        potential2: P2 = P2(),
        potential3: P3 = P3(),
        potential4: P4 = P4(),
        potential5: P5 = P5(),
        T: Float64 = 1e-3,
        dt: Float64 = 1e-8,
        mass: Float64 = M_Yb,
    ):
        self = Self(
            num_qubits,
            potential1,
            potential2,
            potential3,
            potential4,
            potential5,
            int(T / dt),
            dt,
            mass
        )

    fn __init__[P_1: Potential1D, P_2: Potential1D, P_3: Potential1D, P_4: Potential1D, P_5: Potential1D](
        inout self,
        owned other: IonTrapSim1D[P_1, P_2, P_3, P_4, P_5],
        potential1: P1 = P1(),
        potential2: P2 = P2(),
        potential3: P3 = P3(),
        potential4: P4 = P4(),
        potential5: P5 = P5(),
    ):
        self.potential1 = potential1
        self.potential2 = potential2
        self.potential3 = potential3
        self.potential4 = potential4
        self.potential5 = potential5
        self.x_0 = other.x_0
        self.v_0 = other.v_0
        self.n = other.n
        self.dt = other.dt
        self.n_tsteps = other.n_tsteps
        self.etol = other.etol
        self.mass = other.mass

    fn set_to_equilibrium(inout self, rate: Float64):
        var max_acc = 10000.0
        var x = init_vector[Float64](self.n, 0.0)
        try:
            var xs = Python.evaluate("[]")
            var ys = init_vector[PythonObject](self.n, Python.evaluate("[]"))
            var np = Python.import_module("numpy")
            var plt = Python.import_module("matplotlib.pyplot")

            for k in range(len(ys)):
                ys[k] = Python.evaluate("[]")

            for i in range(100000):
                xs.append(i)
                var min_acc = 0.0
                for k in range(self.n):
                    var acc = self.force(0, self.x_0, self.v_0, k) / self.mass
                    min_acc = max(abs(acc), min_acc)
                    x[k] = acc * rate + self.x_0[k]
                    ys[k].append(x[k])

                max_acc = min(min_acc, max_acc)
                self.x_0 = x
#                 if max_acc < 1e-100:
#                     break
#             print(xs, ys[0])

            for k in range(self.n):
                plt.plot(xs, np.array(ys[k]))

            for k in range(self.n):
                print(ys[k][len(ys[k])-1])

            plt.show()
        except e:
            print(e)
            pass

    fn map_x_0[map: fn(Int, Float64) -> Float64](inout self):
        for k in range(self.n):
            self.x_0[k] = map(k, self.x_0[k])

    fn force(self, t: Float64, x: List[Float64], v: List[Float64], k: Int) -> Float64:
        return self.potential1.force(t, x, v, k, self.mass) +
            self.potential2.force(t, x, v, k, self.mass) +
            self.potential3.force(t, x, v, k, self.mass) +
            self.potential4.force(t, x, v, k, self.mass) +
            self.potential5.force(t, x, v, k, self.mass)

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
        print("sim_er using dt:", self.dt, " nt:", self.n_tsteps)

        for k in range(self.n):
            x[0][k] = self.x_0[k]
            v[0][k] = self.v_0[k]
            a[0][k] = self.force(0, self.x_0, self.v_0, k) / self.mass

        for i in range(self.n_tsteps):
            thalf = t[i] + self.dt/2

            for k in range(self.n):
                xhalf[k] = x[i][k] + v[i][k] * (self.dt/2)
                vhalf[k] = v[i][k] + a[i][k] * (self.dt/2)

            for k in range(self.n):
                ahalf = self.force(thalf, xhalf, vhalf, k) / self.mass
                x[i+1][k] = x[i][k] + vhalf[k] * self.dt
                v[i+1][k] = v[i][k] + ahalf * self.dt
                xerr = abs( ( (v[i][k] - vhalf[k])*self.dt ) / 2)
                verr = abs( ( (a[i][k] - ahalf)*self.dt ) / 2)
                err[i][k] = max(xerr, verr)

            for k in range(self.n):
                a[i+1][k] = self.force(t[i], x[i+1], v[i+1], k) / self.mass

            t[i+1] = t[i] + self.dt

            maxerr = 0
            for k in range(self.n):
                maxerr = max(maxerr, err[i][k])

            # dt = 0.9 * sqrt(etol / maxerr) * dt;

        return (x, v, a, t, err)
