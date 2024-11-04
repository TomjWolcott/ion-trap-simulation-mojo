from mojo_utils import *
from math import pi, sin, cos, sqrt, exp, floor
from python import Python, PythonObject
from potentials_3d import *
import time

alias Z = 1.0
alias e = 1.60217883e-19
alias eps0 = 8.854187817e-12
alias M_Yb = 2.8733965e-25
alias hbar: Float64 = 1.05457182e-34

fn main():
    fn map_x0(i: Int, x: Vec3) -> Vec3:
        return x / 1.1

    var ion_trap_sim = IonTrapSim3D(2, RingPsuedoPotential(), CoulombPotential(), T=1e-6, dt=1e-9)
    ion_trap_sim.set_potential_equilibrium(ion_trap_sim.potential1)
    ion_trap_sim.optimize_equilibrium(1e-18)
    ion_trap_sim.map_x_0[map_x0]()

    var xs: List[List[Vec3]]

    xs,_,_,_,_ = ion_trap_sim.sim_er()

    animated_plot(xs, 1000)

fn animated_plot(points: List[List[Vec3]], nt: Int):
    try:
        Python.add_to_path("./")
        animation_callback = Python.import_module("animation_callback")

        var np = Python.import_module("numpy")
        var plt = Python.import_module("matplotlib.pyplot")
        var animation = Python.import_module("matplotlib.animation")
        var fig = plt.figure()
        var ax = fig.add_subplot(projection='3d')
        var line: PythonObject
        var xs = Python.evaluate("[]")
        var ys = Python.evaluate("[]")
        var zs = Python.evaluate("[]")
        var artists = Python.evaluate("[]")

        for i in range(nt):
            var j = i * int(floor((len(points) - 1) / (nt - 1)))
            var x = Python.evaluate("[]")
            var y = Python.evaluate("[]")
            var z = Python.evaluate("[]")
            for k in range(len(points[j])):
                x.append(points[j][k][0])
                y.append(points[j][k][1])
                z.append(points[j][k][2])

            xs.append(x)
            ys.append(y)
            zs.append(z)
            artists.append(Python.evaluate("lambda x: [x]")(ax.scatter(xs[i], ys[i], zs[i], c="red")))

#         var callback = animation_callback.callback(xs, ys, zs, np, ax)
#
#         print(callback(2))
#
#         var anim = animation.FuncAnimation(
#             fig,
#             callback,
#             Python.evaluate("lambda nt: range(nt)")(nt),
#             interval=100,
#             repeat=True,
#             blit=True
#         )

        var anim = animation.ArtistAnimation(
            fig,
            artists,
            20,
            1000,
            repeat=True
        )

        plt.show()
        print(anim)
    except e:
        print(e)

fn convert_list_to_python(xs: List[Vec3]) raises -> PythonObject:
    var xs_py = Python.evaluate("[]")

    for i in range(len(xs)):
        xs_py.append(Python.evaluate("[]"))

        for j in range(3):
            xs_py[i].append(xs[i][j])

    return xs_py

struct IonTrapSim3D[
    P1: Potential3D = NoPotential,
    P2: Potential3D = NoPotential,
    P3: Potential3D = NoPotential,
    P4: Potential3D = NoPotential,
    P5: Potential3D = NoPotential
]:
    var x_0: List[Vec3]
    var v_0: List[Vec3]
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
        num_ions: Int,
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

        self.x_0 = init_vector[Vec3](num_ions, Vec3(0, 0, 0))
        self.v_0 = init_vector[Vec3](num_ions, Vec3(0, 0, 0))
        self.n = num_ions

        self.n_tsteps = n_tsteps
        self.dt = dt
        self.etol = 1e-9
        self.mass = mass

    fn __init__(
        inout self,
        num_ions: Int,
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
            num_ions,
            potential1,
            potential2,
            potential3,
            potential4,
            potential5,
            int(T / dt),
            dt,
            mass
        )

    fn __init__[P_1: Potential3D, P_2: Potential3D, P_3: Potential3D, P_4: Potential3D, P_5: Potential3D](
        inout self,
        owned other: IonTrapSim3D[P_1, P_2, P_3, P_4, P_5],
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

    fn set_potential_equilibrium[P: GetEquilibriumPositions](inout self, potential: P):
        self.x_0 = potential.get_equilibrium_positions(self.n)
        self.v_0 = init_vector[Vec3](self.n, Vec3(0, 0, 0))

    fn optimize_equilibrium(inout self, rate: Float64):
        var max_acc = 10000.0
        var x = init_vector[Vec3](self.n, Vec3(0, 0, 0))
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
                    min_acc = max(abs(length(acc)), min_acc)
                    x[k] = acc * rate + self.x_0[k]
                    ys[k].append(length(x[k]))

                max_acc = min(min_acc, max_acc)
                self.x_0 = x
                if max_acc < 1e-100:
                    break
#             print(xs, ys[0])

            for k in range(self.n):
                plt.plot(xs, np.array(ys[k]))

            for k in range(self.n):
                print(ys[k][len(ys[k])-1])

            plt.show()
        except e:
            print(e)
            pass

    fn map_x_0[map: fn(Int, Vec3) -> Vec3](inout self):
        for k in range(self.n):
            self.x_0[k] = map(k, self.x_0[k])

    fn force(self, t: Float64, x: List[Vec3], v: List[Vec3], k: Int) -> Vec3:
        return self.potential1.force(t, x, v, k, self.mass) +
            self.potential2.force(t, x, v, k, self.mass) +
            self.potential3.force(t, x, v, k, self.mass) +
            self.potential4.force(t, x, v, k, self.mass) +
            self.potential5.force(t, x, v, k, self.mass)

    fn sim_er(
        self
    ) -> (List[List[Vec3]], List[List[Vec3]], List[List[Vec3]], List[Float64], List[List[Float64]]):
        var t_start = time.perf_counter_ns()
        var x = init_vector[List[Vec3]](self.n_tsteps + 1, init_vector[Vec3](self.n, vec3(0, 0, 0)))
        var v = init_vector[List[Vec3]](self.n_tsteps + 1, init_vector[Vec3](self.n, vec3(0, 0, 0)))
        var a = init_vector[List[Vec3]](self.n_tsteps + 1, init_vector[Vec3](self.n, vec3(0, 0, 0)))
        var t = init_vector[Float64](self.n_tsteps + 1, 0.0)
        var err = init_vector[List[Float64]](self.n_tsteps + 1, init_vector[Float64](self.n, 0.0))
        var xhalf = init_vector[Vec3](self.n, vec3(0, 0, 0))
        var vhalf = init_vector[Vec3](self.n, vec3(0, 0, 0))
        var ahalf: Vec3
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
                xerr = dot(v[i][k] - vhalf[k],v[i][k] - vhalf[k]) * self.dt / 2
                verr = dot(a[i][k] - ahalf, a[i][k] - ahalf) * self.dt / 2
                err[i][k] = max(xerr, verr)

            for k in range(self.n):
                a[i+1][k] = self.force(t[i], x[i+1], v[i+1], k) / self.mass

            t[i+1] = t[i] + self.dt

            maxerr = 0
            for k in range(self.n):
                maxerr = max(maxerr, err[i][k])

            # dt = 0.9 * sqrt(etol / maxerr) * dt;

        var t_end = time.perf_counter_ns()
        print("sim_er took", (t_end-t_start)/1e9, "seconds with dt =", self.dt, "and n_tsteps =", self.n_tsteps)
        return (x, v, a, t, err)