from mojo_utils import *
from math import pi, sin, cos, sqrt, exp, floor, clamp, atan2
from python import Python, PythonObject
from potentials_3d import *
import time
import random

alias Z = 1.0
alias e = 1.60217883e-19
alias eps0 = 8.854187817e-12
alias M_Yb = 2.8733965e-25
alias hbar: Float64 = 1.05457182e-34

struct VerticalLabel(Copyable, Movable, CollectionElement):
    var t: Float64
    var name: String

    fn __init__(inout self, t: Float64, name: String):
        self.t = t
        self.name = name

    fn __copyinit__(inout self, other: Self):
        self.t = other.t
        self.name = other.name

    fn __moveinit__(inout self, owned other: Self):
        self.t = other.t
        self.name = other.name

# Zigzag characterization
fn main2() raises:
    var np = Python.import_module("numpy")
    var plt = Python.import_module("matplotlib.pyplot")
    var cm = Python.import_module("matplotlib.cm")
    var fig = plt.figure()
    var ax = fig.add_subplot(projection='3d')

    var ns = np.linspace(10, 300, 10)
    var ds = np.linspace(1e-7, 1e-6, 30)
    var drs = np.zeros((len(ns), len(ds)))
    var xs = Python.evaluate("[]")
    var ys = Python.evaluate("[]")
    var zs = Python.evaluate("[]")

    var rs = Python.evaluate("[]")
#     var drs = Python.evaluate("[]")
    var colors = Python.evaluate("[]")
#
    for i in range(len(ns)):
        for j in range(len(ds)):
            var n = int(float(ns[i]))
            var radius = float(ds[j]) / (2*pi) * n

            var ion_trap_sim = IonTrapSim3D(
                n,
                RingPsuedoPotential(radius),
                CoulombPotential(),
                T=1e-5,
                dt=5e-9
            )
            ion_trap_sim.set_equilibrium_positions(ion_trap_sim.get_equilibrium_positions(ion_trap_sim.potential1))
            var converged = ion_trap_sim.optimize_equilibrium(1e-16, show_convergence=False)

            var d = 2*pi*radius / n
            var dr = abs(length(ion_trap_sim.x_0[0]) - length(ion_trap_sim.x_0[1]))

            dr = clamp(dr, 1e-19, 1e-5)
            print("i:", i, "n:", n, "radius:", radius, "convergence:", converged, "dr:", dr)

#             rs.append(radius)
#             drs.append(dr)
            xs.append(radius)
            ys.append()
            np.put(drs, i+len(ns)*j, dr)

            if converged:
                colors.append("red")
            else:
                colors.append("blue")

#     var ns = np.linspace(-2, 2, 30)
#     var ds = np.linspace(-2, 3, 40)
    print("hi")
    var tup = np.meshgrid(ns, ds)
#     print("hi again")
#     var qs = np.zeros((len(ns), len(ds)))
#
#     for i in range(len(ns)):
#         for j in range(len(ds)):
#             np.put(qs, i+len(ns)*j, sin(ns[i].to_float64())*ds[j].to_float64())
#             print("(",i,",",j,")", qs[i][j])

    ns, ds = tup[0], tup[1] # tuple (and its info) is dropped here unless you reference it below
    print("yo watup")
#     print(qs)
#     ax.set_zscale("log")
    ax.set_zlim(1e-19, 1e-5)
#     ax.plot_trisurf(ns.flatten(), ds.flatten(), drs.flatten(), linewidth = 1, antialiased = True, cmap=cm.Blues)
    ax.scatter(ns.flatten(), ds.flatten(), drs.flatten(), c=colors)
#     ax.plot(rs, drs)
    ax.set_xlabel("number of ions")
    ax.set_ylabel("ion-ion distance (m)")
    ax.set_zlabel("Adjacent ion absolute radius difference (m)")

    plt.title("Radial length difference vs ion-ion arc distance")
    plt.show()
    print(tup)

fn main4() raises:
    fn map_x0(i: Int, x: Vec3) -> Vec3:
        return vec3(x[0] + random.random_float64(-5e-7, 5e-7), x[1] + random.random_float64(-5e-7, 5e-7), random.random_float64(-1e-6, 1e-6))

    var radius = 50.0 * 2e-6 / (2*pi)
#     var defect_potential = CoulombDefect(List(0.5*e, 2*e), List(vec3(radius, 0, 5e-6), vec3(0, radius, 2e-6)))
    var ring_potential = RingPsuedoPotential(radius = radius)

    var ion_trap_sim = IonTrapSim3D(
        200,
        ring_potential,
        CoulombPotential(),
        T=1e-5,
        dt=1e-9
    )
    ion_trap_sim.set_equilibrium_positions(ion_trap_sim.get_equilibrium_positions(ion_trap_sim.potential1))
    ion_trap_sim.map_x_0[map_x0]()
    _ = ion_trap_sim.optimize_equilibrium(1e-16)

    var xs: List[List[Vec3]]
    var ts: List[Float64]
    xs,_,_,ts,_ = ion_trap_sim.sim_er()

    fn extra_points(t: Float64) -> List[Vec3]:
        return List(vec3(0, 0, 1e-6))

    var nt = 1000
    var nqubits = 4000
    animated_plot[extra_points](xs, ion_trap_sim.dt, nt, nqubits, List[VerticalLabel](), radius=ring_potential.radius)

# Animated ring simulations
fn main() raises:
    fn map_x0(i: Int, x: Vec3) -> Vec3:
        return vec3(x[0] + random.random_float64(-5e-8, 5e-8), x[1] + random.random_float64(-5e-8, 5e-8), random.random_float64(-1e-7, 1e-7))

    alias tilted_freq = 2*pi*2e5
    alias t_acc = 2e-5
    alias t_mag1 = 5e-5
    alias t_mag2 = 7e-5
    alias magnitude = 3e7

    fn edirection(t: Float64) -> Vec3:
        var freq = (t/t_acc) * tilted_freq / 2 if t < t_acc else tilted_freq
        var angle = freq * t + 0.5123123

        return vec3(cos(angle), sin(angle), 0)

    fn emagnitude(t: Float64) -> Float64:
        return magnitude if t < t_mag1 else (magnitude * (t_mag2 - t) / (t_mag2 - t_mag1) if t < t_mag2 else 0)

    fn efield(t: Float64, x: Vec3) -> Vec3:
        var dir = edirection(t)
        var perp_dir = vec3(-dir[1], dir[0], 0)

        var E = emagnitude(t) * (dot(x, dir) * dir - dot(x, perp_dir) * perp_dir)
#         print("            efield:", E)
        return E

    var epotential = ElectricFieldPotential[efield]()

    var radius = 10.0 * 2e-6 / (2*pi)
#     var defect_potential = CoulombDefect(List(0.5*e, 2*e), List(vec3(radius, 0, 5e-6), vec3(0, radius, 2e-6)))
    var ring_potential = RingPsuedoPotential(radius = radius)

#     ring_potential.adapt_to_structure_constant(0.95/3, 30, M_Yb)

    fn tilted_efield(t: Float64, x: Vec3) -> Vec3:
        return vec3(10, 0, 0)

    var potential_tilt = ElectricFieldPotential[tilted_efield]()

    var ion_trap_sim_psuedo = IonTrapSim3D(
        6,
        ring_potential,
        CoulombPotential(),
        epotential,
        potential_tilt,
#         defect_potential,
        T=1e-4,
        dt=1e-11
    )
    ion_trap_sim_psuedo.set_equilibrium_positions(ion_trap_sim_psuedo.get_equilibrium_positions(ion_trap_sim_psuedo.potential1))
    ion_trap_sim_psuedo.map_x_0[map_x0]()
    _ = ion_trap_sim_psuedo.optimize_equilibrium(3e-16, show_convergence=False)

    var micromotion_potential = ion_trap_sim_psuedo.potential1#.into_micro_motion_potential()
    var ion_trap_sim = IonTrapSim3D(
        ion_trap_sim_psuedo^,
        micromotion_potential,
        CoulombPotential(),
        epotential,
        potential_tilt,
#         defect_potential,
    )

    var xs: List[List[Vec3]]
    var vs: List[List[Vec3]]
    var ts: List[Float64]
#     ion_trap_sim.v_0[0] = vec3(0, 50, 0)
    xs,vs,_,ts,_ = ion_trap_sim.sim_er()

    var np = Python.import_module("numpy")

#     var csv = String(List[UInt8](capacity=len(xs)*(len(xs[0]) * 6 + 1) * 51 + 300))
#     var skip = 1
#
#     print("starting printing, n_tsteps:", len(xs), "n_ions:", len(xs[0]), "capacity: ", (len(xs)/skip)*(len(xs[0]) * 3 + 1) * 21 + 300)
#     for j in range(-1, int(len(ts)/skip)):
#         if j % 1000 == 0:
#             print("done with row", j)
#         for i in range(len(xs[0])*6+1):
#             var i1 = int((i-1)/6)
#             var i2 = int((i-1)/2) % 3
#             var axis = List("x", "y", "z")[i2]
#
# #             if j >= 0:
# #                 print(j, i1, i2, axis, xs[j*skip][i1], vs[j*skip][i1])
#
#             if j == -1 and i == 0:
#                 csv += String("Time (s)")
#             elif j == -1 and i % 2 == 0:
#                 csv += String(", qubit #{} {}-position (m)").format(i1, axis)
#             elif j == -1:
#                 csv += String(", qubit #{} {}-velocity (m/s)").format(i1, axis)
#             elif i == 0:
#                 csv += String("\n {}").format(ts[j*skip])
#             elif i % 2 == 0:
#                 csv += String(", {}").format(xs[j*skip][i1][i2])
#             else:
#                 csv += String(", {}").format(vs[j*skip][i1][i2])
#
#     print("done w/ csv")
#     var f = open("ring_trap_sim.csv", "rw")
#
#     f.write(csv)
#     f.close()
#     print("written!!")

    fn extra_points(t: Float64) -> List[Vec3]:
        return List(vec3(0, 0, 1e-6), 4e-6 * edirection(t))

    var nt = 2000
    var nqubits = 4000
    animated_plot[extra_points](
        xs,
        ion_trap_sim.dt,
        nt,
        nqubits,
        List(
            VerticalLabel(t_acc, "Fully accelerated"),
            VerticalLabel(t_mag1, "Pringle starts decreasing"),
            VerticalLabel(t_mag2, "Pringle gone")
        ),
        radius=ring_potential.radius
    )

fn to_pylist(xs: List[List[Vec3]]) raises -> PythonObject:
    var pylist = Python.evaluate("[]")

    for i in range(len(xs[0])):
        var x = Python.evaluate("[]")
        var y = Python.evaluate("[]")
        var z = Python.evaluate("[]")

        for j in range(len(xs)):
            x.append(xs[j][i][0])
            y.append(xs[j][i][1])
            z.append(xs[j][i][2])

        pylist.append(x)
        pylist.append(y)
        pylist.append(z)
    print("len(pylist) for List[List[Vec3]]", len(pylist), "inner", len(pylist[0]))
    return pylist

fn to_pylist(xs: List[Float64]) raises -> PythonObject:
    var pylist = Python.evaluate("[]")

    for i in range(len(xs)):
        pylist.append(xs[i])

    print("len(pylist) for List[Float64]", len(pylist))
    return pylist

fn animated_plot[extra_points: fn(Float64) -> List[Vec3]](
    points: List[List[Vec3]],
    inout dt: Float64,
    inout nt: Int,
    inout nqubits: Int,
    vertical_lines: List[VerticalLabel],
    radius: Float64 = 0.0
):
    nt = min(len(points), nt)
    nqubits = min(len(points[0]), nqubits)
    print(nt, len(points), "|", nqubits, len(points[0]))
    try:
        var np = Python.import_module("numpy")
        var plt = Python.import_module("matplotlib.pyplot")
        var animation = Python.import_module("matplotlib.animation")
        var fig = plt.figure()
        var ax = fig.add_subplot(projection='3d')
        var artists = Python.evaluate("[]")
#         if len(extra_points) > 0:
#             var x2 = Python.evaluate("[]")
#             var y2 = Python.evaluate("[]")
#             var z2 = Python.evaluate("[]")
#
#             for i in range(len(extra_points)):
#                 x2.append(extra_points[i][0])
#                 y2.append(extra_points[i][1])
#                 z2.append(extra_points[i][2])
#
#             print("1qq", x2, y2, z2)
#             ax.scatter(x2, y2, z2, c="blue", label="defects")
#             print("2qq", x2, y2, z2)
        ax.scatter(Python.evaluate("[0]"), Python.evaluate("[0]"), Python.evaluate("[1e-9]"), c="green")
        if radius > 0:
            var cs = Python.evaluate("[]")
            var ss = Python.evaluate("[]")
            var n = 200

            for i in range(n+1):
                cs.append(radius*cos(2*pi*i/n))
                ss.append(radius*sin(2*pi*i/n))
            ax.plot(cs, ss, 0)

        var ts = Python.evaluate("[]")
        var angle_vels = Python.evaluate("[]")
        var rads = Python.evaluate("[]")
        for i1 in range(nt-1):
            var i = i1 * int(floor((len(points) - 1) / (nt - 1)))
            var x = Python.evaluate("[]")
            var y = Python.evaluate("[]")
            var z = Python.evaluate("[]")
            var c = Python.evaluate("[]")
            var angle0 = atan2(points[i][0][1], points[i][0][0])
            var avr_freq = 0.0
            var avr_rad = 0.0

            for j1 in range(nqubits):
                var j = j1 * int(floor(len(points[i]) / nqubits))
                var angle = atan2(points[i][j][1], points[i][j][0])
                var radius = length(vec3(points[i][j][0], points[i][j][1], 0))
                var xNew = points[i][j][0]
                var yNew = points[i][j][1]

#                 xNew = radius * cos(angle - angle0)
#                 yNew = radius * sin(angle - angle0)
                x.append(xNew)
                y.append(yNew)
                z.append(points[i][j][2])
                if j1 == 0:
                    c.append("black")
                elif j1 == 1:
                    c.append("darkred")
                elif j1 % 5 == 0:
                    c.append("magenta")
                else:
                    c.append("red")

                var after_angle = atan2(points[i+1][j][1], points[i+1][j][0])
                var w = (after_angle-angle) / dt
                avr_freq += w / (2*pi*nqubits)
                avr_rad += radius / nqubits

            var t = i * dt
            if abs(avr_freq) < 1e8:
                ts.append(t)
                angle_vels.append(avr_freq)
                rads.append(avr_rad)

            for point in extra_points(t):
                x.append(point[][0])
                y.append(point[][1])
                z.append(point[][2])
                c.append("blue")

            artists.append(Python.evaluate("""lambda x, plt, t,n,avr_freq: [
                x,
                plt.gcf().text(0.1,0.1, 't={:.04e}s, avr_freq={:.04e} Hz'.format(t, avr_freq), style='italic', bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10}),
                plt.gcf().text(0.1,0.9, '#ions={}, ion_0=black, ion_1=darkred, ion_5n=magenta'.format(n), style='italic', bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10}),
            ]""")(ax.scatter(x, y, z, c=c), plt, t,nqubits, avr_freq))


        var anim = animation.ArtistAnimation(
            fig,
            artists,
            50,
            1000,
            repeat=True
        )

        ax.set_aspect('equal', share=True, adjustable="datalim")

        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_zlabel("z (m)")
        plt.show()
#         plt.clf()
        print(anim)

        plt.plot(ts, angle_vels)
        plt.xlabel("time (s)")
        plt.ylabel("average rotational frequency (Hz)")
        plt.title("Average rotational frequency vs time")

        for label in vertical_lines:
            plt.axvline(label[].t, color='r', ls=':')
            plt.text(label[].t, 0, label[].name, color='r', ha='right', va='bottom', rotation=90)

        plt.show()

        plt.plot(ts, rads, label="Average radius")
        plt.plot(ts, np.ones(len(ts)) * radius, label="Trap radius")
        plt.xlabel("time (s)")
        plt.ylabel("Radius (m)")
        plt.legend()
        plt.title("Average radius vs time")

        for label in vertical_lines:
            plt.axvline(label[].t, color='r', ls=':')
            plt.text(label[].t, radius, label[].name, color='r', ha='right', va='bottom', rotation=90)

        plt.show()
    except e:
        print(e)

fn convert_list_to_python(xs: List[Vec3]) raises -> (PythonObject, PythonObject, PythonObject):
    var x = Python.evaluate("[]")
    var y = Python.evaluate("[]")
    var z = Python.evaluate("[]")

    for i in range(len(xs)):
        x.append(xs[i][0])
        y.append(xs[i][1])
        z.append(xs[i][2])

    return x,y,z

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

    fn get_equilibrium_positions[P: GetEquilibriumPositions](self, potential: P) -> List[Vec3]:
        return potential.get_equilibrium_positions(self.n)

    fn set_equilibrium_positions(inout self, owned positions: List[Vec3]):
        self.x_0 = positions
        self.v_0 = init_vector[Vec3](self.n, Vec3(0, 0, 0))

    fn optimize_equilibrium(inout self, owned rate: Float64, show_convergence: Bool = True) -> Bool:
        var converged = False
        var updating_rate = rate
        var prev_total_energy = 0.0
        try:
            var x = init_vector[Vec3](self.n, Vec3(0, 0, 0))

            var xs = Python.evaluate("[]")
            var ys = init_vector[PythonObject](self.n, Python.evaluate("[]"))
            var np = Python.import_module("numpy")
            var plt = Python.import_module("matplotlib.pyplot")

            for k in range(len(ys)):
                ys[k] = Python.evaluate("[]")

            for i in range(100000):
                xs.append(i)
                var max_acc = 0.0
                var total_energy = 0.0
                for k in range(self.n):
                    var acc = self.force(0, self.x_0, self.v_0, k) / self.mass
                    x[k] = acc * updating_rate + self.x_0[k]
                    max_acc = max(length(acc), max_acc)
                    total_energy += self.total_energy(0, x, self.v_0, k)
#                     ys[k].append(length(x[k]))
#                     ys[k].append(atan2(x[k][1], x[k][0]))
                    ys[k].append(x[k][2])
#                 if prev_total_energy > total_energy:
#                     updating_rate *= 1.1
#                 else:
#                     updating_rate /= 1.1

                self.x_0 = x
                print(max_acc, "<- max_acc, total_energy ->", total_energy, "  (", updating_rate, ")")
                if max_acc < 1e-3:
                    converged = True
                    break

                if max_acc > 1e20 or updating_rate < 1e-24:
                    break

                prev_total_energy = total_energy

            if show_convergence:
                for k in range(self.n):
                    plt.plot(xs, np.array(ys[k]))

            if show_convergence:
                plt.show()
        except e:
            print(e)
            pass

        return converged

    fn map_x_0[map: fn(Int, Vec3) -> Vec3](inout self):
        for k in range(self.n):
            self.x_0[k] = map(k, self.x_0[k])

    fn force(self, t: Float64, x: List[Vec3], v: List[Vec3], k: Int) -> Vec3:
#         print("    x[k]", x[k])
#         print("    p1.force", k, self.potential1.force(t, x, v, k, self.mass))
#         print("    p2.force", k, self.potential2.force(t, x, v, k, self.mass))
#         print("    p3.force", k, self.potential3.force(t, x, v, k, self.mass))
#         print("    p4.force", k, self.potential4.force(t, x, v, k, self.mass))
#         print("    p5.force", k, self.potential5.force(t, x, v, k, self.mass))

        return self.potential1.force(t, x, v, k, self.mass) +
            self.potential2.force(t, x, v, k, self.mass) +
            self.potential3.force(t, x, v, k, self.mass) +
            self.potential4.force(t, x, v, k, self.mass) +
            self.potential5.force(t, x, v, k, self.mass)

    fn total_energy(self, t: Float64, x: List[Vec3], v: List[Vec3], k: Int) -> Float64:
#         print("    x[k]", x[k])
#         print("    p1.energy", k, self.potential1.energy(t, x, v, k, self.mass))
#         print("    p2.energy", k, self.potential2.energy(t, x, v, k, self.mass))
#         print("    p3.energy", k, self.potential3.energy(t, x, v, k, self.mass))
#         print("    p4.energy", k, self.potential4.energy(t, x, v, k, self.mass))
#         print("    p5.energy", k, self.potential5.energy(t, x, v, k, self.mass))

        return self.potential1.energy(t, x, v, k, self.mass) +
            self.potential2.energy(t, x, v, k, self.mass) +
            self.potential3.energy(t, x, v, k, self.mass) +
            self.potential4.energy(t, x, v, k, self.mass) +
            self.potential5.energy(t, x, v, k, self.mass)

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

        for k in range(self.n):
            x[0][k] = self.x_0[k]
            v[0][k] = self.v_0[k]
            a[0][k] = self.force(0, self.x_0, self.v_0, k) / self.mass
            print(a[0][k])

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