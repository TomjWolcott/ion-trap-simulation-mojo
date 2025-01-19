from mojo_utils import *
from math import pi, sin, cos, sqrt, exp, floor, clamp, atan2
from python import Python, PythonObject
from potentials_3d import *
from utils import Variant
import time
import random
import os

alias Z = 1.0
alias e = 1.60217883e-19
alias eps0 = 8.854187817e-12
alias M_Yb = 2.8733965e-25
alias hbar: Float64 = 1.05457182e-34

fn main() raises:
    var sim = RingTrapSim(
        6,
        30.0 * 1e-6,
        2 * pi * 9e6,
        2 * pi * 10e6,
        n_tsteps=100000,
        dt=1e-10
    )

    sim.initialize_rotational_frequency(1e6)

    var result = sim.sim_er()
    result.show_animation(1000)
    result.save_to_csv("temp_ring_save.csv")
    var result2 = RingTrapResult.load_from_csv("temp_ring_save.csv")
    result2.show_animation(1000)

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

            for i in range(400000):
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
                    ys[k].append(math.sqrt(x[k][0]*x[k][0]+x[k][1]*x[k][1]))
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
        return self.potential1.force(t, x, v, k, self.mass) +
            self.potential2.force(t, x, v, k, self.mass) +
            self.potential3.force(t, x, v, k, self.mass) +
            self.potential4.force(t, x, v, k, self.mass) +
            self.potential5.force(t, x, v, k, self.mass)

    fn total_energy(self, t: Float64, x: List[Vec3], v: List[Vec3], k: Int) -> Float64:
        return self.potential1.energy(t, x, v, k, self.mass) +
            self.potential2.energy(t, x, v, k, self.mass) +
            self.potential3.energy(t, x, v, k, self.mass) +
            self.potential4.energy(t, x, v, k, self.mass) +
            self.potential5.energy(t, x, v, k, self.mass)

    fn sim_er(
        self
    ) -> GenericIonTrap3dResult:
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
        return GenericIonTrap3dResult(t, x, v, a, self.mass)

struct RingTrapSim[
    P3: Potential3D = NoPotential,
    P4: Potential3D = NoPotential,
    P5: Potential3D = NoPotential
]:
    var sim: IonTrapSim3D[RingPsuedoPotential, CoulombPotential, P3, P4, P5]

    fn __init__(
        inout self,
        num_ions: Int,
        radius: Float64,
        wr: Float64,
        wz: Float64,
        n_tsteps: Int = 10000,
        dt: Float64 = 1e-8,
        mass: Float64 = M_Yb
    ):
        var ring_potential = RingPsuedoPotential(radius, wr, wz)
        self.sim = IonTrapSim3D(
            num_ions,
            ring_potential,
            CoulombPotential(),
            P3(),
            P4(),
            P5(),
            n_tsteps=n_tsteps,
            dt=dt,
            mass=mass
        )

        _ = self.sim.set_equilibrium_positions(ring_potential.get_equilibrium_positions(self.sim.n))
        _ = self.sim.optimize_equilibrium(2e-16, True)

    fn initialize_rotational_frequency(inout self, rf_freq: Float64):
        for i in range(self.sim.n):
            var r = self.sim.x_0[i]
            r[2] = 0
            r = normalize(r)

            self.sim.v_0[i] = (rf_freq * (2 * pi * self.sim.potential1.radius)) * vec3(r[1], -r[0], 0)

    fn sim_er(self) -> RingTrapResult:
        return RingTrapResult(
            self.sim.sim_er(),
            self.sim.potential1.radius,
            self.sim.potential1.wr,
            self.sim.potential1.wz,
            Variant[NoneType, Float64](NoneType())
        )

trait IonTrap3dResult(Movable):
    fn save_to_csv(self, filename: String, skip: Int, bunch_len: Int) raises: ...

    @staticmethod
    fn load_from_csv(filename: String, chunk_len: Int) raises -> Self: ...

alias DEFAULT_SKIP = 10
alias DEFAULT_BUNCH_LEN = 1000
alias DEFAULT_CHUNK_LEN = 1000000

struct GenericIonTrap3dResult(IonTrap3dResult):
    var t: List[Float64]
    var x: List[List[Vec3]]
    var v: List[List[Vec3]]
    var a: List[List[Vec3]]
    var mass: Float64

    fn __init__(inout self, t: List[Float64], x: List[List[Vec3]], v: List[List[Vec3]], a: List[List[Vec3]], mass: Float64):
        self.t = t
        self.x = x
        self.v = v
        self.a = a
        self.mass = mass

    fn __moveinit__(inout self, owned existing: Self):
        self.t = existing.t
        self.x = existing.x
        self.v = existing.v
        self.a = existing.a
        self.mass = existing.mass

    # Empty header for the generic trap
    fn save_to_csv(self, filename: String, skip: Int = DEFAULT_SKIP, bunch_len: Int = DEFAULT_BUNCH_LEN) raises:
        self.save_to_csv_with_header(filename, "", skip, bunch_len)

    fn save_to_csv_with_header(self, filename: String, header: String, skip: Int, bunch_len: Int) raises:
        var f = open(filename, "rw")
        f.write(header)

        f.write("\nmass, {}\n".format(self.mass))
        f.write("Data:\n")

        for k1 in range(0, int(len(self.t)/skip/bunch_len)):
            var csv = String("")
            print("Done with line #", k1 * bunch_len)

            for k2 in range(0, bunch_len):
                var j = k1 * bunch_len + k2 - 1

                for i in range(len(self.x[0])*6+1):
                    var i1 = int((i-1)/6)
                    var i2 = int((i-1)/2) % 3
                    var axis = List("x", "y", "z")[i2]

                    if j == -1 and i == 0:
                        csv += String("Time (s)")
                    elif j == -1 and i % 2 == 0:
                        csv += String(", qubit #{} {}-position (m)").format(i1, axis)
                    elif j == -1:
                        csv += String(", qubit #{} {}-velocity (m/s)").format(i1, axis)
                    elif i == 0:
                        csv += String("\n {}").format(self.t[j*skip])
                    elif i % 2 == 0:
                        csv += String(", {}").format(self.x[j*skip][i1][i2])
                    else:
                        csv += String(", {}").format(self.v[j*skip][i1][i2])

            f.write(csv)

        f.close()

    @staticmethod
    fn load_parts_from_csv(filename: String, chunk_len: Int = DEFAULT_CHUNK_LEN) raises -> (String, List[Float64], List[List[Vec3]], List[List[Vec3]], List[List[Vec3]]):
        var f = open(filename, "r")
        var file_chunk = f.read(chunk_len)
        var t = List[Float64]()
        var x = List[List[Vec3]]()
        var v = List[List[Vec3]]()
        var a = List[List[Vec3]]()
        var mass = -1.0
        var num_lines = 0

        var data_start_index = file_chunk.find("Data:\n")

        # If "Data:" is in the csv, you can assume it has all other relevant info aswell
        if data_start_index == -1:
            raise Error("\"Data:\n\" is not found in the file")


        var header = file_chunk[:data_start_index]
        file_chunk = file_chunk[data_start_index + 6:]

        while file_chunk.byte_length() > 0:
            var lines = file_chunk.splitlines(False)

            # The last line is likely only partially captured, so remove it and read it next time
            if lines.size > 1:
                var last_line = lines.pop()
                _ = f.seek(-last_line.byte_length(), os.SEEK_CUR)

            num_lines += lines.size

            for line in lines:
                var cells = List[Float64]()

                for cell in line[].split(","):
                    try:
                        cells.append(string_to_float(cell[]))
                    except e:
                        continue

                # 7 is the minimum # of columns this csv could have
                if cells.size < 7:
                    continue

                if (cells.size - 1) % 6 != 0:
                    raise Error(String("Row has an incorrect number of cells: {}").format(cells.size))

                var num_ions = int((cells.size - 1) / 6)
                var xs = List[Vec3]()
                var vs = List[Vec3]()

                for i in range(num_ions):
                    xs.append(Vec3(cells[6 * i + 2], cells[6 * i + 4], cells[6 * i + 6]))
                    vs.append(Vec3(cells[6 * i + 1], cells[6 * i + 3], cells[6 * i + 5]))

                t.append(cells[0])
                x.append(xs)
                v.append(vs)
                a.append(List[Vec3]())

            file_chunk = f.read(chunk_len)

            print("Finished {} lines".format(num_lines))

        return (header, t, x, v, a)

    @staticmethod
    fn load_from_csv(filename: String, chunk_len: Int = DEFAULT_CHUNK_LEN) raises -> Self:
        var header: String
        var t: List[Float64]
        var x: List[List[Vec3]]
        var v: List[List[Vec3]]
        var a: List[List[Vec3]]
        (header, t, x, v, a) = GenericIonTrap3dResult.load_parts_from_csv(filename, chunk_len)

        var mass = find_data(header, "mass")

        return Self(t, x, v, a, mass)


    fn plot_kinetic_energy_vs_time(self) raises:
        var plt = Python.import_module("matplotlib.pyplot")
        var ts = Python.evaluate("[]")
        var energies = Python.evaluate("[]")

        for i in range(self.v.size):
            var energy = 0.0

            for j in range(self.v[i].size):
                energy += 0.5 * self.mass * dot(self.v[i][j], self.v[i][j])

            energies.append(6.242e+18 * energy)
            ts.append(self.t[i])

        plt.plot(ts, energies)
        plt.xlabel("Time (s)")
        plt.ylabel("Energy (eV)")
        plt.show()

struct RingTrapResult(IonTrap3dResult):
    var result: GenericIonTrap3dResult
    var radius: Float64
    var wr: Float64
    var wz: Float64
    # If NoneType, the result was made using the psuedopotential, otherwise it was with micromotion
    var rf_freq_opt: Variant[NoneType, Float64]

    fn __init__(inout self, owned generic_result: GenericIonTrap3dResult, radius: Float64, wr: Float64, wz: Float64, rf_freq_opt: Variant[NoneType, Float64]):
        self.result = generic_result^
        self.radius = radius
        self.wr = wr
        self.wz = wz
        self.rf_freq_opt = rf_freq_opt

    fn __moveinit__(inout self, owned existing: Self):
        self.result = existing.result^
        self.radius = existing.radius
        self.wr = existing.wr
        self.wz = existing.wz
        self.rf_freq_opt = existing.rf_freq_opt

    fn save_to_csv(self, filename: String, skip: Int = DEFAULT_SKIP, bunch_len: Int = DEFAULT_BUNCH_LEN) raises:
        var rf_freq_tag = String("")

        # Only write out the rf_freq if the result is from micromotion
        if self.rf_freq_opt.isa[Float64]():
            rf_freq_tag = "rf_freq,{}\n".format(self.rf_freq_opt.unsafe_get[Float64]())

        self.result.save_to_csv_with_header(
            filename,
            "ring-trap-result\nradius,{}\nwr,{}\nwz,{}\n{}".format(
                self.radius,
                self.wr,
                self.wz,
                rf_freq_tag
            ),
            skip,
            bunch_len
        )

    @staticmethod
    fn load_from_csv(filename: String, chunk_len: Int = DEFAULT_CHUNK_LEN) raises -> Self:
        var header: String
        var t: List[Float64]
        var x: List[List[Vec3]]
        var v: List[List[Vec3]]
        var a: List[List[Vec3]]
        (header, t, x, v, a) = GenericIonTrap3dResult.load_parts_from_csv(filename, chunk_len)
        print("doing other header stuff\n\nheader:\n", header, "\n\n")
        var mass = find_data(header, "mass")
        var radius = find_data(header, "radius")
        var wr = find_data(header, "wr")
        var wz = find_data(header, "wz")
        var rf_freq_opt: Variant[NoneType, Float64] = Variant[NoneType, Float64](NoneType())

        print("reading rf_freq?")
        try:
            var rf_freq = find_data(header, "rf_freq")
            rf_freq_opt = Variant[NoneType, Float64](rf_freq)
        except e:
            pass
        print("I readed it")

        return Self(GenericIonTrap3dResult(t, x, v, a, mass), radius, wr, wz, rf_freq_opt)

#     fn get_mode_energies(self) -> List[List[Float64]]:
#         var n = self.result.x[0].size
#         var n_tsteps = self.result.t.size
#         var efreqs: List[Float64]
#         var evecs: List[Float64]
#         efreqs, evecs = get_ring_eigenmodes(ensemble_properties, potentials)
#
#         mode_energies = np.zeros((n_tsteps, dims*n))
#         for i in range(n_tsteps):
#             for j in range(n_zero_freq_modes, dims*n):
#                 mode_energies[i][j] = (1/2) * ensemble_properties["mass"] * (efreqs[j] ** 2) * ( ((evecs[j].T @ r_sim[i]) ** 2) + (( (evecs[j].T @ v_sim[i]) / efreqs[j]) ** 2) )
#
#         return mode_energies

    fn show_animation(self, owned nt: Int) raises:
        var dt = self.result.t[1] - self.result.t[0]
        nt = min(len(self.result.x), nt)
        nqubits = len(self.result.x[0])
        print(nt, len(self.result.x), "|", nqubits, len(self.result.x[0]))
        var np = Python.import_module("numpy")
        var plt = Python.import_module("matplotlib.pyplot")
        var animation = Python.import_module("matplotlib.animation")
        var fig = plt.figure()
        var ax = fig.add_subplot(projection='3d')
        var artists = Python.evaluate("[]")
        ax.scatter(Python.evaluate("[0]"), Python.evaluate("[0]"), Python.evaluate("[1e-9]"), c="green")
        if self.radius > 0:
            var cs = Python.evaluate("[]")
            var ss = Python.evaluate("[]")
            var n = 200

            for i in range(n+1):
                cs.append(self.radius*cos(2*pi*i/n))
                ss.append(self.radius*sin(2*pi*i/n))
            ax.plot(cs, ss, 0)

        var ts = Python.evaluate("[]")
        var angle_vels = Python.evaluate("[]")
        var rads = Python.evaluate("[]")
        for i1 in range(nt-1):
            var i = i1 * int(floor((len(self.result.x) - 1) / (nt - 1)))
            var x = Python.evaluate("[]")
            var y = Python.evaluate("[]")
            var z = Python.evaluate("[]")
            var c = Python.evaluate("[]")
            var angle0 = atan2(self.result.x[i][0][1], self.result.x[i][0][0])
            var avr_freq = 0.0
            var avr_rad = 0.0

            for j1 in range(nqubits):
                var j = j1 * int(floor(len(self.result.x[i]) / nqubits))
                var angle = atan2(self.result.x[i][j][1], self.result.x[i][j][0])
                var ion1_radius = length(vec3(self.result.x[i][j][0], self.result.x[i][j][1], 0))
                var xNew = self.result.x[i][j][0]
                var yNew = self.result.x[i][j][1]
                print(xNew, yNew)

    #             xNew = ion1_radius * cos(angle - angle0)
    #             yNew = ion1_radius * sin(angle - angle0)
                x.append(xNew)
                y.append(yNew)
                z.append(self.result.x[i][j][2])
                if j1 == 0:
                    c.append("black")
                elif j1 == 1:
                    c.append("darkred")
                elif j1 % 5 == 0:
                    c.append("magenta")
                else:
                    c.append("red")

                var after_angle = atan2(self.result.x[i+1][j][1], self.result.x[i+1][j][0])
                var w = (after_angle-angle) / dt
                avr_freq += w / (2*pi*nqubits)
                avr_rad += self.radius / nqubits

            var t = i * dt
            if abs(avr_freq) < 1e8:
                ts.append(t)
                angle_vels.append(avr_freq)
                rads.append(avr_rad)

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

        print(anim)

# fn get_ring_eigenmodes() -> ():
#     eq_positions = get_ring_eq_pos(ensemble_properties, potentials=potentials)
#     hess = _hessian(eq_positions, ensemble_properties, potentials=potentials)
#
#     eigenvalues, eigenvectors = np.linalg.eigh(hess)
#     eigenfrequencies = np.sqrt(eigenvalues / ensemble_properties["mass"])
#     eigenvectors = eigenvectors.T
#     return eigenfrequencies, eigenvectors

fn string_to_float(str: String) raises -> Float64:
    if "e" in str:
        var parts = str.split("e")
        return atof(parts[0]) * pow(10.0, atof(parts[1]))
    else:
        return atof(str)

# Expects "{data_tag},{value}"
fn find_data(file_chunk: String, data_tag: String) raises -> Float64:
    var tag_index = file_chunk.find(data_tag + ",")
    if tag_index == -1:
        raise Error("Could not find \"{},\" tag in the file_chunk".format(data_tag))
    var value_start_index = tag_index + data_tag.byte_length() + 1

    var value_end_index: Int
    var next_comma = file_chunk[value_start_index:].find(",")
    var next_newline = file_chunk[value_start_index:].find("\n")

    if next_comma == -1 and next_newline == -1:
        value_end_index = file_chunk.byte_length() - 1
    elif next_newline == -1:
        value_end_index = value_start_index + next_comma
    elif next_newline < next_comma or next_comma == -1:
        value_end_index = value_start_index + next_newline
    else:
        value_end_index = value_start_index + next_comma

    return string_to_float(file_chunk[value_start_index:value_end_index])