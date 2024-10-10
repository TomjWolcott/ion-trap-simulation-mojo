from ion_trap_3d_lib import *
from mojo_utils import *
from collections.list import List
from python import Python, PythonObject
import benchmark
import time

fn main() raises:
    Python.add_to_path("src")
    ion_trap_3d = Python.import_module("ion_trap_3d")

#     dist_to_laser_test(10)
    sim_leapfrog_test(40)
#     sim_leapfrog_plot(3)
    sim_er_test(40)
    sim_er_plot(3)

fn inputs(num_qubits: Int) -> (Int, Float64, Float64, List[Vec3], List[Vec3]):
    r_0 = List[Vec3](capacity=num_qubits)
    v_0 = List[Vec3](capacity=num_qubits)

    for i in range(num_qubits):
        r_0.append(vec3(0, 0, 5e-6 * (i - (num_qubits - 1.0) / 2.0 + 1e-7)))
        v_0.append(vec3(0.0, 0.0, 0.0))

    # nt_steps, dt, etol, r_0, v_0
    return (400, 1e-8, 1e-9, r_0, v_0)

fn dist_to_laser_test(num_qubits: Int) raises:
    ion_trap_3d = Python.import_module("ion_trap_3d")
    var tup = inputs(num_qubits)
    var k = min(num_qubits - 1, 2)

    print("\n\ndist_to_laser\n---------------------")

    var mojo_time = time.perf_counter()
    var mojo_res = dist_to_laser(tup[3], k)
    mojo_time = time.perf_counter() - mojo_time

    print("Mojo: ", mojo_res, "", mojo_time, "s")

    var c_time = time.perf_counter()
    var c_res = ion_trap_3d.dist_to_laser(to_pylist(tup[3]), k)
    c_time = time.perf_counter() - c_time

    print("C: ", c_res, "", c_time, "s\nErr: ", mojo_res - c_res)

fn sim_leapfrog_test(num_qubits: Int) raises:
    ion_trap_3d = Python.import_module("ion_trap_3d")
    var tup = inputs(num_qubits)
    var T = tup[1] * tup[0]

    print("\n\nsim_leapfrog\n---------------------")

    var mojo_time = time.perf_counter()
    var mojo_res = sim_leapfrog(T, tup[1], tup[3], tup[4])
    mojo_time = time.perf_counter() - mojo_time
    var mojo_last = len(mojo_res[0])-1

    print("Mojo: ", mojo_res[0][mojo_last][1][2], "", mojo_time, "s")

    var c_r = to_pylist(tup[3])
    var c_v = to_pylist(tup[4])
    var c_time = time.perf_counter()
    var c_res = ion_trap_3d.sim_leapfrog(len(tup[3]), T, tup[1], c_r, c_v)
    c_time = time.perf_counter() - c_time
    var c_last = len(c_res[0])-1

    print("C: ", c_res[0][c_last][5], "", c_time, "s\nErr: ", mojo_res[0][mojo_last][1][2] - c_res[0][c_last][5])

fn sim_leapfrog_plot(num_qubits: Int) raises:
    ion_trap_3d = Python.import_module("ion_trap_3d")
    var tup = inputs(num_qubits)
    var T = tup[1] * tup[0]
    var np = Python.import_module("numpy")
    var plt = Python.import_module("matplotlib.pyplot")

    print("\n\nsim_leapfrog plot\n---------------------")

    var mojo_res = sim_leapfrog(T, tup[1], tup[3], tup[4])
    var c_res_ = ion_trap_3d.sim_leapfrog(len(tup[3]), T, tup[1], to_pylist(tup[3]), to_pylist(tup[4]))
    var c_res = (from_pylist(c_res_[0]), from_pylist(c_res_[1]))

#     var ax = plt.figure().add_subplot(projection='3d')

    for i in range(len(tup[3])):
        var mojo_xs = Python.evaluate("[]")
        var mojo_ys = Python.evaluate("[]")
        var mojo_zs = Python.evaluate("[]")
        var c_xs = Python.evaluate("[]")
        var c_ys = Python.evaluate("[]")
        var c_zs = Python.evaluate("[]")
        var ts = Python.evaluate("[]")

        for j in range(len(mojo_res[0])):
            mojo_xs.append(mojo_res[0][j][i][0])
            mojo_ys.append(mojo_res[0][j][i][1])
            mojo_zs.append(mojo_res[0][j][i][2])
            c_xs.append(c_res[0][j][i][0])
            c_ys.append(c_res[0][j][i][1])
            c_zs.append(c_res[0][j][i][2])
            ts.append(j*tup[1])

        plt.plot(ts, mojo_zs, label=str("q_m").__add__(str(i)))
        plt.plot(ts, c_zs, label=str("q_c").__add__(str(i)))
#         ax.plot(mojo_xs, mojo_ys, mojo_zs, label=str("q_m").__add__(str(i)))
#         ax.plot(c_xs, c_ys, c_zs, label=str("q_c").__add__(str(i)))

#     ax.legend()
    plt.legend()
    plt.show()

fn sim_er_test(num_qubits: Int) raises:
    ion_trap_3d = Python.import_module("ion_trap_3d")
    var tup = inputs(num_qubits)

    print("\n\nsim_er\n---------------------")

    var mojo_time = time.perf_counter()
    var mojo_res = sim_er(tup[0], tup[1], tup[2], tup[3], tup[4])
    mojo_time = time.perf_counter() - mojo_time
    var mojo_last = len(mojo_res[0])-1

    print("Mojo: ", mojo_res[0][mojo_last][1][2], "", mojo_time, "s")

    var c_r = to_pylist(tup[3])
    var c_v = to_pylist(tup[4])
    var c_time = time.perf_counter()
    var c_res = ion_trap_3d.sim_er(len(tup[3]), tup[0], tup[1], tup[2], c_r, c_v)
    c_time = time.perf_counter() - c_time
    var c_last = len(c_res[0])-1

    print("C: ", c_res[0][c_last][5], "", c_time, "s\nErr: ", mojo_res[0][mojo_last][1][2] - c_res[0][c_last][5])

fn sim_er_plot(num_qubits: Int) raises:
    ion_trap_3d = Python.import_module("ion_trap_3d")
    var tup = inputs(num_qubits)
    var np = Python.import_module("numpy")
    var plt = Python.import_module("matplotlib.pyplot")

    print("\n\nsim_er plot\n---------------------")

    var mojo_res = sim_er(tup[0], tup[1], tup[2], tup[3], tup[4])
    var c_res_ = ion_trap_3d.sim_er(len(tup[3]), tup[0], tup[1], tup[2], to_pylist(tup[3]), to_pylist(tup[4]))
    var c_res = (from_pylist(c_res_[0]), from_pylist(c_res_[1]))

    var ax = plt.figure().add_subplot(projection='3d')

    for i in range(len(tup[3])):
        var mojo_xs = Python.evaluate("[]")
        var mojo_ys = Python.evaluate("[]")
        var mojo_zs = Python.evaluate("[]")
        var c_xs = Python.evaluate("[]")
        var c_ys = Python.evaluate("[]")
        var c_zs = Python.evaluate("[]")
        var ts = Python.evaluate("[]")

        for j in range(len(mojo_res[0])):
            mojo_xs.append(mojo_res[0][j][i][0])
            mojo_ys.append(mojo_res[0][j][i][1])
            mojo_zs.append(mojo_res[0][j][i][2])
            c_xs.append(c_res[0][j][i][0])
            c_ys.append(c_res[0][j][i][1])
            c_zs.append(c_res[0][j][i][2])
            ts.append(j*tup[1])

#         plt.plot(ts, mojo_zs, label=str("q_m").__add__(str(i)))
#         plt.plot(ts, c_zs, label=str("q_c").__add__(str(i)))
        ax.plot(mojo_xs, mojo_ys, mojo_zs, label=str("q_m").__add__(str(i)))
        ax.plot(c_xs, c_ys, c_zs, label=str("q_c").__add__(str(i)))

    ax.legend()
#     plt.legend()
    plt.show()

fn to_pylist(lst: List[Vec3]) raises -> PythonObject:
    var py_lst = Python.evaluate("[]")

    for i in range(len(lst)):
        py_lst.append(lst[i][0])
        py_lst.append(lst[i][1])
        py_lst.append(lst[i][2])

    return py_lst

fn from_pylist(lst: PythonObject) raises -> List[List[Vec3]]:
    var nt: Int = len(lst)
    var n: Int = int(len(lst[0])/3)

    var mojo_lst = List[List[Vec3]](capacity = nt)

    for i in range(nt):
        mojo_lst.append(List[Vec3](capacity = n))
        for j in range(n):
            mojo_lst[i].append(vec3(lst[i][j*3].to_float64(), lst[i][j*3+1].to_float64(), lst[i][j*3+2].to_float64()))

    return mojo_lst