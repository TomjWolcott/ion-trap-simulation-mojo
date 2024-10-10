from ion_trap_1d_lib import *
from mojo_utils import *
from collections.list import List
from python import Python, PythonObject
import benchmark
import time

fn main() raises:
    Python.add_to_path("src")
    ion_trap_1d = Python.import_module("ion_trap_1d")

    total_energy_test(10)
    neg_grad_energy_test(10)
    force1_test(10)
    simleapfrog_test(10)
    simleafrog_plot(5)
    sim_er_test(10)
    sim_er_plot(5)
    sim_leapfrog_dless_test(10)
    sim_leapfrog_dless_plot(5)

fn inputs(num_qubits: Int) -> (Int, Float64, Float64, Float64, Float64, List[Float64], List[Float64]):
    var xs = List[Float64](capacity=num_qubits)
    var vs = List[Float64](capacity=num_qubits)

    for i in range(num_qubits):
        # These are intentionally oddly spaced x values, I wanted to make sure that the c++ & mojo code match
        # even when the positions are spaced sort of randomly
        xs.append(4e-6 * ((i*i / 100.0) + i - (num_qubits - 1.0) / 2.0) + 1e-6)
        vs.append(0.0)

    # (n_tsteps, dt, etol, M, b, x_0, v_0)
    return (1000000, 1e-10, 1e-8, M_Yb, 3e-20, xs, vs)

fn total_energy_test(num_qubits: Int) raises:
    ion_trap_1d = Python.import_module("ion_trap_1d")
    var tup = inputs(num_qubits)

    print("\n\ntotal_energy\n---------------------")

    var mojo_time = time.perf_counter()
    var mojo_res = total_energy(tup[5])
    mojo_time = time.perf_counter() - mojo_time

    print("Mojo: ", mojo_res, "", mojo_time, "s")

    var c_time = time.perf_counter()
    var c_res = ion_trap_1d.total_energy(to_pylist(tup[5]))
    c_time = time.perf_counter() - c_time

    print("   C: ", c_res, "", c_time, "s\nErr: ", mojo_res - c_res)

fn neg_grad_energy_test(num_qubits: Int) raises:
    ion_trap_1d = Python.import_module("ion_trap_1d")
    var tup = inputs(num_qubits)
    var k = min(num_qubits - 1, 2)
    var dx = 0.04

    print("\n\nneg_grad_energy\n---------------------")

    var mojo_time = time.perf_counter()
    var mojo_res = neg_grad_energy(tup[5], k, dx)
    mojo_time = time.perf_counter() - mojo_time

    print("Mojo: ", mojo_res, "", mojo_time, "s")

    var c_time = time.perf_counter()
    var c_res = ion_trap_1d.neg_grad_energy(to_pylist(tup[5]), k, dx)
    c_time = time.perf_counter() - c_time

    print("C: ", c_res, "", c_time, "s\nErr: ", mojo_res - c_res)

fn force1_test(num_qubits: Int) raises:
    ion_trap_1d = Python.import_module("ion_trap_1d")
    var tup = inputs(num_qubits)
    var k = min(num_qubits - 1, 2)

    print("\n\nforce1\n---------------------")

    var mojo_time = time.perf_counter()
    var mojo_res = force(tup[5], k)
    mojo_time = time.perf_counter() - mojo_time

    print("Mojo: ", mojo_res, "", mojo_time, "s")

    var c_time = time.perf_counter()
    var c_res = ion_trap_1d.force(to_pylist(tup[5]), k)
    c_time = time.perf_counter() - c_time

    print("C: ", c_res, "", c_time, "s\nErr: ", mojo_res - c_res)

# The other force function is not public/available from the c++ code
# fn force2_test(num_qubits: Int) raises:
#     ion_trap_1d = Python.import_module("ion_trap_1d")
#     var tup = inputs(num_qubits)
#     var k = min(num_qubits - 1, 2)
#     var b = 0.1
#
#     print("\n\nforce2\n---------------------")
#
#     var mojo_time = time.perf_counter()
#     var mojo_res = force(tup[5], tup[6], k, b)
#     mojo_time = time.perf_counter() - mojo_time
#
#     print("Mojo: ", mojo_res, "", mojo_time, "s")
#
#     var c_time = time.perf_counter()
#     var c_res = ion_trap_1d.force(len(tup[5]), to_pylist(tup[5]), to_pylist(tup[6]), k, b)
#     c_time = time.perf_counter() - c_time
#
#     print("C: ", c_res, "", c_time, "s\nErr: ", mojo_res - c_res)

fn simleapfrog_test(num_qubits: Int) raises:
    ion_trap_1d = Python.import_module("ion_trap_1d")
    var tup = inputs(num_qubits)
    var T = tup[1] * tup[0]

    print("\n\nsim_leapfrog\n---------------------")

    var mojo_time = time.perf_counter()
    var mojo_res = sim_leapfrog(T, tup[1], tup[3], tup[5], tup[6])
    var mojo_last = len(mojo_res[0])-1
    mojo_time = time.perf_counter() - mojo_time

    print("Mojo: ", mojo_res[0][mojo_last][1], "", mojo_time, "s")

    var c_r = to_pylist(tup[5])
    var c_v = to_pylist(tup[6])
    var c_time = time.perf_counter()
    var c_res = ion_trap_1d.sim_leapfrog(len(tup[5]), T, tup[1], tup[3], c_r, c_v)
    var c_last = len(c_res[0])-1
    c_time = time.perf_counter() - c_time

    print("C: ", c_res[0][c_last][1], "", c_time, "s\nErr: ", mojo_res[0][mojo_last][1] - c_res[0][c_last][1])

fn simleafrog_plot(num_qubits: Int) raises:
    ion_trap_1d = Python.import_module("ion_trap_1d")
    var np = Python.import_module("numpy")
    var plt = Python.import_module("matplotlib.pyplot")
    var tup = inputs(num_qubits)
    var T = tup[1] * tup[0]

    print("\n\nsim_leapfrog_plot")

    var mojo_tup = sim_leapfrog(T, tup[1], tup[3], tup[5], tup[6])
    var c_tup = ion_trap_1d.sim_leapfrog(len(tup[5]), T, tup[1], tup[3], to_pylist(tup[5]), to_pylist(tup[6]))

    for i in range(len(tup[5])):
        var mojo_xs = Python.evaluate("[]")
        var c_xs = Python.evaluate("[]")
        var ts = Python.evaluate("[]")

        for j in range(len(mojo_tup[0])):
            mojo_xs.append(mojo_tup[0][j][i])
            c_xs.append(c_tup[0][j][i])
            ts.append(j*tup[1])

        plt.plot(ts, np.array(mojo_xs), label=str("q_m").__add__(str(i)))
        plt.plot(ts, np.array(c_xs), label=str("q_c").__add__(str(i)))

    plt.legend()
    plt.show()

fn sim_er_test(num_qubits: Int) raises:
    ion_trap_1d = Python.import_module("ion_trap_1d")
    var tup = inputs(num_qubits)

    print("\n\nsim_er\n---------------------")

    var mojo_time = time.perf_counter()
    var mojo_res = sim_er(tup[0], tup[1], tup[2], tup[3], tup[4], tup[5], tup[6])
    var mojo_last = len(mojo_res[0])-1
    mojo_time = time.perf_counter() - mojo_time

    print("Mojo: ", mojo_res[0][mojo_last][1], "", mojo_time, "s")

    var c_r = to_pylist(tup[5])
    var c_v = to_pylist(tup[6])
    var c_time = time.perf_counter()
    var c_res = ion_trap_1d.sim_er(len(tup[5]), tup[0], tup[1], tup[2], tup[3], tup[4], c_r, c_v)
    var c_last = len(c_res[0])-1
    c_time = time.perf_counter() - c_time

    print("C: ", c_res[0][c_last][1], "", c_time, "s\nErr: ", mojo_res[0][mojo_last][1] - c_res[0][c_last][1])

fn sim_er_plot(num_qubits: Int) raises:
    ion_trap_1d = Python.import_module("ion_trap_1d")
    var np = Python.import_module("numpy")
    var plt = Python.import_module("matplotlib.pyplot")
    var tup = inputs(num_qubits)

    print("\n\nsim_er_plot")

    var mojo_tup = sim_er(tup[0], tup[1], tup[2], tup[3], tup[4], tup[5], tup[6])
    var c_tup = ion_trap_1d.sim_er(len(tup[5]), tup[0], tup[1], tup[2], tup[3], tup[4], to_pylist(tup[5]), to_pylist(tup[6]))

    for i in range(len(tup[5])):
        var mojo_xs = Python.evaluate("[]")
        var c_xs = Python.evaluate("[]")
        var ts = Python.evaluate("[]")

        for j in range(len(mojo_tup[0])):
            mojo_xs.append(mojo_tup[0][j][i])
            c_xs.append(c_tup[0][j][i])
            ts.append(mojo_tup[3][j][0])

        plt.plot(ts, np.array(mojo_xs), label=str("q_m").__add__(str(i)))
        plt.plot(ts, np.array(c_xs), label=str("q_c").__add__(str(i)))

    plt.legend()
    plt.show()

fn sim_leapfrog_dless_test(num_qubits: Int) raises:
    ion_trap_1d = Python.import_module("ion_trap_1d")
    var tup = inputs(num_qubits)
    var T = tup[1] * tup[0]

    print("\n\nsim_leapfrog_dless\n---------------------")

    var mojo_time = time.perf_counter()
    var mojo_res = sim_leapfrog_dless(T, tup[1], tup[5], tup[6])
    var mojo_last = len(mojo_res[0])-1
    mojo_time = time.perf_counter() - mojo_time

    print("Mojo: ", mojo_res[0][mojo_last][1], "", mojo_time, "s")

    tup = inputs(num_qubits) # because sim_leapfrog_dless modifies its inputs

    var c_time = time.perf_counter()
    var c_res = ion_trap_1d.sim_leapfrog_dless(num_qubits, T, tup[1], to_pylist(tup[5]), to_pylist(tup[6]))
    var c_last = len(c_res[0])-1
    c_time = time.perf_counter() - c_time

    print("C: ", c_res[0][c_last][1], "", c_time, "s\nErr: ", mojo_res[0][mojo_last][1] - c_res[0][c_last][1])

fn sim_leapfrog_dless_plot(num_qubits: Int) raises:
    ion_trap_1d = Python.import_module("ion_trap_1d")
    var np = Python.import_module("numpy")
    var plt = Python.import_module("matplotlib.pyplot")
    var tup = inputs(num_qubits)
    var T = tup[1] * tup[0]

    print("\n\nsim_leapfrog_dless_plot")

    var mojo_tup = sim_leapfrog_dless(T, tup[1], tup[5], tup[6])
    tup = inputs(num_qubits) # because sim_leapfrog_dless modifies its inputs
    var c_tup = ion_trap_1d.sim_leapfrog_dless(num_qubits, T, tup[1], to_pylist(tup[5]), to_pylist(tup[6]))

    for i in range(len(tup[5])):
        var mojo_xs = Python.evaluate("[]")
        var c_xs = Python.evaluate("[]")
        var ts = Python.evaluate("[]")

        for j in range(len(mojo_tup[0])):
            mojo_xs.append(mojo_tup[0][j][i])
            c_xs.append(c_tup[0][j][i])
            ts.append(j*tup[1])

        plt.plot(ts, np.array(mojo_xs), label=str("q_m").__add__(str(i)))
        plt.plot(ts, np.array(c_xs), label=str("q_c").__add__(str(i)))

    plt.legend()
    plt.show()

fn to_pylist(xs: List[Float64]) raises -> PythonObject:
    var xs_py = Python.evaluate("[]")

    for i in range(len(xs)):
        xs_py.append(xs[i])

    return xs_py