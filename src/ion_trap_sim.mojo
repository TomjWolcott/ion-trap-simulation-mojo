from mojo_utils import *
from math import pi

#   This file does not yet compile because traits cannot yet have parameters. Once
# this feature is made available this would be a nice way of structuring things.

alias Z = 1.0
alias e = 1.60217883e-19
alias eps0 = 8.854187817e-12
alias M_Yb = 2.8733965e-25
alias nu = 0.25 * 2*pi*1e6

fn main():
    var trap_sim = IonTrapSim[1, PsuedoPotential1D](5)

    print(trap_sim.total_energy())

struct PsuedoPotential1D(Potential[1]):
    @staticmethod
    fn energy(xs: List[Float64]) -> Float64:
        var energy: Float64 = 0

        for i in range(len(xs)):
            energy += 0.5 * M_Yb * (nu*nu) * (xs[i]*xs[i])

trait Potential[dims: Int]:
    @staticmethod
    fn energy(xs: List[SIMD[DType.float64, dims]]) -> Float64: ...

struct IonTrapSim[dims: Int]:
    alias VecN = SIMD[DType.float64, dims]

    var xs: List[Self.VecN]
    var vs: List[Self.VecN]
    var n: Int
    var dt: Float64
    var n_tsteps: Int
    var etol: Float64

    fn __init__(
        inout self,
        num_qubits: Int,
        n_tsteps: Int = 10000,
        dt: Float64 = 1e-8
    ):
        self.xs = List[Self.VecN](capacity=num_qubits)
        self.vs = List[Self.VecN](capacity=num_qubits)
        self.n = num_qubits

        for i in range(num_qubits):
            self.xs.append((i - (num_qubits - 1.0) / 2.0) + 6e-7)
            self.vs.append(0)

        self.n_tsteps = n_tsteps
        self.dt = dt
        self.etol = 1e-9

    fn total_energy(self) -> Float64:
        var total_energy = self.potential.energy(self.xs)

        for i in range(self.n):
            for j in range(i+1, self.n-1):
                total_energy += ((Z*Z * e*e) / (4 * pi * eps0)) * (1 / abs(self.xs[i] - self.xs[j]))

        return total_energy