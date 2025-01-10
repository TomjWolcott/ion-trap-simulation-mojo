from mojo_utils import *
from math import pi, sqrt, sin, cos

alias Z = 1.0
alias e = 1.60217883e-19
alias eps0 = 8.854187817e-12
alias M_Yb = 2.8733965e-25
alias hbar: Float64 = 1.05457182e-34

trait Potential3D(Copyable, Defaultable):
    fn force(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Vec3: ...

    fn energy(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Float64: ...

trait GetEquilibriumPositions():
    fn get_equilibrium_positions(self, num_ions: Int) -> List[Vec3]: ...

struct NoPotential(Potential3D, GetEquilibriumPositions):
    fn __init__(inout self):
        pass

    fn __copyinit__(inout self, other: Self):
        pass

    fn force(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Vec3:
        return vec3(0, 0, 0)

    fn energy(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Float64:
        return 0.0

    fn get_equilibrium_positions(self, num_ions: Int) -> List[Vec3]:
        return init_vector[Vec3](num_ions, vec3(0, 0, 0))

struct LaserPotential(Potential3D):
    var kvector: Vec3
    var origin: Vec3
    var laserWidth: Float64
    var s: Float64
    var decayRate: Float64
    var detuning: Float64

    fn __init__(inout self):
        self = Self(369 * 1e-9, vec3(0, 0, 0), 1e-6, 1, 19.6 * 1e6 * 2*pi, -19.6 * 1e6 * pi)

    fn __init__(inout self, wavelength: Float64, origin: Vec3, laserWidth: Float64, s: Float64, decayRate: Float64, detuning: Float64):
        var k: Float64 = (2*pi) / laserWavelength
        self.kvector = -vec3(k, k, k) / sqrt(3)
        self.origin = origin
        self.laserWidth = laserWidth
        self.s = s
        self.decayRate = decayRate
        self.detuning = detuning

    fn __copyinit__(inout self, other: Self):
        pass

    fn dist_to_laser(x: Vec3) -> Float64:
        var dist = r[k] - self.origin
        var kvector_length = sqrt(dot(kvector, kvector));
        var cross_product_length = length(cross(rminuso, kvector))

        return cross_product_length / kvector_length

    fn force(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Vec3:
        var kdotv = dot(self.kvector, v[k])
        var satParam = self.s * exp( -2 * pow(self.dist_to_laser(xs[k]), 2) / (self.laserWidth * self.laserWidth))
        var scattering_rate = (self.decayRate * satParam) / (1 + 2*satParam + pow( (2*(self.detuning - kdotv)) / self.decayRate, 2) )
        return self.kvector * (hbar * scattering_rate) / mass

    fn energy(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Float64:
        # TBD
        return 0.0

struct RingPsuedoPotential(Potential3D, GetEquilibriumPositions):
    var radius: Float64
    # Secular frequencies for radial and normal directions
    var wr: Float64
    var wz: Float64

    fn __init__(inout self):
        self = Self(40 * 1e-6 / (2*pi))

    fn __init__(inout self, radius: Float64, wr: Float64 = 9 * 2*pi*1e6, wz: Float64 = 10 * 2*pi*1e6):
        self.radius = radius
        self.wr = wr
        self.wz = wz

    fn __copyinit__(inout self, other: Self):
        self.radius = other.radius
        self.wr = other.wr
        self.wz = other.wz

    fn adapt_to_structure_constant(inout self, linear_particle_density: Float64, num_ions: Int, mass: Float64):
        self.wz = self.wr
        self.radius = num_ions / (2 * pi * linear_particle_density) * pow(3.0 * e*e / (8.0 * pi * eps0 * self.wr*self.wr * mass), 1.0 / 3.0)

    fn force(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Vec3:
        var r: Float64 = sqrt(xs[k][0]*xs[k][0] + xs[k][1]*xs[k][1])

        return -mass * vec3(
            (self.wr*self.wr) * (r - self.radius) * (xs[k][0] / r),
            (self.wr*self.wr) * (r - self.radius) * (xs[k][1] / r),
            (self.wz*self.wz) * xs[k][2]
        )

    fn energy(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Float64:
        var r: Float64 = sqrt(xs[k][0]*xs[k][0] + xs[k][1]*xs[k][1])

        return mass * (self.wr*self.wr) * pow(r - self.radius, 2) / 2 + (self.wz*self.wz) * pow(xs[k][2], 2) / 2

    fn get_equilibrium_positions(self, num_ions: Int) -> List[Vec3]:
        var xs = init_vector[Vec3](num_ions, vec3(0, 0, 0))

        for i in range(num_ions):
            var theta = 2*pi*i / num_ions

            xs[i] = self.radius * vec3(cos(theta), sin(theta), 0)

        return xs

    fn into_micro_motion_potential(self, owned rf_freq: Float64 = -1.0) -> RingMicroMotionPotential:
        if rf_freq == -1.0:
            rf_freq = 10 * self.wr

        return RingMicroMotionPotential(self.radius, rf_freq, self.wr, self.wz)

struct RingMicroMotionPotential(Potential3D):
    var radius: Float64
    var rf_freq: Float64
    # Secular frequencies for radial and normal directions
    var wr: Float64
    var wz: Float64

    fn __init__(inout self):
        self = Self(40 * 1e-6 / (2*pi))

    fn __init__(inout self, radius: Float64, rf_freq: Float64 = 100 * 2*pi*1e6, wr: Float64 = 9 * 2*pi*1e6, wz: Float64 = 10 * 2*pi*1e6):
        self.radius = radius
        self.rf_freq = rf_freq
        self.wr = wr
        self.wz = wz

    fn __copyinit__(inout self, other: Self):
        self.radius = other.radius
        self.rf_freq = other.rf_freq
        self.wr = other.wr
        self.wz = other.wz

    fn force(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Vec3:
        var r: Float64 = sqrt(xs[k][0]*xs[k][0] + xs[k][1]*xs[k][1])

        return -sqrt(2) * mass * self.rf_freq * cos(self.rf_freq * t) * vec3(
            self.wr * (r - self.radius) * (xs[k][0] / self.radius),
            self.wr * (r - self.radius) * (xs[k][1] / self.radius),
            -self.wz * xs[k][2]
        )

    fn energy(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Float64:
        var r: Float64 = sqrt(xs[k][0]*xs[k][0] + xs[k][1]*xs[k][1])

        return mass * self.rf_freq * sin(self.rf_freq * t) * (self.wr * pow(r - self.radius, 2) / 2 + self.wz * pow(xs[k][2], 2) / 2)

struct ElectricFieldPotential[field: fn(Float64, Vec3) -> Vec3](Potential3D):
    fn __init__(inout self):
        pass

    fn __copyinit__(inout self, other: Self):
        pass

    fn force(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Vec3:
        var field = self.field(t, xs[k])
#         print("        force:", field, xs[k])
        return field * e * Z

    fn energy(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Float64:
        var field = self.field(t, xs[k])
#         print("        energy:", field, xs[k])
        return dot(field, xs[k]) * e * Z


struct CoulombDefect(Potential3D):
    var qs: List[Float64]
    var points: List[Vec3]

    fn __init__(inout self):
        self.qs = List[Float64]()
        self.points = List[Vec3]()

    fn __init__(inout self, qs: List[Float64], points: List[Vec3]):
        self.qs = qs
        self.points = points

    fn __copyinit__(inout self, other: Self):
        self.qs = other.qs
        self.points = other.points

    fn force(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Vec3:
        var force = vec3(0, 0, 0)

        for i in range(len(self.qs)):
            var dr_mag = length(xs[k] - self.points[i])
            var dr = (xs[k] - self.points[i]) / dr_mag
            force += dr * ((Z * e * self.qs[i]) / (4*pi*eps0)) / (dr_mag * dr_mag)

        return force

    fn energy(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Float64:
        var energy = 0.0

        for i in range(len(self.qs)):
            var dr_mag = length(xs[k] - self.points[i])
            var dr = (xs[k] - self.points[i]) / dr_mag
            energy += (Z * e * self.qs[i]) / (8*pi*eps0 * dr_mag)

        return energy

struct CoulombPotential(Potential3D):
    fn __init__(inout self):
        pass

    fn __copyinit__(inout self, other: Self):
        pass

    fn force(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Vec3:
        var force = vec3(0, 0, 0)

        for i in range(len(xs)):
            if i != k:
                var dr_mag = length(xs[k] - xs[i])
                var dr = (xs[k] - xs[i]) / dr_mag
                force += dr * ((Z * Z * e * e) / (4*pi*eps0)) / (dr_mag * dr_mag)

        return force

    fn energy(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Float64:
        var energy = 0.0

        for i in range(len(xs)):
            if i != k:
                var dr_mag = length(xs[k] - xs[i])
                var dr = (xs[k] - xs[i]) / dr_mag
                energy += (Z * Z * e * e) / (8*pi*eps0 * dr_mag)

        return energy

struct WrappedPartialCoulombPotential(Potential3D):
    var num_ions: Int

    fn __init__(inout self):
        self = Self(100000000)

    fn __init__(inout self, num_ions: Int):
        self.num_ions = num_ions
        pass

    fn __copyinit__(inout self, other: Self):
        self.num_ions = num_ions
        pass

    fn force(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Vec3:
        var force = 0.0
        var num_ions = max(self.num_ions, len(xs))

        for i in range(-num_ions/2, num_ions/2+1):
            if i != 0:
                var j = i % len(xs)
                force += (1 if xs[j] < xs[k] else -1) * ((Z * Z * e * e) / (4*pi*eps0)) * (1 / ((xs[k] - xs[j]) * (xs[k] - xs[j])))

        return force

    fn energy(self, t: Float64, xs: List[Vec3], vs: List[Vec3], k: Int, mass: Float64) -> Float64:
        var energy = 0
        var num_ions = max(self.num_ions, len(xs))

        for i in range(-num_ions/2, num_ions/2+1):
            if i != 0:
                var j = (i+k) % len(xs)
                var dr_mag = length(xs[k] - xs[j])
                var dr = (xs[k] - xs[j]) / dr_mag
                energy += (Z * Z * e * e) / (8*pi*eps0 * dr_mag)

        return energy