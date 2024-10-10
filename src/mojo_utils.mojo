import math

alias Vec3 = SIMD[DType.float64, 4]

fn vec3(x: Float64, y: Float64, z: Float64) -> Vec3:
    return Vec3(x, y, z, 0)

fn length(v: Vec3) -> Float64:
    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

fn normalize(v: Vec3) -> Vec3:
    var l = length(v)
    return Vec3(v[0]/l, v[1]/l, v[2]/l)

fn dot(v1: Vec3, v2: Vec3) -> Float64:
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

fn cross(v1: Vec3, v2: Vec3) -> Vec3:
    return Vec3(
        v1[1]*v2[2] - v1[2]*v2[1],
        v1[2]*v2[0] - v1[0]*v2[2],
        v1[0]*v2[1] - v1[1]*v2[0]
    )

fn init_vector[T: CollectionElement](n: Int, t: T) -> List[T]:
    var v = List[T](capacity = n)

    for i in range(n):
        v.append(t)

    return v