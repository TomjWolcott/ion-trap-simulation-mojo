fn main():
    var abc = ABC[](5)

    abc.print_me()

struct NoXYZ(XYZ):
    fn xyz(self, x: Float64) -> Int:
        return 0

struct X(XYZ):
    fn xyz(self, x: Float64) -> Int:
        return 4 * int(x)

trait XYZ:
    fn xyz(self, x: Float64) -> Int: ...

struct ABC[x: XYZ = NoXYZ, y: XYZ = NoXYZ, z: XYZ = NoXYZ]:
    var a: Int

    fn __init__(inout self, a: Int):
        self.a = a

    fn print_me(self):
        print(self.a)