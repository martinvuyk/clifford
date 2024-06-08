from surface import Surface
from vec3d import Vec3D


@value
struct FluidVolume[T: DType]:
    """Represents a fluid body."""

    var surface: Surface[T]
    """The FluidVolume's surface."""
    var mass: Scalar[T]
    """The FluidVolume's mass."""
    var density: Scalar[T]
    """The FluidVolume's density."""
    var viscosity: Scalar[T]
    """The FluidVolume's viscosity."""
    var velocity: Vec3D[T]
    """The FluidVolume's velocity."""
