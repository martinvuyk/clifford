from surface import Surface
from vec3d import Vec3D


@value
struct SoftVolume[T: DType]:
    """Represents a soft body."""

    var surface: Surface[T]
    """The SoftVolume's surface."""
    var mass: Scalar[T]
    """The SoftVolume's mass."""
    var density: Scalar[T]
    """The SoftVolume's density."""
    var velocity: Vec3D[T]
    """The SoftVolume's velocity."""
