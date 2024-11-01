from forge_tools.complex import Quaternion

from surface import Surface
from vec3d import Vec3D


@value
struct FluidVolume[T: DType]:
    """Represents a fluid body.

    Notes:
        Methods assume incompressible flow and uniform pressure
        & density. Most use laminar flow. And viscocity is
        ignored when it would otherwise have a big performance hit
        to compute Navier-Stokes or an approximation.
    """

    var surface: Surface[T]
    """The FluidVolume's surface."""
    var mass: Scalar[T]
    """The FluidVolume's mass."""
    var density: Scalar[T]
    """The FluidVolume's density. It is assumed to be the same for
    the whole volume."""
    var viscosity: Scalar[T]
    """The FluidVolume's kinematic viscosity (v = μ / ρ). It is assumed to be the same for
    the whole volume."""
    var pressure: Scalar[T]
    """The FluidVolume's pressure. It is assumed to be the same for
    the whole volume."""
    var velocity: Vec3D[T]
    """The FluidVolume's velocity. It is assumed to be the same for
    the whole volume."""
    var divergence: Vec3D[T]
    """The FluidVolume's divergence. It is assumed to be the same for
    the whole volume."""

    fn __init__(
        inout self,
        surface: Surface[T],
        mass: Scalar[T],
        density: Scalar[T] = 1,
        viscosity: Scalar[T] = 0,
        pressure: Scalar[T] = 0,
        velocity: Vec3D[T] = Vec3D[T](0, 0, 0),
        divergence: Vec3D[T] = Vec3D[T](0, 0, 0),
    ):
        self.surface = surface
        self.mass = mass
        self.density = density
        self.viscosity = viscosity
        self.pressure = pressure
        self.velocity = velocity
        self.divergence = divergence

    fn acceleration(self, f: Vec3D[T]) -> Vec3D[T]:
        """Calculate the value of the derivative of the FluidVolume's
        velocity with given parameters.

        Notes:
            The Incompressible Navier-Stokes equation with constant
            density and viscosity is used. Which ends up becoming
            the Incompressible Euler equation with constant and
            uniform density + a viscosity term.

        Args:
            f: The external forces applied to the FluidVolume.

        Returns:
            The acceleration of the FluidVolume.
        """

        var divergence = self.divergence.qr.vec
        var v_div2_u = self.viscosity * divergence**2 * self.velocity.qt.vec
        var div_p_over_rho = divergence * self.pressure / self.density
        var u2 = self.velocity.qt.vec * divergence * self.velocity.qt.vec
        var q = Quaternion[T](f.qt.vec + v_div2_u - (div_p_over_rho + u2))
        return Vec3D[T](q.normalized(), q)

    # TODO: Hagen–Poiseuille equation for flow through pipe
    # fn flow(self):
    #     ...

    # TODO: Shallow water equations
    # fn recieve_impact(self) -> Self:
    #     ...
