from math import sqrt

from complex_2 import Quaternion

from vec3d import Vec3D
from surface import Surface


@value
struct BendingBreakingImpact[T: DType]:
    var transform: List[List[Vec3D[T]]]
    """List of tuples of start, end state."""
    var break_vec: Optional[Vec3D[T]]
    """The vector of the place where it breaks."""
    var break_surface: Optional[Surface[T]]
    """The surface that will go through self when breaking."""
    var reciprocal_force: Scalar[T]
    """The force applied to the incoming object."""
    var end_velocity_other: Vec3D[T]
    """The end velocity of the other after impact."""


@value
struct InelasticImpact[T: DType]:
    var transform: List[List[Vec3D[T]]]
    """List of tuples of start, end state."""
    var reciprocal_force: Scalar[T]
    """The force applied to the incoming object."""
    var end_velocity_other: Vec3D[T]
    """The end velocity of the other after impact."""


@value
struct ElasticImpact[T: DType]:
    var transform: List[List[Vec3D[T]]]
    """List of tuples of start, end state."""
    var reciprocal_force: Scalar[T]
    """The force applied to the incoming object."""
    var end_velocity_other: Vec3D[T]
    """The end velocity of the other after impact."""


@value
struct RigidVolume[T: DType]:
    """Represents a rigid body."""

    var surface: Surface[T]
    """The RigidVolume's surface."""
    var mass: Scalar[T]
    """The RigidVolume's mass."""
    var density: Scalar[T]
    """The RigidVolume's density."""
    var elasticity: Scalar[T]
    """The RigidVolume's Young's modulus of elasticity."""
    var max_tension: Scalar[T]
    """The RigidVolume's maximum axial tension pressure."""
    var max_compression: Scalar[T]
    """The RigidVolume's maximum axial compression pressure."""
    var velocity: Vec3D[T]
    """The RigidVolume's velocity."""

    fn __init__(
        inout self,
        surface: Surface[T],
        mass: Scalar[T],
        density: Scalar[T],
        elasticity: Scalar[T] = 0,
        max_tension: Scalar[T] = 0,
        max_compression: Scalar[T] = 0,
        velocity: Vec3D[T] = Vec3D[T](0, 0, 0),
    ):
        self.surface = surface
        self.mass = mass
        self.density = density
        self.elasticity = elasticity
        self.max_tension = max_tension
        self.max_compression = max_compression
        self.velocity = velocity

    fn volume(self) -> Scalar[T]:
        return self.mass / self.density

    fn center_of_mass(self) -> Vec3D[T]:
        return Vec3D.center(self.surface.vecs)

    fn outgoing_vector(self, incoming: Vec3D[T]) -> Vec3D[T]:
        """Calculate the outgoing vector from an incoming one.

        Args:
            incoming: The incoming vector.

        Returns:
            The outgoing vector.
        """
        # TODO: calculate distance to other side
        return Vec3D[T](0, 0, 0)

    fn force_to_break(
        self, incoming: Vec3D[T], outgoing: Vec3D[T]
    ) -> Scalar[T]:
        """Calculate the force to break the RigidVolume coming from an incoming
        vector to an outgoing one.

        Args:
            incoming: The incoming vector.
            outgoing: The outgoing vector.

        Returns:
            The force needed to break.
        """
        # TODO: use density and self.max_compression to estimate force needed
        return 0

    fn get_bending_breaking_impact(
        self, other: Self
    ) -> BendingBreakingImpact[T]:
        """Calculate the impact from another RigidVolume into self.

        Notes:
            This function in particular will take a long time and
                many more tests to actually work.

        Args:
            other: The other RigidVolume.

        Returns:
            The resulting BendingBreakingImpact.
        """

        var velocity = Vec3D(
            other.velocity.qr - self.velocity.qr,
            other.velocity.qt.vec - self.velocity.qt.vec,
        )
        var com_o = other.center_of_mass()
        var res = com_o.intersect(self.surface.vecs)
        if not res:
            var items = List[List[Vec3D[T]]]()
            for i in self.surface.vecs:
                items.append(List(i[], i[]))
            return BendingBreakingImpact(
                items, None, None, Scalar[T](0), other.velocity
            )
        var inters = res.value()
        var incoming = inters.vec.intersect(other.surface.vecs).value().vec
        var incidence = Vec3D[T](incoming.qr, inters.vec.qt)
        var out_vec = self.outgoing_vector(incidence)
        var break_force = self.force_to_break(incidence, out_vec)
        var kinetic_energy = 0.5 * other.mass * velocity.qt.__abs__() ** 2
        var distance_out = incidence.distance(out_vec)
        var breaking_energy = break_force * distance_out
        velocity.qt /= velocity.qt.__abs__()
        if breaking_energy < kinetic_energy:
            var new_magn = sqrt(
                (kinetic_energy - breaking_energy) * 2 / other.mass
            )
            velocity.qt *= new_magn
            var items = List[List[Vec3D[T]]]()
            for i in self.surface.vecs:
                items.append(List(i[], i[]))
            return BendingBreakingImpact(
                items, incidence, other.surface, break_force, velocity
            )

        var items = List[List[Vec3D[T]]]()

        @parameter
        fn break_bend(vec: Vec3D[T]) -> BendingBreakingImpact[T]:
            var break_vec = vec
            break_vec.qr.vec *= -1
            var bending_distance = incidence.distance(break_vec)
            var bending_surface_break = Surface(
                incidence.horizontal(self.surface.vecs)
            )
            return BendingBreakingImpact[T](
                items,
                break_vec,
                bending_surface_break,
                kinetic_energy / bending_distance,
                Vec3D[T](0, 0, 0),
            )

        var axes = incidence.get_orthonormal_set(self.surface.vecs)
        """A coordinate system around the vectors orthogonal to the 
        incidence vector."""
        var vert = Vec3D(incidence.qr * axes[2].qr, incidence.qt).orthogonal(
            self.surface.vecs
        )
        var center_vert = Vec3D.center(vert)
        """Center of mass for the vertical slice."""
        var i_xx: Scalar[T] = 0
        """Surface moment of inertia xx."""
        var i_yy: Scalar[T] = 0
        """Surface moment of inertia yy."""
        var i_xy: Scalar[T] = 0
        """Surface moment of inertia xy."""
        for i in vert:
            var x = i[]
            x.qt.vec[2] = center_vert.qt.vec[2]
            var x_dist = center_vert.distance(x)
            i_xx += x_dist**2
            var y = i[]
            y.qt.vec[1] = center_vert.qt.vec[1]
            var y_dist = center_vert.distance(y)
            i_yy += y_dist**2
            i_xy -= x_dist * y_dist

        var force = Vec3D(
            incidence.qr, incidence.qr * (kinetic_energy / distance_out)
        )
        for i in range(len(self.surface.vecs)):
            # rod bending equations for coordinate system unequal to
            # the main inertia axes
            var vec = self.surface.vecs[i]
            var z_surface = Surface(
                Vec3D(axes[2].qr, vec.qt).orthogonal(self.surface.vecs)
            )
            var z_area = z_surface.area()
            var divisor = i_xx * i_yy - i_xy**2
            var divisors = SIMD[T, 4](1, divisor, divisor, z_area)
            var distances = (incidence.qt - vec.qt).vec.__abs__()
            distances[3] = 1
            var momentums = force.qt.vec * distances
            var mb_same = momentums * SIMD[T, 4](0, i_yy, -i_xx, 1)
            var mb_diff = momentums * SIMD[T, 4](0, i_xy, -i_xy, 0)
            var locations = vec.qt.vec - Vec3D.center(z_surface.vecs).qt.vec
            locations[3] = 1
            var total_vec = ((mb_same + mb_diff) / divisors) * locations
            var sigma_zz = (total_vec**2).reduce_add()
            if sigma_zz > self.max_tension or sigma_zz < -self.max_compression:
                for _ in range(i, len(self.surface.vecs)):
                    items.append(List(vec, vec))
                return break_bend(vec)
            else:
                var displacement = vec
                displacement.qt += displacement.qt * sigma_zz / self.elasticity
                items.append(List(vec, displacement))

        var bending_distance = 0
        var bent_force = 0
        var remainder_force = kinetic_energy * bending_distance - bent_force
        var new_magn = sqrt(
            (remainder_force / bending_distance) * 2 / other.mass
        )
        velocity.qt *= new_magn
        return BendingBreakingImpact(
            items, None, None, remainder_force, velocity
        )

    fn break_form(self, surface: Surface[T]) -> List[Self]:
        """Break the RigidVolume, possibly into many smaller pieces
        or make a hole through.

        Args:
            surface: The Surface whose front facing normals
                will slice through the Rigidvolume.

        Returns:
            A list of all new Rigidvolumes.
        """

        # TODO
        return List[Self](self)

    fn get_inelastic_impact(self, other: Self) -> InelasticImpact[T]:
        """Calculate the inelastic impact from another RigidVolume into self.

        Args:
            other: The other RigidVolume.

        Returns:
            The resulting InelasticImpact.
        """
        var items = List[List[Vec3D[T]]]()
        for i in self.surface.vecs:
            items.append(List(i[], i[]))
        # TODO
        return InelasticImpact(items, 0, other.velocity)

    fn get_elastic_impact(self, other: Self) -> ElasticImpact[T]:
        """Calculate the elastic impact from another RigidVolume into self.

        Args:
            other: The other RigidVolume.

        Returns:
            The resulting ElasticImpact.
        """
        var items = List[List[Vec3D[T]]]()
        for i in self.surface.vecs:
            items.append(List(i[], i[]))
        # TODO
        return ElasticImpact(items, 0, other.velocity)
