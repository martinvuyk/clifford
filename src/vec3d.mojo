"""Defines Vec3D."""


from math import sqrt

from complex_2.quaternion import Quaternion


@value
@register_passable("trivial")
struct Ray[T: DType = DType.float64]:
    var vec: Vec3D[T]
    var light_intensity: Scalar[T]
    """The Vec3D's light intensity."""


@value
@register_passable("trivial")
struct Intersection[T: DType = DType.float64]:
    var vec: Vec3D[T]
    var dot_n: Scalar[T]


@register_passable("trivial")
struct Vec3D[T: DType = DType.float64]:
    """A structure used to represent 3D vectors with Quaternion
    logic underneath.
    """

    alias _scalar_type = Scalar[T]
    var qr: Quaternion[T]
    """The Vec3D's rotation Quaternion."""
    var qt: Quaternion[T]
    """The Vec3D's translation Quaternion."""
    alias _color_type = SIMD[T, 4]
    var color: Self._color_type
    """The Vec3D's color, can be: R, G, B, A."""

    fn __init__(
        inout self,
        x: Scalar[T],
        y: Scalar[T],
        z: Scalar[T],
        color: Self._color_type = Self._color_type(1, 1, 1, 1),
    ):
        # TODO calculate rotation
        # TODO calculate translation
        self.qr = Quaternion[T]()
        self.qt = Quaternion[T]()
        self.color = color

    fn __init__(
        inout self,
        qr: Quaternion[T],
        qt: Quaternion[T],
        color: Self._color_type = Self._color_type(1, 1, 1, 1),
    ):
        self.qr = qr
        self.qt = qt
        self.color = color

    fn distance(self, other: Self) -> Scalar[T]:
        """Calculate the euclidian distance between the two vectors.

        Args:
            other: The other vector.

        Returns:
            The distance.
        """
        return sqrt(((self.qt.vec - other.qt.vec) ** 2).reduce_add())

    fn intersect(self, values: List[Self]) -> Optional[Intersection[T]]:
        """Calculate the first intersection between a Vec3D and a
        List[Vec3D].

        The dot product with the rotation vector will be negative
        when doing dot product with a vector going in it's same
        direction but opposite sense. Magnitude of the rotational
        quaternions will always be 1. So this is calculating the
        cosine, and filtering < 1  means it's 90 degrees or less.
        The items are then sorted based on the translation distance.

        Returns:
            An Optional Vec3D representing the reflexion angle
                and the location of the intersection.
        """
        # FIXME: This might cause artifacts for coords around
        # the origin.
        var magn_ray_x2 = self.qt.__abs__() * 2

        var closest_vec = values[0]
        var closest_vec_angle_dot: Scalar[T] = 0
        var closest_dist = Scalar[T].MAX_FINITE
        for i in values:
            var item = (self.qr.vec * i[].qr.vec).reduce_add()
            if item < 0:
                # make a parabola around the distance to the ray
                var dot = (self.qt.vec * i[].qt.vec).reduce_add()
                var dist = ((dot) / magn_ray_x2 - 1) ** 2
                if dist < closest_dist:
                    closest_dist = dist
                    closest_vec = i[]
                    closest_vec_angle_dot = item
        if closest_dist == Scalar[T].MAX_FINITE:
            return None
        return Intersection(closest_vec, closest_vec_angle_dot)

    fn orthogonal(self, values: List[Self]) -> List[Self]:
        """Get the vectors orthogonal to self.

        Args:
            values: The vectors.

        Returns:
            The orthogonal vectors.
        """
        var items = List[Self]()
        for i in values:
            var item = (self.qr.vec * i[].qr.vec).reduce_add()
            if item == 0:
                items.append(i[])
        return items

    fn vertical(self, values: List[Self]) -> List[Self]:
        """Get the vectors parallel to the z axis at
        the x and y direction of self.

        Args:
            values: The vectors.

        Returns:
            The orthogonal vectors.
        """

        var items = List[Self]()
        for i in values:
            if self.qt.i == i[].qt.i and self.qt.j == i[].qt.j:
                items.append(i[])
        return items

    fn horizontal(self, values: List[Self]) -> List[Self]:
        """Get the x,y plane at the z cordinate of self.

        Args:
            values: The vectors.

        Returns:
            The orthogonal vectors.
        """
        var items = List[Self]()
        for i in values:
            if self.qt.k == i[].qt.k:
                items.append(i[])
        return items

    @staticmethod
    fn center(values: List[Self]) -> Self:
        """Calculate the center of the vectors.

        Args:
            values: The vectors.

        Returns:
            The average coordinate of the vectors.
        """
        var avg = Vec3D[T](0, 0, 0)
        var amnt = 0
        for vec in values:
            amnt += 1
            avg.qr.vec += vec[].qr.vec
            avg.qt.vec += vec[].qt.vec
        avg.qr.vec /= amnt
        avg.qt.vec /= amnt
        return avg

    fn furthest(self, values: List[Self]) -> Self:
        """Get the furthest vector from self."""
        # TODO
        return self

    fn cross(self, other: Self) -> Self:
        """Get the cross product from self to other."""
        # TODO
        return Self(0, 0, 0)

    fn get_orthonormal_set(self, values: List[Self]) -> (Self, Self, Self):
        """Get 3 axes from the vector. The function looks for the
        two furthest away vectors to use as the z axis, uses the existing
        vector as the x axis and finds the y axis by cross product.

        Args:
            values: The vectors.

        Returns:
            The new coordinate system.
        """

        var ort = self.orthogonal(values)
        var center = Self.center(ort)
        var furthest_pos_vec = values[0]
        var furthest_pos_vec_dist: Scalar[T] = 0
        var furthest_neg_vec = values[0]
        var furthest_neg_vec_dist: Scalar[T] = 0

        for i in values:
            var item = (center.qr.vec * i[].qr.vec).reduce_add()
            # make a parabola around the distance to the center
            var dot = (center.qt.vec * i[].qt.vec).reduce_add()
            var dist = ((dot) / 2 - 1) ** 2
            if item > 0:
                if dist > furthest_pos_vec_dist:
                    furthest_pos_vec_dist = dist
                    furthest_pos_vec = i[]
            else:
                if dist > furthest_neg_vec_dist:
                    furthest_neg_vec_dist = dist
                    furthest_neg_vec = i[]

        var furthest_pos_out = furthest_pos_vec.furthest(values)
        var dist_pos = furthest_pos_vec.distance(furthest_pos_out)
        var furthest_neg_out = furthest_neg_vec.furthest(values)
        var dist_neg = furthest_neg_vec.distance(furthest_neg_out)

        var direction: Quaternion[T]
        if dist_pos > dist_neg:
            direction = furthest_pos_vec.qt - furthest_pos_out.qt
        else:
            direction = furthest_neg_vec.qt - furthest_neg_out.qt
        direction.normalize()

        var z = Vec3D(direction, center.qt)
        var x = Vec3D(self.qr, center.qt)
        var y = x.cross(z)
        return x, y, z
