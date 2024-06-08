from math import sqrt

from complex_2 import Quaternion
from vec3d import Vec3D, Ray


@value
struct Surface[T: DType]:
    var vecs: List[Vec3D[T]]
    """The vectors that describe the surface."""
    var diffuse_reflectivity: Scalar[T]
    """Value between 0 and 1 that represents the proportion of light
    that gets scattered by the Surface."""
    var refraction_index: Scalar[T]
    """The refraction index of the Surface."""

    fn __init__(
        inout self,
        vecs: List[Vec3D[T]],
        diffuse_reflectivity: Scalar[T] = 0,
        refraction_index: Scalar[T] = 1,
    ):
        debug_assert(len(vecs) > 2, "surfaces need more than 2 vertices")
        self.vecs = vecs
        self.diffuse_reflectivity = diffuse_reflectivity
        self.refraction_index = refraction_index

    fn raycast(
        self, ray: Ray[T], prev_refraction_idx: Scalar[T]
    ) -> (Ray[T], Optional[Ray[T]], Scalar[T]):
        """Calculate the reflected, refracted, and diffused portions.

        Args:
            ray: The vector.
            prev_refraction_idx: The refraction index of the previous
                material.

        Returns:
            The reflected Ray according to its path, the refraction Ray,
                and the diffuse reflectance intensity.
        """
        alias tu = Tuple[Ray[T], Optional[Ray[T]], Scalar[T]]
        var res = ray.vec.intersect(self.vecs)
        if not res:
            return tu(ray, None, 0)
        var inters = res.value()

        var rot = ray.vec.qr.vec - 2 * (inters.dot_n) * inters.vec.qr.vec
        var specular_reflection = Vec3D(
            Quaternion[T](rot), inters.vec.qt, ray.vec.color * inters.vec.color
        )
        var diffuse_reflectance = (
            self.diffuse_reflectivity * ray.light_intensity * inters.dot_n
        )
        var n_r = prev_refraction_idx / self.refraction_index
        var sin_out = (n_r**2) * (1.0 - inters.dot_n**2)
        if sin_out < 1:
            var cos_out = sqrt(1.0 - sin_out)
            var norm_prod = (n_r * -inters.dot_n - cos_out) * ray.vec.qt.vec
            var refr = (n_r * ray.vec.qt.vec + norm_prod)
            var refraction = Vec3D(
                Quaternion[T](refr),
                inters.vec.qt,
                ray.vec.color * inters.vec.color,
            )
            var delta_reflectance = ray.light_intensity - diffuse_reflectance
            var r0 = (
                (prev_refraction_idx - self.refraction_index)
                / (prev_refraction_idx + self.refraction_index)
            ) ** 2
            var r_phi = r0 + (1 - r0) * (1 + inters.dot_n) ** 5
            var refl_intensity = delta_reflectance * r_phi
            var refr_intensity = delta_reflectance - refl_intensity
            return tu(
                Ray(specular_reflection, refl_intensity),
                Ray(refraction, refr_intensity),
                diffuse_reflectance,
            )
        # Total Internal Reflection
        var refl_intensity = ray.light_intensity - diffuse_reflectance
        return tu(
            Ray(specular_reflection, refl_intensity), None, diffuse_reflectance
        )

    fn area(self) -> Scalar[T]:
        # TODO
        return 0
