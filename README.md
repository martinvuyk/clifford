# Clifford
This is a hobby project trying to setup a Physics engine in Mojo using Quaternions and DualQuaternions as the geometric basis, because the computational advantage over homogeneous coordinates and many other processes in the graphics rendering pipeline can accelerate and simplify the calculations.

### Example advantages:
- Transformation between 3D world space and 2D view screen can be done using DualQuaternions much more efficiently. 
- Occlussion mapping can be done in the engine itself by using the included vector normals (`Vec3D().qr`) and translation coordinates (`Vec3D().qt`) for the depth buffer.
- Intersection tests can be accelerated by just using the dot product to the vector normals (`Vec3D().qr`) before doing the dot product to the vector translation coordinates (`Vec3D().qt`). Avalable using `Vec3D().intersect(values: List[Vec3D])`
- Ray tracing with specular, diffusion, and refraction is available using `Surface.raycast(ray: Ray)` where Ray is a Vec3D plus a light intensity component.
- Bounding volume calculation is as simple as doing 
```mojo
var vol = RigidVolume()
var center = vol.center_of_mass()
var radius = center.furthest(vol.surface.vecs)
```

### Structures
- Vec3D: Can represent the location and surface normal, can represent a ray of light.
    - Has a rotation Quaternion.
    - Has a translation Quaternion.
    - Has an RGBA component.
- Surface: Represents a 3 dimensional surface.
    - vecs: A list of Vec3D that compose it.
    - diffuse_reflectivity: The percentage of light diffused when coliding.
    - refraction_index: The refraction index of the Surface.
- RigidVolume: Represents a rigid body.
    - surface: The RigidVolume's surface.
    - mass: The RigidVolume's mass.
    - density: The RigidVolume's density.
    - elasticity: The RigidVolume's Young's modulus of elasticity.
    - max_tension: The RigidVolume's maximum axial tension pressure.
    - max_compression: The RigidVolume's maximum axial compression pressure.
    - velocity: The RigidVolume's velocity.
- SoftVolume: Represents a soft body.
    - surface: The SoftVolume's surface.
    - mass: The SoftVolume's mass.
    - density: The SoftVolume's density.
    - velocity: The SoftVolume's velocity.
- FluidVolume: Represents a fluid body.
    - surface: The FluidVolume's surface.
    - mass: The FluidVolume's mass.
    - density: The FluidVolume's density.
    - viscosity: The FluidVolume's viscosity.
    - pressure: The FluidVolume's pressure.
    - velocity: The FluidVolume's velocity.
    - divergence: The FluidVolume's divergence.
