# Clifford
This is a hobby project trying to setup a Physics engine in Mojo using Quaternions and DualQuaternions as the geometric basis wherever possible because... why not.


### Structures
- Vec3D: Can represent the location and surface normal, can represent a ray of light.
    - Has a rotation Quaternion.
    - Has a translation Quaternion.
    - Has an RGBA component.
- Surface: Represents a 3 dimensional surface.
    - Has a list of Vec3D that compose it.
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
    - velocity: The FluidVolume's velocity.