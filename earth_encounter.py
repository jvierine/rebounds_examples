import numpy as np
import matplotlib.pyplot as plt
import assist
import rebound

# -------------------------------------------------------------------
# Load planetary ephemerides (DE440 + satellites)
# -------------------------------------------------------------------
#https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440
#https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp
ephem = assist.Ephem("data/linux_p1550p2650.440", "data/sb441-n16.bsp")
print("Reference JD:", ephem.jd_ref)

# Create REBOUND simulation in AU, days, solar masses
sim = rebound.Simulation()
sim.G = 0.00029591220828559104  # Gauss constant (AU^3 / M_sun / day^2)
sim.integrator = "ias15"        # high accuracy for close encounters
sim.dt = 0.000001                    # integration step in days

# Add Sun
sim.add(m=1.0, hash="SUN")

# Add planets with masses (not just test particles!)
planets = {
    "Earth": 3.00348959632e-6,   # M_sun
    "Venus": 2.4478383e-6,
    "Mars": 3.227151e-7,
    "Jupiter": 9.5458e-4
}
# You can add more planets if needed, but for clarity I use 4.

for name, m in planets.items():
    p = ephem.get_particle(name, 0)  # position/velocity at epoch
    sim.add(m=m, x=p.x, y=p.y, z=p.z, vx=p.vx, vy=p.vy, vz=p.vz, hash=name)

# -------------------------------------------------------------------
# Add a spacecraft (massless) set to encounter Earth
# -------------------------------------------------------------------
# Start near Earth but with a hyperbolic excess velocity relative to Earth
earth = sim.particles["Earth"]
print(earth)
sc_offset = np.array([2*6.68459e-6, 0.0, 0.0])    # 0.01 AU away from Earth along +x
sc_velocity = np.array([0.0, 0.025, 0.0]) # ~0.025 AU/day ~ 43 km/s relative
sim.add(m=0.0,
        x=earth.x+sc_offset[0], y=earth.y+sc_offset[1], z=earth.z+sc_offset[2],
        vx=earth.vx+sc_velocity[0], vy=earth.vy+sc_velocity[1], vz=earth.vz+sc_velocity[2],
        hash="SC")

print("Spacecraft added with initial offset and velocity relative to Earth.")

# -------------------------------------------------------------------
# Integrate forward and collect trajectories
# -------------------------------------------------------------------
Noutputs = 20000
times = np.linspace(-50, 50, Noutputs)  # integrate 1 year
xyz_sc = np.zeros((Noutputs, 3))
xyz_earth = np.zeros((Noutputs, 3))

for i, t in enumerate(times):
    sim.integrate(t)
    sc = sim.particles["SC"]
    earth = sim.particles["Earth"]
    xyz_sc[i] = [sc.x, sc.y, sc.z]
    xyz_earth[i] = [earth.x, earth.y, earth.z]

# -------------------------------------------------------------------
# Plot trajectory (XY plane)
# -------------------------------------------------------------------
plt.figure(figsize=(8,8))
plt.plot(xyz_sc[:,0], xyz_sc[:,1], ".",label="Spacecraft")
plt.plot(xyz_earth[:,0], xyz_earth[:,1], label="Earth")
plt.scatter(0,0,c="orange",marker="*",s=120,label="Sun")
plt.xlabel("x [AU]")
plt.ylabel("y [AU]")
plt.title("Near-Earth gravity assist trajectory")
plt.axis("equal")
plt.grid(ls="--",alpha=0.6)
plt.legend()
plt.show()