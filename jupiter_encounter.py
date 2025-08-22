import numpy as np
import matplotlib.pyplot as plt
#import assist
import rebound
import pykep as pk
# -------------------------------------------------------------------
# Load planetary ephemerides (DE440 + satellites)
# -------------------------------------------------------------------
#https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440
#https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp
#ephem = assist.Ephem("data/linux_p1550p2650.440", "data/sb441-n16.bsp")
#print("Reference JD:", ephem.jd_ref)

# Create REBOUND simulation in AU, days, solar masses
sim = rebound.Simulation()
sim.units=("s","m","kg")
#sim.G = 6.6743e-11 # SI
sim.integrator = "ias15"        # high accuracy for close encounters
#sim.dt = 100

# Add Sun
t0=pk.epoch_from_string("2024-09-01 00:00:00").jd
print(t0)
epoch="JD%1.6f"%(t0)
sim.add("Sun",date=epoch)
sim.add("Earth",date=epoch)
sim.add("Jupiter",date=epoch)
sim.add("Mars",date=epoch)
sim.add("Venus",date=epoch)
sim.add("Mercury",date=epoch)
sim.add("Pluto",date=epoch)
n_part=len(sim.particles)
#earth = sim.particles["Earth"]
#print(earth)
#sc_offset = np.array([2*6.68459e-6, 0.0, 0.0])    # 0.01 AU away from Earth along +x
# normalized unit vector towards earth's velocity vector
#ve=n.array([earth.vx,earth.vy,earth.vz])
#ve0=ve=ve/np.linalg.norm(ve)
#sc_velocity = np.array([0.0, 0.025, 0.0]) # ~0.025 AU/day ~ 43 km/s relative

Noutputs = 10000
times = np.linspace(0, 50*365*24*3600, Noutputs)  # integrate 50 years

xyz = np.zeros((n_part,Noutputs, 3))
for i, t in enumerate(times):
    sim.integrate(times[i])
    for pi in range(n_part):
        sc = sim.particles[pi]
        xyz[pi,i,:] = [sc.x, sc.y, sc.z]
for pi in range(n_part):        
    plt.plot(xyz[pi,:,0], xyz[pi,:,1],label="stuff")
for pi in range(n_part):        
    plt.scatter(xyz[pi,0,0], xyz[pi,0,1],s=20)
print(xyz[1,0,:])
plt.show()
        
