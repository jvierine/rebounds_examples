#!/usr/bin/env python3
"""
Demonstration of how solving the Lambert problem can be used to find a minimum delta v launch 
so that a rocket laden with fusion bombs can be sent to hit a newly discovered comet that would collide with Earth.

Uses ESAs awesome pykep library. 

Juha Vierinen
"""
import numpy as np
import matplotlib.pyplot as plt
import pykep as pk
from pykep.orbit_plots import plot_planet, plot_lambert

# ------------------------ Configuration ------------------------ #

# Start epoch (MJD2000). You can set this to "now" if desired.
# For a simple run, set t0 = 0 (J2000). We'll define the impact at +365 d.
t0 = 0.0  # MJD2000
earth=pk.planet.jpl_lp("earth")

# this many days after epoch, there will be an impact
t_impact=2200
t = pk.epoch(t_impact, 'mjd2000')

# Get position and velocity
r, v = earth.eph(t)

# This is the velocity that the comet impacts the Earth
v_comet = np.array([0,-30e3,30e3])
# create a keplerian element for the comet
# position at earth and velocity the specified velocity
comet = pk.planet.keplerian(t,r, v_comet, pk.MU_SUN, 10, 10, 10,  'Comet')

# search for optimal launch with minimum delta v requirement
n_tof=1000
departure_time = np.linspace(0,t_impact,num=n_tof)
arrival_time = np.linspace(0,t_impact,num=n_tof)
delta_v = np.zeros([2,n_tof,n_tof])
delta_v[:,:]=np.nan
best_l=None
best_dv=1e99
best_res={}
for i in range(n_tof):
    # when do we depart
    # launch from earth at t=0
    re,ve=earth.eph(pk.epoch(departure_time[i], 'mjd2000'))

    for j in range(n_tof):
        # reach target 30 days before impact
        if (arrival_time[j]>(departure_time[i]+10)) and arrival_time[j]<(t_impact-30) :
            rc,vc=comet.eph(pk.epoch(arrival_time[j], 'mjd2000'))
            tof=arrival_time[j]-departure_time[i]
            dt=tof*pk.DAY2SEC
            l = pk.lambert_problem(r1 = re, r2 = rc, tof = dt, mu = pk.MU_SUN, max_revs=10,cw=False)
            v1=l.get_v1()
            n_sol=len(v1)
            best_dv_this=1e99
            for si in range(n_sol):
                exit_dv_cw=np.linalg.norm(np.array(ve)-np.array(v1[si]))
                if exit_dv_cw < best_dv_this:
                    best_dv_this=exit_dv_cw
                if exit_dv_cw < best_dv:
                    best_l=l
                    best_dv=exit_dv_cw
                    best_res={"l":l,"t0":departure_time[i],"t1":arrival_time[j],"tof":dt,"re":re,"sol":si}
                    print("found better delta v %1.2f (km/s) si=%d tof=%1.2f years"%(best_dv/1e3,si,dt/24/3600/365))
            delta_v[0,i,j]=best_dv_this

# porkchop plot
plt.pcolormesh(departure_time,arrival_time,delta_v[0,:,:].T/1e3,cmap="turbo",vmin=0,vmax=50)
plt.title("Prograde orbits")
cb=plt.colorbar()
cb.set_label(r"$\Delta v$ (km/s)")
plt.xlabel("Departure time (days)")
plt.ylabel("Arrival time (days)")
plt.show()

# plot one lambert solution
t_exit=best_res["t0"]
t_arrival=best_res["t1"]
# position of earth and comet at launch and impact
re,ve=earth.eph(pk.epoch(t_exit, 'mjd2000'))
rc,vc=comet.eph(pk.epoch(t_arrival, 'mjd2000'))

# solve the Lambert problem for the best solution
l = pk.lambert_problem(r1 = re, r2 = rc, tof = best_res["tof"], mu = pk.MU_SUN, max_revs=10,cw=False)
rocket_v=l.get_v1()[best_res["sol"]]
rocket_r=l.get_x()[best_res["sol"]]
rocket = pk.planet.keplerian(pk.epoch(t_exit, 'mjd2000'),re, rocket_v, pk.MU_SUN, 10, 10, 10,  'Rocket')

# plot orbits of Earth and the impacting comet until impact time
times = np.linspace(-100, t_impact, 10000)
rs=[]
crs=[]
for t in times:
    r,v=earth.eph(pk.epoch(t,"mjd2000"))
    rs.append(r)
    r,v=comet.eph(pk.epoch(t,"mjd2000"))
    crs.append(r)

# figure out where the rocket is between launch and impact times
times = np.linspace(t_exit, t_arrival, 1000)
rrocket=[]
for t in times:
    r,v=rocket.eph(pk.epoch(t,"mjd2000"))
    rrocket.append(r)
# convert to numpy arrays
rs = np.array(rs)
crs = np.array(crs)
rrocket = np.array(rrocket)
# plot positions
plt.figure(figsize=(6,6))
plt.plot(rs[:,0]/pk.AU, rs[:,1]/pk.AU, label="Earth orbit")
plt.plot(crs[:,0]/pk.AU, crs[:,1]/pk.AU, label="Comet")
plt.plot(rrocket[0,0]/pk.AU, rrocket[0,1]/pk.AU, "*", label="Rocket is lauched",zorder=5)
plt.plot(rrocket[:,0]/pk.AU, rrocket[:,1]/pk.AU, label="Rocket")
plt.plot(crs[-1,0]/pk.AU, crs[-1,1]/pk.AU, "*", label="Comet impacts Earth",zorder=5)
plt.plot(rrocket[-1,0]/pk.AU, rrocket[-1,1]/pk.AU, "*", label="Rocket impacts comet",zorder=5)
plt.scatter([0],[0],c='orange',marker='*',s=200,label="Sun")
plt.axis('equal')
plt.xlabel("x [AU]")
plt.ylabel("y [AU]")
plt.title("Planetary defense")
plt.legend()
plt.grid(True)
plt.show()

