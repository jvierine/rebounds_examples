#!/usr/bin/env python3
"""
Demonstration of how solving the Lambert problem can be used to find a minimum delta v launch 
for a earth-mars transfer. This should approximately correspond to the Hohmann transfer window

Demonstrates how often there is a low delta v transfer window

Uses ESAs awesome pykep library. 

Juha Vierinen
"""
import numpy as np
import matplotlib.pyplot as plt
import pykep as pk
from pykep.orbit_plots import plot_planet, plot_lambert
import matplotlib.dates as mdates
date_form = mdates.DateFormatter('%Y-%m-%d')

mjd_epoch = np.datetime64('1858-11-17T00:00:00')

# use jpl ephemeris for earth and mars
earth=pk.planet.jpl_lp("earth")
mars=pk.planet.jpl_lp("mars")

# start search at this date (Modified Julian Date)
t_start=pk.epoch_from_string("2024-09-01 00:00:00").mjd

# 10 years into the future (days)
max_delay=6*365
max_dur=2*364

# search for optimal launch with minimum delta v requirement
# how many 
n_tof=500
departure_time = np.linspace(0,max_delay,num=n_tof)
travel_time = np.linspace(30,max_dur,num=n_tof)

mjd_np=[]
for i in range(n_tof):
    #delta_s=np.array((t_start+departure_time) * 24 * 3600,dtype=int)
    mjd_timedelta = np.timedelta64(int(departure_time[i]*24*3600+t_start*24*3600), 's')
    mjd_np.append(mjd_epoch + mjd_timedelta)
mjd_np=np.array(mjd_np)

# delta v required to exit and enter orbit. 
delta_v = np.zeros([n_tof,n_tof])
delta_v[:,:]=np.nan

best_l=None
best_dv=1e99
best_res={}
for i in range(n_tof):
    # when do we depart
    # launch from earth at t=0
#    print(departure_time[i]+t_start)
 #   print(pk.epoch(departure_time[i]+t_start,"mjd"))

    re,ve=earth.eph(pk.epoch(departure_time[i]+t_start, 'mjd'))

    for j in range(n_tof):
        arrival_time=travel_time[j]+departure_time[i]+t_start
        rc,vc=mars.eph(pk.epoch(arrival_time, 'mjd'))
        tof=travel_time[j]
        dt=tof*pk.DAY2SEC

        l = pk.lambert_problem(r1 = re, r2 = rc, tof = dt, mu = pk.MU_SUN, max_revs=5,cw=False)
        v1=l.get_v1()
        v2=l.get_v2()

        n_sol=len(v1)
        best_dv_this=1e99
        for si in range(n_sol):
            exit_dv=np.linalg.norm(np.array(ve)-np.array(v1[si]))
            entry_dv=np.linalg.norm(np.array(vc)-np.array(v2[si]))
            dv_this=exit_dv+entry_dv
            if dv_this < best_dv_this:
                best_dv_this=dv_this

            if dv_this < best_dv:
                best_l=l
                best_dv=dv_this
                best_res={"l":l,"t0":departure_time[i]+t_start,"arrival_time":arrival_time,"tof":dt,"re":re,"sol":si}
                print("found better delta v %1.2f (km/s) sol=%d travel_time %d days"%(best_dv/1e3,si,travel_time[j]))
        delta_v[i,j]=best_dv_this

# porkchop plot
fig,ax=plt.subplots()
im=ax.pcolormesh(mjd_np,travel_time,delta_v[:,:].T/1e3,cmap="turbo",vmin=0,vmax=50)
ax.set_title("Prograde orbits")
fig.colorbar(im,ax=ax,label=r"$\Delta v$ (km/s)")
#cb=plt.colorbar()
#cb.set_label()
ax.set_xlabel("Departure time (UTC)")
ax.set_ylabel("Travel time (days)")
ax.xaxis.set_major_formatter(date_form)
plt.xticks(rotation=20, ha='right')
plt.tight_layout()
plt.show()


# Create the figure and axis
fig = plt.figure(figsize = (16,5))
ax1 = fig.add_subplot(1, 3, 1, projection='3d')
ax1.scatter([0], [0], [0], color=['y'])

ax2 = fig.add_subplot(1, 3, 2, projection='3d')
ax2.scatter([0], [0], [0], color=['y'])
ax2.view_init(90, 0)

ax3 = fig.add_subplot(1, 3, 3, projection='3d')
ax3.scatter([0], [0], [0], color=['y'])
ax3.view_init(0,0)

for ax in [ax1, ax2, ax3]:
    # Plot the planet orbits
    plot_planet(earth, t0=pk.epoch(best_res["t0"],"mjd"),  legend=True, units=pk.AU, axes=ax)
    plot_planet(mars, t0=pk.epoch(best_res["arrival_time"],"mjd"),  legend=True, units=pk.AU, axes=ax)
    # Plot the Lambert solutions
    axis = plot_lambert(best_res["l"], sol=best_res["sol"], legend=True, units=pk.AU, axes=ax)

plt.show()
