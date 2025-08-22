# Rebound and pykep examples

The rebound numerical propagator and the pykep library are useful tools for celestial mechanics and orbit transfers. Here are some simple demonstrations. See the software for more sophisticated examples.

## Smashing into a comet before it is too late

Here is an example of using the Lambert problem solver in pykep can be used to search for minimum delta v launch windows when in a pinch. Here there is an imminent impact by a comet at epoch 2200 days. We need to find a way to hit it with a rocket before this.
<code>
> conda install pykep
> python3 porkchop_planetary_defense.py
</code>

<img width="577" height="449" alt="Screenshot 2025-08-20 at 21 17 13" src="https://github.com/user-attachments/assets/9f54b41d-32c6-452c-b34f-d4eafca95dfc" />

Here is the lowest $\Delta v$ launch (there are quite many different options with similar requirement). In this case, the rocket makes a few orbits before impact, taking about 5 years from launch to impact at comet:

<img width="534" height="537" alt="Screenshot 2025-08-21 at 08 31 33" src="https://github.com/user-attachments/assets/36c390d3-0e43-4238-8cb3-56c9bc4391b9" />

## Earth-Mars Hohmann transfer

Here is a more standard example. Transferring from Earth to Mars. 

<code>
> python3 porkchop_earth_mars.py
</code>

<img width="1577" height="553" alt="mars" src="https://github.com/user-attachments/assets/d275c56c-ee10-4c22-9db3-0cb1fd2db050" />



<img width="1305" height="391" alt="Screenshot 2025-08-21 at 08 05 22" src="https://github.com/user-attachments/assets/7891c92e-ab63-4a61-8b2a-eb1d3e4fea00" />


Dario Izzo. (2019). esa/pykep: Bug fixes and more support on Equinoctial Elements (v2.3). Zenodo. https://doi.org/10.5281/zenodo.2575462
