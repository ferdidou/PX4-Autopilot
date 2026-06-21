# hexa_tilting Gazebo-classic model

Six-arm omnidirectional tilting hexarotor (one rotor + one tilt servo per arm) for
the differential-allocation airframe `10020_gazebo-classic_hexa_tilting`
(`SYS_AUTOSTART=10020`, `CA_AIRFRAME=16`, `CA_METHOD=3`).

PX4 keeps Gazebo-classic models in the `sitl_gazebo-classic` submodule, so this copy
is kept here to live in the same repository as the firmware. To run the simulation,
make the model visible to Gazebo by linking it into the submodule's `models` folder:

```sh
ln -s ../../hexa_tilting_model/hexa_tilting \
  Tools/simulation/gazebo-classic/sitl_gazebo-classic/models/hexa_tilting
```

Then:

```sh
make px4_sitl gazebo-classic_hexa_tilting
```
