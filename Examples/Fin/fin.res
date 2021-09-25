Element type................................    quadrilateral
Number of nodes per element.................                4
Number of elements..........................                9
Number of nodes.............................               20
Number of integrating points................                4
Number of dimensions........................                2
Calculation time step.......................       0.1000E+02
Time-integration weighting parameter........       0.5000E+00
Initial time................................ 2021-07-05 00:00
Final time.................................. 2021-07-05 00:10
Time units..................................                s
Earth's latitude............................       0.9999E+04
Azimuth of the reference axes...............       0.9999E+04
Number of different property types..........                1


Thermal properties
Material         kx             ky            rho             cp            htco            ab
       1     0.1500E+02     0.1500E+02     0.1000E+03     0.3000E+03     0.1500E+02     0.0000E+00


Node coordinates
    Node          x              y
       1     0.0000E+00     0.0000E+00
       2     0.3703E-01     0.0000E+00
       3     0.7407E-01     0.0000E+00
       4     0.1111E+00     0.0000E+00
       5     0.1481E+00     0.0000E+00
       6     0.1852E+00     0.0000E+00
       7     0.2222E+00     0.0000E+00
       8     0.2592E+00     0.0000E+00
       9     0.2963E+00     0.0000E+00
      10     0.3333E+00     0.0000E+00
      11     0.0000E+00    -0.8330E-01
      12     0.3703E-01    -0.8330E-01
      13     0.7407E-01    -0.8330E-01
      14     0.1111E+00    -0.8330E-01
      15     0.1481E+00    -0.8330E-01
      16     0.1852E+00    -0.8330E-01
      17     0.2222E+00    -0.8330E-01
      18     0.2592E+00    -0.8330E-01
      19     0.2963E+00    -0.8330E-01
      20     0.3333E+00    -0.8330E-01


Element connectivities
 Element      Nodes
       1         11       1       2      12
       2         12       2       3      13
       3         13       3       4      14
       4         14       4       5      15
       5         15       5       6      16
       6         16       6       7      17
       7         17       7       8      18
       8         18       8       9      19
       9         19       9      10      20


There is no reservoir temperatures boundary conditions


Other temperatures boundary conditions
Fixed temperature

It is applied to
    Node      Value
       1 0.1100E+04
      11 0.1100E+04


There is no radiation boundary conditions


Convection boundary conditions

Harmonic function
Yearly mean temperature                        =  0.1000E+03
Amplitude of annual temperature variation      =  0.0000E+00
Yearly temperature phase difference            =  0.0000E+00
Yearly mean amplitude                          =  0.0000E+00
Amplitude of the daily amplitude variation     =  0.2803E-44
Yearly phase difference of the daily amplitude =  0.0000E+00
Daily temperature phase differencee            =  0.0000E+00

It is applied to
   Element     iside
         1         4
         2         4
         3         4
         4         4
         5         4
         6         4
         7         4
         8         4
         9         4
         1         3
         2         3
         3         3
         4         3
         5         3
         6         3
         7         3
         8         3
         9         3
         9         2


Initial temperatures
Constant value =  0.0000E+00


The results are printed every    60 time steps


The number of required point results is =    1


There are        20  equations and the skyline storage is        138
