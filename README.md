# Examples in Power Systems Analysis

## Steady-State Power Flow Analysis
Implementation of Power flow solver with
1. Newton-Raphson
2. Backward-Forward sweep (Gauss-Seidal method) for distribution system (Radial Network)

39 bus system data from MATPOWER

References:
1. Shirmohammadi, Dariush, et al. "A compensation-based power flow method for weakly meshed distribution and transmission networks." IEEE Transactions on power systems 3.2 (1988): 753-762.
2. Zimmerman, Ray Daniel, Carlos Edmundo Murillo-Sánchez, and Robert John Thomas. "MATPOWER: Steady-state operations, planning, and analysis tools for power systems research and education." IEEE Transactions on power systems 26.1 (2011): 12-19.

## Single-Machine Infinite-Bus example
References:
1. Andersson, Göran. "Dynamics and control of electric power systems." Lecture notes (2012): 227-0528.
2. Sauer, Peter W., Mangalore A. Pai, and Joe H. Chow. Power System Dynamics and Stability: With Synchrophasor Measurement and Power System Toolbox. John Wiley & Sons, 2017.


## Boundary Controlling UEP
The code runs Transient Stability Assessment based on the following methods
1. Time domain simulation
2. Potential Energy Boundary Surface (PEBS)
3. Boundary Controlling u.e.p. Method (BCU)

The data is imported in PSAT format. The code closely follows the example in Sauer et.al. [1].

References
1. Sauer, Peter W., Mangalore A. Pai, and Joe H. Chow. Power System Dynamics and Stability: With Synchrophasor Measurement and Power System Toolbox. John Wiley & Sons, 2017.
2. Chiang, Hsiao-Dong. Direct methods for stability analysis of electric power systems: theoretical foundation, BCU methodologies, and applications. John Wiley & Sons, 2011.

## Automatic Generation Control
Simulink model and parameter setting by looking at the root locus.

Reference:
1. Andersson, Göran. "Dynamics and control of electric power systems." Lecture notes (2012): 227-0528. (Chapter 4)

