# Introduction to the Linked Cell Method
I created this repository to showcase an implementation of the Linked-Cell Method used for molecular dynamics (MD) simulations in which short range interactions dominate.
This implementation is based off of my reading of a very informative textbook written by Michael Griebel, Stephan Knapek and Gerhard Zumbusch entitled 
"Numerical Simulation in Molecular Dynamics". In this document I will summarize this implementation disucssing both the mathematical foundations as well as 
the algorithmic considerations of the program. I will also provide a simple output simulation. 

# Time Discretization - Transformation from a Continuos to Discrete Problem 
Time discretization is an essential first step in casting a problem in a format amenable to MD simulation. All MD simulations proceed through a series of iterative steps in which a
series of calculations are done to compute the forces, velocities and positions of the particles making up the system under study. After a complete cycle of calculations is complete, 
a variable storing the current timepoint is updated by a small increment to advance the simulation. Once the problem is changed from continuous to discrete, the differential operator is
replaced with the difference operator as shown below.  
<p align="center">
  <img src="https://github.com/viktorboris-1/1.-Linked-Cell-Method-for-Short-Range-Potentials/assets/93276956/ac0625bf-b84e-4256-aa57-c26baeb0e356"/>
</p>
Where the continous problem required the integration of continuos differential equations (In this case Newton's equations of motion) the descrtized problem requires the solution of a system of algebraic equations instead.
This simplification is part of what makes MD simulations a powerful tool for probing complex systems.

# Discretization of Newton's Equations of Motion
The particular equations for a particle's position and velocity solved for in this implementation of the Linked-Cell method are summarized below. These are derived from 
Newton's equations of motion which are then descritized and solved for the position and velocity at the next time point (n+1) for each particle i. The force is calculated using 
the Lennard-Jones potential. 
<p align="center">
  <img src="https://github.com/viktorboris-1/1.-Linked-Cell-Method-for-Short-Range-Potentials/assets/93276956/9c8330c1-45cc-4d42-aace-4099896ef845"/>
</p>

<p align="center">
  <img src="https://github.com/viktorboris-1/1.-Linked-Cell-Method-for-Short-Range-Potentials/assets/93276956/4b378858-c0e8-4558-8e53-0f3a9ed4564a"/>
</p>

# Integration Method of Stormer-Verlet
In this particular implementation of the so called velocity Stormer-Verlet, the current timepoint's position, force and velocity are used to calculate the next time point's 
position. With the updated positions the new force's are calculated in a pairwise fashion within a user specified cutoff radius. With the current positions and forces the 
new time point's velocities can be calculated.

<p align="center">
  <img src="https://github.com/viktorboris-1/1.-Linked-Cell-Method-for-Short-Range-Potentials/assets/93276956/d01f784f-74d6-4786-acee-bd0da8875d2f"/>
</p>
Below is a pseudo-code summary of the basic MD simulation implemented with the Linked-Cell method.

<p align="center">
  <img src="https://github.com/viktorboris-1/1.-Linked-Cell-Method-for-Short-Range-Potentials/assets/93276956/d5731a08-88ab-4ceb-b404-a09d83861a8f"/>
</p>

# The Cutoff Radius and the Linked Cell Method for Computing Short Range Forces
<p align="center">
  <img src="https://github.com/viktorboris-1/1.-Linked-Cell-Method-for-Short-Range-Potentials/assets/93276956/443f874c-7a8e-473d-8f1b-26b3314f6966"/>
</p>

# Example Simulation 
<p align="center">
  <img src="https://github.com/viktorboris-1/1.-Linked-Cell-Method-for-Short-Range-Potentials/assets/93276956/8d5a8b47-a67e-4588-91d1-822fb4573880"/>
</p>
