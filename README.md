# Double-pendulum
This code produces data for a double rod pendulum of variable length and weight distribution in an attempt to study its chaotic behavior. The differential equations are solved using Ruge-Kutta of 4th order. 

Let (x1,y1) and (x2,y2) be the positions of the center of mass of rod 1 and 2 respectively. We define α to be the ratio of the distance of the center of mass (measured from the origin) with the length l1 of the rod. β is defined in the same way as measured from the junction point and with the rod length l2. Thus, for a constant mass distribution, α=β=1/2. 

We can describe the center of mass’s position the following way:                                                       
<img align="center" src="https://github.com/guillaumecote/Double-pendulum/blob/master/math/eq1.png" alt="COM">

From this, some substituting and differentiating later, we find: 
<img align="center" src="https://github.com/guillaumecote/Double-pendulum/blob/master/math/eq2.png" alt="Lagrangian">

Where L is the Lagrangian, Ii the moment of inertia around the center of mass and V0 is a constant potential adjustment. V0 disappears in the equations of motions and only has an impact in the energy calculation.  We find the momenta p1 and p2 by differentiating the Lagrangian with respect to the generalized velocities, from which we find dtheta1/dt and dtheta2/dt. 
By using the Hamiltonian one can then find p1/dt and p2/dt. The exact equations are in the program, in the system function. 

The following Poincare sections visualizations were done with for m1=m2=1 and l1=l2=1, from theta1=0, theta2=PI/2 to theta1=theta2=PI.
<img align="center" src="https://github.com/guillaumecote/Double-pendulum/blob/master/sample%20output/img1.png">
<img align="center" src="https://github.com/guillaumecote/Double-pendulum/blob/master/sample%20output/img2.png">
<img align="center" src="https://github.com/guillaumecote/Double-pendulum/blob/master/sample%20output/img3.png">
<img align="center" src="https://github.com/guillaumecote/Double-pendulum/blob/master/sample%20output/img4.png">
<img align="center" src="https://github.com/guillaumecote/Double-pendulum/blob/master/sample%20output/img5.png">

Here are two high energy cases using theta1=theta2=PI and ptheta1=ptheta2=5 (left) and ptheta1=ptheta2=10 (right).
<img align="center" src="https://github.com/guillaumecote/Double-pendulum/blob/master/sample%20output/img6.png">

We see from the above that low energies and high energies are similar in the sense that they recover their predictability and resemble the behavior of a simple pendulum. We can also deduce that energies around E=15 are highly chaotic since the plot is uniform and no disctinct pattern can be seen.Around those energies, a small initial change in angle can lead to drastic differences in the long run, as shown in the following figure:

<img align="center" src="https://github.com/guillaumecote/Double-pendulum/blob/master/sample%20output/img7.png">

Please feel free to play around and study cases where the relative masses of the rods changes. Here's a few edge cases: 

<img align="center" src="https://github.com/guillaumecote/Double-pendulum/blob/master/sample%20output/img8.png">
