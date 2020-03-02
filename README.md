# Double-pendulum
This code produces data for a double rod pendulum of variable length and weight distribution in an attempt to study its chaotic behavior. The differential equations are solved using Ruge-Kutta of 4th order. 

Let (x1,y1) and (x2,y2) be the positions of the center of mass of rod 1 and 2 respectively. We define α to be the ratio of the distance of the center of mass (measured from the origin) with the length l1 of the rod. β is defined in the same way as measured from the junction point and with the rod length l2. Thus, for a constant mass distribution, α=β=1/2. 

We can describe the center of mass’s position the following way:                                                       


From this, some substituting and differentiating later, we find: 
      
 
 
                     
 
 
       
 
 
         
Where L is the Lagrangian, Ii the moment of inertia around the center of mass and V0 is a constant potential adjustment. V0 disappears in the equations of motions and only has an impact in the energy calculation.  We find the momenta p1 and p2 by differentiating the Lagrangian with respect to the generalized velocities, from which we find dtheta1/dt and dtheta2/dt. 
By using the Hamiltonian one can then find p1/dt and p2/dt. The exact equations are in the program, in the system function. 

The following Poincare sections visualizations were done with for m1=m2=1 and l1=l2=1
