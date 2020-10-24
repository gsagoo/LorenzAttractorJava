# LorenzAttractorJava
Java code that produces a GUI window to illustrate the Lorenz attractor and vary its parameters.

# Compiling source code and running the class.


1. Compile the Lorenz.java source code by typing the following at the command prompt 

  ~~~ 
  javac Lorenz.java 
  ~~~

(please take note that there is a "c" on the end of javac)

2. The Lorenz.class executable file will then be compiled which is then run by typing 

  ~~~
java Lorenz
  ~~~

this will produce a window with the Lorenz button and buttons to control the view.


# Lorenz Java Class

Parameters
^^^^^^^^^^
(a, b, c) a parameters in the coupled Lorenz equations.
Changing them will result in different attractors/
fixed points/unstable spirals.


n-plane
^^^^^^^
The Lorenz attractor is projected onto the n-plane.
The n-plane is the unique 2D plane that contains the point
(n0. n1, n2) and is perpendicular to this vector too.

Iterations
^^^^^^^^^^
Number of points plotted.

Rotation
^^^^^^^^
Anticlockwise rotation about at axis that is perpendicular
to the screen. Measured in degrees.

Redraw
^^^^^^
Draws the Lorenz trajectories according to the above values.

Euler/R-K
^^^^^^^^^
Changes the scheme used in the finite differencing. The choices are:
1. Runge-Kutta (4th order accurate) and 
2. Euler's method ( 1st order accurate). 

Orbit 90
^^^^^^^^
Generates a small animation (10-20 frames) which shows the
attractor as observed from the postion (n0, n1, n2) when the 
orbit axis throught the orbit point P is being circled . 
In this case the attractor is projected onto a 2D plane that:
1. Contains the point (n0. n1, n2) and
2. is parrallel to the orbit axis passing the through the P and
3. is perpendicular to the foot of the perpendicular from (n0. n1, n2) to 
   the orbit axis which passes through the P.

The orbit axis is automatically
normalised to a unit vector by the program. 
The Orbit point is point about which the orbit takes place, the
default value is set to the middle of the bounding box that 
encapsulates the attractor.

NOTE: The values of orbit-axis and orbit point have no affect on 'redraw'
