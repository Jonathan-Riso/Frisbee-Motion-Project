### FrisbeeMotionOptimization

FrisbeeMotionOptimization when compiled and then ran will look for the best throwing angle in order to get maximum horizontal distance.

#### Wind Speed and Direction
If you desire to alter wind speed adjust the following constants:
speedOfWind:Sets the speed of the wind(mph)
epsilon:    Sets the direction of the wind regarding the frisbee (shortest angle) 



### FrisbeeMotionGraphing

When compiled and then run, it will output a 5 graphs based on the throwing Angle you decide (Angle must be radians and by default it is set to 10&deg;). The Graphs outputted will be in 5 separate windows and are; Horizontal Distance vs Vertical Distance, Time vs Velocity in x and y directions, and Time vs Acceleration in the x and y directions.

##### FrisbeeMotionData.txt

FrisbeeMotionGraphing when compiled will output the text file *FrisbeeMotionData.txt*. This text file will hold all the results in the arrays used to plot the graphs, which are; Time, Horizontal Distance, Vertical Distance, Velocity in X direction, Velocity in Y direction, Acceleration in X direction and Acceleration in Y direction. When the program is run the file is wiped and the new results are added, should you desire to save these results move them into the data folder and save them.
