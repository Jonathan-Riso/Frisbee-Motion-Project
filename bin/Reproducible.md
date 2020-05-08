# Reproducing Most Optimal Throwing Angle

#### Without Windspeed

No adjustments need to be made to the code if you desire you get the most optimal throwing angle without wind speed.

#### With Windspeed

```Java
   public static final double speedOfWind = 32;     // the speed of the wind is in miles per hour (mph)
   public static final double epsilon = .79;  // the direction of the wind vector in rads

```

The following parameters are found on lines **56** and **57** of the code. Adjust them to how ever you desire.
	    
```Java	    
   If there's no wind let speedOfWind and epsilon equal to 0
	
   Input the speed of the Wind and its direction regarding the frisbee (shortest angle)
   
   If the wind vector is directed in the third octant (between the south and the east) 
     - the user would have to replace the equation for windx by: windx = Wind(angleOfAttack)*Math.cos(Math.PI + epsilon);
       and windy by: windy = Wind(angleOfAttack)*Math.sin(Math.PI + epsilon);
       
   If the wind vector is directed in the fourth octant (bwtween the west and the south) 
     - the user would have to replace the equation for windx by: windx = Wind(angleOfAttack)*Math.cos(-epsilon);
       and windy by: windy = Wind(angleOfAttack)*Math.sin(-epsilon);
``` 

The direction of the wind can as well be changed by following the instructions from lines **46** to **54**.


#### Bisection Method

```java
while (Math.abs(c-a) > 0.0001)
		{
			if ((c-b) > (b-a))
			{
                x = (b + 0.38197*(c-b));
				
				if (flight(x) > flight(b))
				{
					c = x;
				}
				else
				{ 
                    a = b;
				    b = x;
				}
            }
            else
            {
                x = (b - 0.38197*(b-a));
				 
 				if (flight(x) < flight(b))
				{
				    c = b; 
					b = x; 
				}
				else 
				{
					a = x; 
				}
			}		
		}
```

The following is the Golden Section Search coded in java.

The *Flight(x)* method is calculating the final range of the frisbee thrown at an angle x. 

```java
public static double flight (double throwingAngle)
	{
		
		// Calculate length of arrays
		int imax = (int)(Tmax/DT);		    // Maximal index
		
		// Allocate the arrays	
		double[] alpha = new double[imax];   // Array for the angle of attack which change with time in rad
		double[] beta = new double[imax];    // Array for the direction of the velocity vector in rad
		
	        double[] x = new double[imax];       // x-position in m
		double[] y = new double[imax];	     // y-position in m
		
		double[] t = new double[imax];		 // time in sec 
		
		double[] vx = new double[imax];      // x-velocity in m/s
	        double[] vy = new double[imax];      // y-velocity in m/s
		
		double[] accx = new double[imax];    // x-acceleration in m/s^2
		double[] accy = new double[imax];    // y-acceleration in m/s^2
		
	
	
		//Calculate the initial angle of attack
		beta[0] = throwingAngle;
		alpha[0] = throwingAngle-beta[0]; 

        /**
			throwingAngle corresponds to beta, which is the angle at which the velocity vector points
			when the frisbee directly leaves the hand of the thrower 
			For this simulation we assume that at t = 0 beta and throwingAngle are equal.
        */		
		
		// Initialize the first step in time
		
		t[0] = 0;

		x[0] = xInitial;
		y[0] = yInitial;

		vx[0] = vInitial*Math.cos(beta[0]);
		vy[0] = vInitial*Math.sin(beta[0]);

		accx[0] = accelerationx (throwingAngle, beta[0], alpha[0], vx[0], vy[0]);
		accy[0] = accelerationy (throwingAngle, beta[0], alpha[0], vx[0], vy[0]);
		
		double range=0;
		
		for(int i=1; y[i-1]>=0; i++)   // Need to change size of arrays before continuing - the loop should stop when y=0
		{
			
			// Update time
			t[i] = t[i-1]+DT;

			// Update position
			x[i] = x[i-1] + vx[i-1]*DT;
			y[i] = y[i-1] + vy[i-1]*DT;
		   
			// Update velocity
			vx[i] = vx[i-1]+accx[i-1]*DT;
			vy[i] = vy[i-1]+accy[i-1]*DT;
			
			// Update the angle of attack 
			beta[i] = angleOfBeta(vx[i],vy[i]);
			alpha[i] = throwingAngle - beta[i]; 

			// Update acceleration
			accx[i] = accelerationx(throwingAngle, beta[i], alpha[i], vx[i], vy[i]);
		    accy[i] = accelerationy(throwingAngle, beta[i], alpha[i], vx[i], vy[i]);
			
		    range = x[i]; 
			
		}
		
		return -range; 
		
	}
```



# Reproducing Graphs

All graphs can be reproduced by running the FrisbeeMotionGraphing.java. Information on how to run this and which parameters to change are in the code and README.md file

#### Which Parameters and Equations to Change

```java
public static final double speedOfWind = 32;     // the speed of the wind is in miles per hour (mph)
public static final double epsilon = 1.56;  // the direction of the wind vector in rads
```

The above parameters are for adjusting the windspeed and direction. They are located on lines **57** and **58** of the code

```Java
public static final double throwingAngle = 0.16903; //intial throwing angle in rads
```

This parameter is for adjusting the initial angle the frisbee is thrown at. It is located on line **68** of the code.

```java
   public static double WindX (double angleOfAttack)
   {
        double windx = Wind(angleOfAttack)*Math.cos(Math.PI + epsilon);   // The magnitude of the Wind force in the x direction
        return windx;
   }
   public static double WindY (double angleOfAttack)
   {
	double windy = Wind(angleOfAttack)*Math.sin(Math.PI + epsilon);   // The magnitude of the Wind force in the y direction
	return windy;
   }
```
These two methods found at the lines **267** and **272** should also be modified by following the instructions
#### Working with The Graphs

```Java
//create PlotPanel
Plot2DPanel plot0 = new Plot2DPanel();
// define the legend position
plot0.addLegend("SOUTH");
// add a line plot to the PlotPanel
plot0.addLinePlot("FrisbeeMotion", xP, yP); //plotting the current position on x axis vs current position on y axis

//Creating Axis Titles and Graph Title
plot0.setAxisLabel(0,"Distance(m)");
plot0.getAxis(0).setLabelPosition(0.5,-0.1);
plot0.setAxisLabel(1,"Height(m)");
BaseLabel title0 = new BaseLabel("Frisbee Motion", Color.BLACK, 0.5, 1.1);
title0.setFont(new Font("Courier", Font.BOLD, 14));
plot0.addPlotable(title0);
```

The following is a snippet of the code used to create one the graphs, it is located on lines **199-212** of the code.

```java
//create PlotPanel
Plot2DPanel plot0 = new Plot2DPanel();
// define the legend position
plot0.addLegend("SOUTH");
// add a line plot to the PlotPanel
plot0.addLinePlot("FrisbeeMotion", xP, yP); //plotting the current position on x axis vs current position on y axis
```

The first line as stated creates the Plot Panel Object, specifically a 2d graph.

The second line sets the graphs position

The 3 line creates the function line based on two arrays inputted and a line name.

**Ex:**

```Java
plot.addLinePlot("LINE NAME", Array1,Array2);
```

It is important to note that Array1 will be the x value and Array2 is the y value.

The next section of the code snippet is concerned with labeling the graph. 

```java
//Creating Axis Titles and Graph Title
plot0.setAxisLabel(0,"Distance(m)");
plot0.getAxis(0).setLabelPosition(0.5,-0.1);
plot0.setAxisLabel(1,"Height(m)");
BaseLabel title0 = new BaseLabel("Frisbee Motion", Color.BLACK, 0.5, 1.1);
title0.setFont(new Font("Courier", Font.BOLD, 14));
plot0.addPlotable(title0);
```

The first and second lines sets the label of the X  and Y axis.

The final 3 lines create the title and place it on the graph.

#### Recreating Excel Graph(Figure 3)

To create the graph for Figure 3, we used excel and the horizontal ranges of the different angles. The results are located in date\ThrowingAngleVsDistanceData.txt
