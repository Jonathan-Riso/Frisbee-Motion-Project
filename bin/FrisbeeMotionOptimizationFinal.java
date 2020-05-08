/**
   Motion of a Frisbee in flight
   Author: Charles Frederic Michel and Jonathan Risco
   Current version written: May 2019
   Description: By using the bisection method, Optimizes for the best
                throwing angle for a Frisbee to go the furthest horizontal distance. 
                It takes into account drag, spin rate and lift for a specific
                angle of attack.
*/

public class FrisbeeMotionOptimizationFinal
{
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Declare constants
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
	public static final double DT = 1.e-3;    //time interval in s
        public static final double Tmax = 50.;    //max time in s
	
	/**
		All of the public static final parameters were taken from "The Aerodynamics of Frisbee Flight" by Baumback
		since the world record for the longest frisbee flight is  20.52 s made by Christer Fuglesang (ESA Astronaut) (SWE)
		the maximum time that was set is largely enough
	*/
	
	// Parameters of the equations 
        public static final double mass = 0.175;   // in (kg)
	public static final double gravity = 9.81; // in (m/s^2) 
	public static final double p = 1.23;	   // density of air  in (kg/m^3)
	public static final double d = 0.26;       // diameter of the frisbee in (m)
        public static final double Area = Math.PI*Math.pow(d/2.0,2); // Area of the frisbee in (m^2)     
	public static final double Cd0 = 0.08;     // the form drag 
        public static final double Cda = 2.72;     // the induced drag 
        public static final double alpha0 = 0.0698;  // the angle of attack that produce the least lift (in rad) (-4 degree is 0.0.698 rad)
	public static final double Cl0 = 0.15; 	   // lift  y-intercept when alpha = 0
	public static final double Cla = 1.4;      // lift  the slope of the graph
	
        // Initial position of the frisbee in x and y
        public static final double xInitial = 0;	      // at 0
	public static final double yInitial = 1.5;	      // at which height the frisbee is thrown (1.5 m from the ground)

        // Initial speed of the frisbee 
        public static final double vInitial = 14.0;        // in (m/s^2)   
   
        /**
	    If there's no wind let speedOfWind and epsilon equal to 0
		
            Input the speed of the Wind and it's direction regarding the frisbee (shortest angle)
		If the wind vector is directed in the third octant (between the south and the east) 
		 - the user would have to replace the eqaution for windx by: windx = Wind(angleOfAttack)*Math.cos(Math.PI + epsilon);
		                                               and windy by: windy = Wind(angleOfAttack)*Math.sin(Math.PI + epsilon);
		If the wind vector is directed in the fourth octant (bewteen the west and the south) 
		 - the user would have to replace the eqaution for windx by: windx = Wind(angleOfAttack)*Math.cos(-epsilon);
		                                               and windy by: windy = Wind(angleOfAttack)*Math.sin(-epsilon);
        */
        public static final double speedOfWind = 0;     // the speed of the wind is in miles per hour (mph)
        public static final double epsilon = 0;  // the direction of the wind vector in rads

	
        // Start main method
        public static void main(String[] args)
        {
		
		/**
		  Golden Section Search method from Baywatch
		  the values from a to c correspond to the angle between [0, Pi/2[
		*/
		double a = 0.000; 
		double c = 1.56;
		double b = (a + 0.38197*(c-a));  
		double x; 
		
		
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
	System.out.printf("The angle at which the frisbee should be thrown to have maximum range is %3.5f\n",(a+c)/2);
        System.out.printf("The maximum range is %3.4f\n",-flight((a+c)/2));
		
	} // end of main method 
	
	/**
	@param pitchAngle it's the angle at which the person is throwing the angle calculated from the x-axis
	@return -R which is the value of the range of the frisbee when thrown from thnrowingAngle
	
	We need to return a negative R value in order to find the Golden Section Search which finds the absolute minimum
	It's same principle as in the Baywatch code but the difference is that instead of having an equation we need to 
        numerically find the values of f(x) => Range(angle)
	*/
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
		
    /**
	@param liftf it's the value of the lift's force exerted on the COP
	@return accy it's the value of the acceleration in the y-direction 
	        of the frisbee when it has a speed of vy 
    */	
       public static double accelerationy (double initialAngle, double beta, double angleOfAttack, double vx, double vy)
       {
		double accy = ((dragforceY(angleOfAttack, beta, vx, vy) + liftforceY(angleOfAttack, beta, vx, vy) + WindY(angleOfAttack))/mass) - gravity; 
                return accy;
       }	
	
    /**
	@param dragf it's the value of the drag's force opposing the direction of flight
	@return accx it's the value of the acceleration in the x-direction 
	        of the frisbee when it has a speed of vx 
    */	
	public static double accelerationx (double initialAngle, double beta, double angleOfAttack, double vx, double vy)
	{
		double accx = (dragforceX(angleOfAttack, beta, vx, vy) + liftforceX(angleOfAttack, beta, vx, vy) + WindX(angleOfAttack))/mass;   //Equation to find the acceleration in x direction
                return accx;
	}
	
	public static double dragforceX (double angleOfAttack, double beta, double vx, double vy)
	{
	    double Cd = Cd0 + (Cda*Math.pow((angleOfAttack - alpha0),2));  // Equation for the drag coefficient with change time
		double v = Math.sqrt(Math.pow(vx,2)+ Math.pow(vy,2));          // The norm of the velocity 
		double dragf = 0.5*(Cd*p*Area*(Math.pow(v,2)));               // Equation to find the drag force 
		
		double dragfx = dragf*Math.cos(beta + Math.PI);               // The force of the drag in x-direction 
		return dragfx;
	}
	public static double dragforceY (double angleOfAttack, double beta, double vx, double vy)
	{
	        double Cd = Cd0 + (Cda*Math.pow((angleOfAttack - alpha0),2));  // Equation for the drag coefficient with change time
		double v = Math.sqrt(Math.pow(vx,2)+ Math.pow(vy,2));          // The norm of the velocity 
		double dragf = 0.5*(Cd*p*Area*(Math.pow(v,2)));               // Equation to find the drag force 
		
		double dragfy = dragf*Math.sin(beta + Math.PI);               // The force of the drag in y-direction 
		return dragfy;
	}
	public static double liftforceX (double angleOfAttack, double beta, double vx, double vy)
	{
		double Cl = Cl0 + (Cla*angleOfAttack);                 // Equation for the lift coefficient with change time
		double v = Math.sqrt(Math.pow(vx,2)+ Math.pow(vy,2));  // The norm of the velocity 
		double liftf = 0.5*p*Math.pow(v,2)*Area*Cl;            // Equation to find the lift force 
		
		double liftfx = liftf*Math.cos(beta + (Math.PI/2));         // The force of the lift in x-direction 
		return liftfx;
	}
	public static double liftforceY (double angleOfAttack, double beta, double vx, double vy)
	{
		double Cl = Cl0 + (Cla*angleOfAttack);                             // Equation for the lift coefficient with change time
		double v = Math.sqrt(Math.pow(vx,2)+ Math.pow(vy,2));              // The norm of the velocity 
		double liftf = 0.5*p*Math.pow(v,2)*Area*Cl;                        // Equation to find the lift force 
		
	        double liftfy = liftf*Math.sin(beta + (Math.PI/2));                // The force of the lift in x-direction 
		
		return liftfy;
	}
        public static double angleOfBeta (double vx, double vy)
  	{
		double beta = Math.atan(vy/vx);                                    //Equation to find the angle of the velocity vector in respect to the x-axis 
		return beta;
	} 
	public static double Wind (double angleOfAttack)
	{
		double Cd = Cd0 + (Cda*Math.pow((angleOfAttack - alpha0),2));      // Equation for the drag coefficient with change of angle of attack 
		double pressure = 0.00256*Math.pow(speedOfWind,2);                 // The pressure exerted by the wind on the frisbee
		double force = Area*pressure*Cd;                                   // The force of the wind 
		return force; 
	}
	public static double WindX (double angleOfAttack)
	{
		double windx = Wind(angleOfAttack)*Math.cos(Math.PI + epsilon);     // The magnitude of the Wind force in the x direction
		return windx;
	}
	public static double WindY (double angleOfAttack)
	{
		double windy = Wind(angleOfAttack)*Math.sin(Math.PI + epsilon);              // The magnitude of the Wind force in the y direction
		return windy;
	}
} //end of the program
