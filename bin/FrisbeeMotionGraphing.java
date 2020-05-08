import java.io.*;
import java.awt.*;
import javax.swing.*;
import java.util.Scanner;
import org.math.plot.*;
import org.math.plot.plotObjects.*;
import java.util.Arrays;


/*
   Motion of a Frisbee in flight

   Author: Charles Frederic Michel and Jonathan Risco
   Current version written: May 2019
   Description: Calculates and plot the motion of Frisbee taking into account
                drag, spin rate and lift for a specific angle of attack.
*/

public class FrisbeeMotionGraphing {
	
    //////////////////////
    // Declare constants//
    //////////////////////
	public static final double DT = 1.e-3;    //time interval in s
    public static final double Tmax = 50.;    //max time in s
	
	// Parameters of the equations 
    public static final double mass = 0.175;   // in (kg)
	public static final double gravity = 9.81; // in (m/s^2) 
	public static final double p = 1.23;	   // density of air  in (kg/m^3)
	public static final double d = 0.26;       // diameter of the frisbee in (m)
    public static final double Area = Math.PI*Math.pow(d/2.0,2); // Area of the frisbee in (m^2)     
	public static final double Cd0 = 0.08;     // the form drag 
    public static final double Cda = 2.72;     // the induced drag 
    public static final double alpha0 = 0.07;  // the angle of attack that produce the least lift (in rad) (-4degree is 0.0.698 rad)
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

  
  
  
   ///////////////////////////////
   //EDITING THE INITAL THROWNING ANGLE//
   ///////////////////////////////
   //Change the angle to what ever you desire. It must be in rads
   public static final double throwingAngle = 0.16903; //intial throwing angle in rads
   
	
    // Start main method
    public static void main(String[] args){
		
				
		
		// Calculate length of arrays
		int N = (int)(Tmax/DT);		    // Maximal index
		
		// Allocate the arrays	
		double[] alpha = new double[N];   // Array for the angle of attack which change with time in rad
		double[] beta = new double[N];    // Array for the direction of the velocity vector in rad
		
	    double[] x = new double[N];       // x-position in m
		double[] y = new double[N];	     // y-position in m
		
		double[] t = new double[N];		 // time in sec 
		
		double[] vx = new double[N];      // x-velocity in m/s
	    double[] vy = new double[N];      // y-velocity in m/s
		
		double[] accx = new double[N];    // x-acceleration in m/s^2
		double[] accy = new double[N];    // y-acceleration in m/s^2
	
		
		//Calculate the initial angle of attack
		beta[0] = throwingAngle;
		alpha[0] = throwingAngle - beta[0]; 

        /**
			throwingAngle corresponds to beta, which is the angle at which the velocity vector points
			when the frisbee directly leaves the hand of the thrower 
			For this simulation we assume that at t = 0 beta and throwingAngle are equal.
        */			
		
		
		////////////////
		//Opening File//
		///////////////
		String filename = "FrisbeeFlightData.txt";
        PrintWriter outputFile = null;
        try
        {
            outputFile = new PrintWriter(new FileOutputStream(filename,false));
        }
        catch(FileNotFoundException e)
        {
            System.out.println("File error.  Program aborted.");
            System.exit(0);
        }

		
		
		
		
		
		// Initialize the first step in time
		
		t[0] = 0;

		x[0] = xInitial;
		y[0] = yInitial;

		vx[0] = vInitial*Math.cos(beta[0]);
		vy[0] = vInitial*Math.sin(beta[0]);

		accx[0] = accelerationx (throwingAngle, beta[0], alpha[0], vx[0], vy[0]);
		accy[0] = accelerationy (throwingAngle, beta[0], alpha[0], vx[0], vy[0]);
		
		outputFile.printf("%2.4f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\n" ,t[0],x[0],y[0],vx[0],vy[0],accx[0],accy[0]);
	
		
		for(int i=1; y[i-1]>0 ;i++)   // Need to change size of arrays before continuing - the loop should stop when y=0
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
			outputFile.printf("%2.4f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\n" ,t[i],x[i],y[i],vx[i],vy[i],accx[i],accy[i]);
		}
		
		 ////////////////////
		 //Closing the File//
		 ////////////////////
		 outputFile.close();
		
		int cnt = 0;
		for(int i = 0; y[i]>0; i++)//counting number of value in arrays
			cnt++;

		//Creating Plotting Arrays
		double[] tP = new double [cnt];
		double[] xP = new double[cnt];
		double[] yP = new double[cnt];
		double[] vxP = new double[cnt];
		double[] vyP = new double[cnt];
		double[] accxP = new double[cnt];
		double[] accyP = new double[cnt];
		
        //Filling Arrays
		for(int i = 1; y[i-1]>0; i++){
			tP[i-1] = t[i-1];
			xP[i-1] = x[i-1];
			yP[i-1] = y[i-1];
			vxP[i-1] = vx[i-1];
			vyP[i-1] = vy[i-1];
			accxP[i-1] = accx[i-1];
			accyP[i-1] = accy[i-1];
		}
      /*
		Graphing
		Following section creates 5 plots with data from the plotting arrays above
	  */
	  
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
		

		//create PlotPanel
		Plot2DPanel plot1 = new Plot2DPanel();
		// define the legend position
		plot1.addLegend("SOUTH");
		// add a line plot to the PlotPanel
		plot1.addLinePlot("FrisbeeSpeedX", tP, vxP); // plotting time on x axis vs velocity of x
		
		//Creating Axis Titles and Graph Title
		plot1.setAxisLabel(0,"Time(s)");
        plot1.getAxis(0).setLabelPosition(0.5,-0.1);
        plot1.setAxisLabel(1,"Velocity(m/s)");
        BaseLabel title1 = new BaseLabel("Velocity in x direction vs time", Color.BLACK, 0.5, 1.1);
        title1.setFont(new Font("Courier", Font.BOLD, 14));
        plot1.addPlotable(title1);
		
		
		
		
		
		//create PlotPanel
		Plot2DPanel plot2 = new Plot2DPanel();
		// define the legend position
		plot2.addLegend("SOUTH");
		// add a line plot to the PlotPanel
		plot2.addLinePlot("FrisbeeSpeedY", tP, vyP);// plotting time on x axis vs velocity of y
		
		//Creating Axis Titles and Graph Title
		plot2.setAxisLabel(0,"Time(s)");
        plot2.getAxis(0).setLabelPosition(0.5,-0.1);
        plot2.setAxisLabel(1,"Velocity(m/s)");
        BaseLabel title2 = new BaseLabel("Velocity in y direction vs time", Color.BLACK, 0.5, 1.1);
        title2.setFont(new Font("Courier", Font.BOLD, 14));
        plot2.addPlotable(title2);
		
		
		//create PlotPanel
		Plot2DPanel plot3 = new Plot2DPanel();
		// define the legend position
		plot3.addLegend("SOUTH");
		// add a line plot to the PlotPanel
		plot3.addLinePlot("FrisbeeAccX", tP, accxP);// plotting time on x axis vs acceleration of x

		//Creating Axis Titles and Graph Title
		plot3.setAxisLabel(0,"Time(s)");
        plot3.getAxis(0).setLabelPosition(0.5,-0.1);
        plot3.setAxisLabel(1,"Acceleration(m/s^2)");
        BaseLabel title3 = new BaseLabel("Acceleration in x direction vs time", Color.BLACK, 0.5, 1.1);
        title3.setFont(new Font("Courier", Font.BOLD, 14));
        plot3.addPlotable(title3);
		
		//create PlotPanel
		Plot2DPanel plot4 = new Plot2DPanel();
		// define the legend position
		plot4.addLegend("SOUTH");
		// add a line plot to the PlotPanel
		plot4.addLinePlot("FrisbeeMotionAccY", tP, accyP);// plotting time on x axis vs acceleration of y
		
		//Creating Axis Titles and Graph Title
		plot4.setAxisLabel(0,"Time(s)");
        plot4.getAxis(0).setLabelPosition(0.5,-0.1);
        plot4.setAxisLabel(1,"Acceleration(m/s^2)");
        BaseLabel title4 = new BaseLabel("Acceleration in y direction vs time", Color.BLACK, 0.5, 1.1);
        title4.setFont(new Font("Courier", Font.BOLD, 14));
        plot4.addPlotable(title4);
		
		
		//Creates 5 panels where 5 graphs are displayed
		
		// put the PlotPanel in a JFrame like a JPanel
		JFrame frame1 = new JFrame("Velocity in x direction vs time");
		frame1.setSize(600, 600);
		frame1.setContentPane(plot1);
		frame1.setVisible(true);
		
		JFrame frame2 = new JFrame("Velocity in y direction vs time");
		frame2.setSize(600, 600);
		frame2.setContentPane(plot2);
		frame2.setVisible(true);
		
		JFrame frame3 = new JFrame("Acceleration in x direction vs time");
		frame3.setSize(600, 600);
		frame3.setContentPane(plot3);
		frame3.setVisible(true);
		
		JFrame frame4 = new JFrame("Acceleration in y direction vs time");
		frame4.setSize(600, 600);
		frame4.setContentPane(plot4);
		frame4.setVisible(true);
		
		JFrame frame0 = new JFrame("X Distance vs Y Height");
		frame0.setSize(800, 600);
		frame0.setContentPane(plot0);
		frame0.setVisible(true);
	}//end of main method
	/**
	@param liftf it's the value of the lift's force exerted on the COP
	@return accy it's the value of the acceleration in the y-direction 
	        of the frisbee when it has a speed of vy 
    */	
    public static double accelerationy (double throwingAngle, double beta, double angleOfAttack, double vx, double vy)
    {
		double accy = ((dragforceY(angleOfAttack, beta, vx, vy) + liftforceY(throwingAngle, angleOfAttack, vx, vy) + WindY(angleOfAttack))/mass)-gravity; 
        return accy;
    }	

	/**
	@param dragf it's the value of the drag's force opposing the direction of flight
	@return accx it's the value of the acceleration in the x-direction 
	        of the frisbee when it has a speed of vx 
    */	
	public static double accelerationx (double throwingAngle, double beta, double angleOfAttack, double vx, double vy)
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
		double Cl = Cl0 + (Cla*angleOfAttack);                 // Equation for the lift coefficient with change time
		double v = Math.sqrt(Math.pow(vx,2)+ Math.pow(vy,2));  // The norm of the velocity 
		double liftf = 0.5*p*Math.pow(v,2)*Area*Cl;            // Equation to find the lift force 
		
	    double liftfy = liftf*Math.sin(beta + (Math.PI/2));       // The force of the lift in x-direction 
		
		return liftfy;
	}
    public static double angleOfBeta (double vx, double vy)
	{
		double beta = Math.atan(vy/vx);                                      //Equation to find the angle of the velocity vector in respect to the x-axis 
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