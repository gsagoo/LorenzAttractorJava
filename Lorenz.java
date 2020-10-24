///////////////////////////////////////////////////
///// 
///// Author: G Sagoo
///// Shows the Lorenz Attractor in a window and allows simple explore
///// Compile with: java Lorenz
/////  
///////////////////////////////////////////////////

import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.lang.*;
import java.util.*;
import java.util.Arrays;
import java.io.*;
import java.util.regex.Pattern; 
import java.util.regex.Matcher; 
import java.awt.event.*;
import javax.swing.border.LineBorder;
import java.awt.event.ActionListener;


public class Lorenz extends JPanel implements ActionListener
{
    private long iterations = 26000;
    private long iterationLapse = Math.min(500, iterations/2);  //Must be less than variable iterations
    private double  h = 0.01, a = 10.0, b = 28.0, c = 8.0 / 3.0;
    private static final double TWOPI = 2.0*Math.PI;
    private static final double TINY = 0.000001;
    private double xmax = 0, ymax = 0, xmin = 0, ymin = 0;
    private double [] u = new double[3]; 
    private double [] v = new double[3];
    private double [] n = new double[3];
    private double [] p = new double[3]; 
    private double [] l_max = new double[3];
    private double [] l_min = new double[3];  
    private int hShift = 0,  vShift = 0, orbit = 0;
    private double dZoom = 0, rotationAngle = 0.0;
    private boolean reDraw = false, method_rk = false;
    public int getOrbit(){return(orbit);}
   
    //Frames and Panels
    JFrame frame = new JFrame(); 
    JPanel myContainerPanel = new JPanel();
    JPanel myControlsContainer = new JPanel();
    JPanel myControlPanel = new JPanel();
    JPanel myZoomPanel = new JPanel();
    
    //Make some Buttons
    JButton bZoomIn = new JButton("Zoom+");
    JButton bZoomOut = new JButton("Zoom-");
    JButton bUp = new JButton("Up");
    JButton bDown = new JButton("Down");
    JButton bLeft = new JButton("Left");
    JButton bRight = new JButton("Right");
    JButton bDraw = new JButton("Redraw");
    JButton bOrbit = new JButton("<html>Orbit 90<font size=2><sup>o</sup></font></html>");
    JButton bReset = new JButton("Reset");
    JButton bMethod = new JButton("Euler");
   
    JTextField n0Value = new JTextField(""+n[0], 3); 
    JTextField n1Value = new JTextField(""+n[1], 3);
    JTextField n2Value = new JTextField(""+n[2], 3);
    JTextField aValue = new JTextField(""+a, 3); 
    JTextField bValue = new JTextField(""+b, 3);
    JTextField cValue = new JTextField(""+c, 3);
    JTextField rValue = new JTextField(""+rotationAngle, 3); 
    JTextField iValue = new JTextField(""+iterations, 4);
    JTextField oA0Value = new JTextField(""+v[0], 3); 
    JTextField oA1Value = new JTextField(""+v[1], 3);
    JTextField oA2Value = new JTextField(""+v[2], 3);
    JTextField oP0Value = new JTextField(""+p[0], 3); 
    JTextField oP1Value = new JTextField(""+p[1], 3);
    JTextField oP2Value = new JTextField(""+p[2], 3);
    JTextField orbitAngleValue = new JTextField("0.0", 3);
    
    JLabel nLabel = new JLabel("n-plane vector: ");
    JLabel pLabel = new JLabel("Parameters (a,b,c):");
    JLabel rLabel = new JLabel("<html>Rotate<font size=2><sup>o</sup></font> :</html>"); //Embedded HTML!
    JLabel iLabel = new JLabel("#Iterations: ");
    JLabel oAxisLabel = new JLabel("Orbit-axis: ", SwingConstants.RIGHT); //Right Justify
    JLabel oPointLabel = new JLabel("Orbit-point: ", SwingConstants.RIGHT);
    
    // Constructor
    Lorenz()
    {     
         //Create a frame to hold the Jpanel
        final int FRAME_WIDTH = 500;
        final int FRAME_HEIGHT = 400;
        n[0] = 8; n[1] = -4; n[2] = 3;
        n0Value.setText(""+n[0]); n1Value.setText(""+n[1]);n2Value.setText(""+n[2]);
        orbitAngleValue.setEditable(false);
        
        frame.setSize( FRAME_WIDTH, FRAME_HEIGHT );
        frame.setTitle("Lorenz Attractor - Projected in n-plane:  n.(x, y, z) = |n|");      
        frame.setDefaultCloseOperation( JFrame.EXIT_ON_CLOSE );
        
        //Frame colours are also reset inside paintComponent (to make repainting clean)
        frame.setBackground( Color.white );
        frame.setForeground( Color.white );
        
        //Adds buttons and labels to top control panel
        GridBagConstraints constraintsGridBag = new GridBagConstraints();
        GridBagLayout gridbag = new GridBagLayout();
        myControlPanel.setLayout(gridbag);
        constraintsGridBag.fill = GridBagConstraints.NONE;
        constraintsGridBag.weightx=0.5;
        
        //Row 1
        constraintsGridBag.gridx = 0;
        constraintsGridBag.gridy = 0;
        constraintsGridBag.gridwidth = 3;
        constraintsGridBag.anchor = GridBagConstraints.EAST; //Right Justify
        myControlPanel.add(pLabel, constraintsGridBag);
        constraintsGridBag.gridwidth = 1;
        constraintsGridBag.gridx = 3;
        myControlPanel.add(aValue, constraintsGridBag);
        constraintsGridBag.gridx = 4;
        myControlPanel.add(bValue, constraintsGridBag);
        constraintsGridBag.gridx = 5;
        myControlPanel.add(cValue, constraintsGridBag);
        constraintsGridBag.gridx = 6;
        constraintsGridBag.gridwidth = 2;
        constraintsGridBag.ipady = -9;
        constraintsGridBag.anchor = GridBagConstraints.CENTER;
        myControlPanel.add(bReset, constraintsGridBag);
        constraintsGridBag.anchor = GridBagConstraints.EAST;
        constraintsGridBag.ipady = 0;
        constraintsGridBag.gridwidth = 1;
        constraintsGridBag.gridx = 8;
        myControlPanel.add(oAxisLabel, constraintsGridBag);
        constraintsGridBag.gridx = 9;
        myControlPanel.add(oA0Value, constraintsGridBag);
        constraintsGridBag.gridx = 10;
        myControlPanel.add(oA1Value, constraintsGridBag);
        constraintsGridBag.gridx = 11;
        myControlPanel.add(oA2Value, constraintsGridBag);
        
        //Row2
        constraintsGridBag.gridy = 1;
        constraintsGridBag.gridx = 0;
        constraintsGridBag.gridwidth = 3;
        myControlPanel.add(nLabel, constraintsGridBag);
        constraintsGridBag.gridwidth = 1;
        constraintsGridBag.gridx = 3;
        myControlPanel.add(n0Value, constraintsGridBag);
        constraintsGridBag.gridx = 4;
        myControlPanel.add(n1Value, constraintsGridBag);
        constraintsGridBag.gridx = 5;
        myControlPanel.add(n2Value, constraintsGridBag);
        
        constraintsGridBag.gridx = 6;
        constraintsGridBag.gridwidth = 2;
        constraintsGridBag.ipady = -9;
        constraintsGridBag.anchor = GridBagConstraints.CENTER;
        myControlPanel.add(bMethod, constraintsGridBag);
        constraintsGridBag.anchor = GridBagConstraints.EAST;
        constraintsGridBag.ipady = 0;
        constraintsGridBag.gridwidth = 1;
        
        constraintsGridBag.gridx = 8;
        myControlPanel.add(oPointLabel, constraintsGridBag);
        constraintsGridBag.gridx = 9;
        myControlPanel.add(oP0Value, constraintsGridBag);
        constraintsGridBag.gridx = 10;
        myControlPanel.add(oP1Value, constraintsGridBag);
        constraintsGridBag.gridx = 11;
        myControlPanel.add(oP2Value, constraintsGridBag);
        
        //Row 3
        constraintsGridBag.gridy = 2;
        constraintsGridBag.gridx = 0;
        constraintsGridBag.gridwidth = 2;
        constraintsGridBag.anchor = GridBagConstraints.EAST;
        myControlPanel.add(iLabel, constraintsGridBag);
        constraintsGridBag.gridwidth = 1;
        constraintsGridBag.gridx = 2;
        myControlPanel.add(iValue, constraintsGridBag);
        constraintsGridBag.gridx = 3;
        constraintsGridBag.gridwidth = 2;
        myControlPanel.add(rLabel, constraintsGridBag);
        constraintsGridBag.gridwidth = 1;
        constraintsGridBag.gridx = 5;
        myControlPanel.add(rValue, constraintsGridBag);
        constraintsGridBag.gridx = 6;
        constraintsGridBag.gridwidth = 2;
        constraintsGridBag.ipady = -9;
        constraintsGridBag.anchor = GridBagConstraints.CENTER;
        myControlPanel.add(bDraw, constraintsGridBag);
        constraintsGridBag.ipady = 0;
        constraintsGridBag.gridwidth = 1;   
        constraintsGridBag.gridx = 9;
        constraintsGridBag.gridwidth = 2;
        constraintsGridBag.ipady = -9;
        myControlPanel.add(bOrbit, constraintsGridBag);
        constraintsGridBag.ipady = 0;
        constraintsGridBag.gridwidth = 1;
        constraintsGridBag.gridx = 11;
        myControlPanel.add(orbitAngleValue, constraintsGridBag);
        
            
        //Adding components to zoompanel
        GridBagConstraints gbc = new GridBagConstraints();
        GridBagLayout gb = new GridBagLayout();
        myZoomPanel.setLayout(gb);
        gbc.fill = GridBagConstraints.NONE;
        gbc.weightx = 0;
        gbc.ipady = -10;
        myZoomPanel.add(bUp, gbc);
        myZoomPanel.add(bDown, gbc);
        myZoomPanel.add(bLeft, gbc);
        myZoomPanel.add(bRight, gbc);
        myZoomPanel.add(bZoomIn, gbc);
        myZoomPanel.add(bZoomOut, gbc);
        
        //Adding Listeners on Buttons
        bUp.addActionListener(this);
        bDown.addActionListener(this);
        bLeft.addActionListener(this);
        bRight.addActionListener(this);
        bZoomIn.addActionListener(this);
        bZoomOut.addActionListener(this);
        bDraw.addActionListener(this);
        bOrbit.addActionListener(this);
        bReset.addActionListener(this);
        bMethod.addActionListener(this);
        
        myControlPanel.setOpaque(true);
        myContainerPanel.setLayout(new BorderLayout());
        myContainerPanel.add(myControlPanel, BorderLayout.NORTH);
        myContainerPanel.add(myZoomPanel, BorderLayout.SOUTH);
        myContainerPanel.add(this, BorderLayout.CENTER);
        
        myContainerPanel.setDoubleBuffered(true);
        
        frame.add(myContainerPanel);
        frame.setVisible(true);
    }
   
   private void rotateVectorAboutAxis(double[] vec, double rA, double [] ax)
   {
    
    //vec = vector to be rotated
    //rA = rotation angle
    //ax = unit vector representing the axis to be rotated about
    double [][] rotMatrix = new double[3][3];
    double [] rotatedVector = new double[3];
    int i = 0, j = 0;
    
    rotMatrix[0][0] = ax[0]*ax[0] + (1 - ax[0]*ax[0] )*Math.cos(rA);
    rotMatrix[0][1] = ax[0]*ax[1]*(1 - Math.cos(rA)) + ax[2]*Math.sin(rA);
    rotMatrix[0][2] = ax[0]*ax[2]*(1 - Math.cos(rA)) - ax[1]*Math.sin(rA);
    
    rotMatrix[1][0] = ax[0]*ax[1]*(1 - Math.cos(rA)) - ax[2]*Math.sin(rA); 
    rotMatrix[1][1] = ax[1]*ax[1] + (1 - ax[1]*ax[1] )*Math.cos(rA);
    rotMatrix[1][2] = ax[1]*ax[2]*(1 - Math.cos(rA)) + ax[0]*Math.sin(rA); 
    
    rotMatrix[2][0] = ax[0]*ax[2]*(1 - Math.cos(rA)) + ax[1]*Math.sin(rA);
    rotMatrix[2][1] = ax[1]*ax[2]*(1 - Math.cos(rA)) - ax[0]*Math.sin(rA);
    rotMatrix[2][2] = ax[2]*ax[2] + (1 - ax[2]*ax[2] )*Math.cos(rA);
    
    //Check if pointers are used
    //Multiply the vector by matrix
    for( i = 0; i < 3; i++)
       for( j = 0, rotatedVector[i] = 0; j < 3; j++)
            rotatedVector[i] +=rotMatrix[i][j]*vec[j];
    for( i = 0; i < 3; i++)
        vec[i] = rotatedVector[i];
   }
   
    private void lorenzBoundingBox(double x, double y, double z)
    {
        double[] dtemp = new double[3];
        dtemp[0] = x; dtemp[1] = y; dtemp[2] = z;
    
        for( int i = 0; i < 3; i++) 
        if( l_max[i] < dtemp[i])
            l_max[i] = dtemp[i];
        else if (l_min[i] > dtemp[i] ) 
            l_min[i] = dtemp[i];
    }
   
   private double roundD( double a, int places)
   {
        return( ((double)((long)(Math.pow(10,places)*a)))/Math.pow(10,places) );
   }
   
   private double vectorCrossProduct(double [] a, double [] b, double [] c)
   {
      c[0] = a[1]*b[2]-a[2]*b[1];
      c[1] = a[2]*b[0]-a[0]*b[2];
      c[2] = a[0]*b[1]-a[1]*b[0];
      return(c[0]+c[1]+c[2]);
   } 
  
   private double vectorMag(double [] a)
   {
    double mag = 0;
    for (int i = 0; i<3; i++)
        mag += a[i]*a[i];
    mag = Math.abs( Math.sqrt(mag) );
    return(mag);    
   }
   
   private double vectorDot(double [] a, double [] b)
   {
    double dotProduct = 0;
    for (int i = 0; i<3; i++)
        dotProduct += a[i]*b[i];
    return(dotProduct);    
   }
   
   private void print_n_oA_oP()
   {
        n0Value.setText(""+roundD(n[0],3)); n1Value.setText(""+roundD(n[1],3));n2Value.setText(""+roundD(n[2],3));
        oP0Value.setText(""+roundD(p[0],3)); oP1Value.setText(""+roundD(p[1],3));oP2Value.setText(""+roundD(p[2],3));
        oA0Value.setText(""+roundD(v[0],3)); oA1Value.setText(""+roundD(v[1],3));oA2Value.setText(""+roundD(v[2],3));
   }
   
    //The driver
   public static void main(String[] args)
   {   	   	            
        Lorenz myLorenz = new Lorenz();
        //Need this infinite loop to check the orbit button status.
        //If the orbit button is used to trigger the orbit animation the listener on the 
        //orbit button will  causes the screen to freeze whilst the GUI repaints 
        //each animation frame; this is a quick work around.
        
        while(true) 
        {
            try {Thread.sleep(1200);} 
                    catch (InterruptedException e){}
            if(myLorenz.getOrbit() == 1)                      
                myLorenz.orbitButtonPress();
        }
    }
   
   private void h_f_RungeKutta(double x, double y, double z, double[] k_store)
   {
     k_store[0] =  h*a*(y -x);
     k_store[1] =  h*(x*(b -z) - y);
     k_store[2] =  h*(x*y-c*z);
   }
    
    private void calculatePoints( Graphics2D gContext )
    {
        int i = 0, j = 0, fixedPointsCount = 0;
        double [] fixedPoint = new double[3];
        double x0,y0,z0,x1,y1,z1;
        double nxmin = 0, nxmax = 0, nymin = 0, nymax = 0;
        double pxold, pyold, nx,ny;
        double px = 0, py = 0;
        Dimension dimJPanel = getSize();
        double pHeight = dimJPanel.height; 
        double pWidth = dimJPanel.width;
        double dtemp = 0;
        boolean rangesSet = false;
        double [] k1  = new double[3];
        double [] k2  = new double[3];
        double [] k3  = new double[3];
        double [] k4  = new double[3];
                
        //Find new Coordinate system which is normal to (n0,n1,n2)
        //Find a vector in n-plane
        for( i = 0; i < 3 && orbit == 0; i++)
        {     
            if(Math.abs(n[i]) > TINY)
            {                
                //Make all parallel planes equivalent by switching to unit n-plane that
                //has the first non-zero component which is positive 
                dtemp = 1.0;
                if(n[i] < TINY)
                    dtemp = -1;     
                //Make u a position vector ON the UNIT-n-plane
                u[(i+1)%3] = 1 +  dtemp*n[(i+1)%3]/vectorMag(n);
                u[(i+2)%3] = 1;
                u[i] =  vectorMag(n) - (u[(i+1)%3]*n[(i+1)%3] + u[(i+2)%3]*n[(i+2)%3])*dtemp;
                u[i] /= n[i]*dtemp;
                //Now make u a direction vector IN the n-plane
                for( j = 0; j<3; j++)
                    u[j] = u[j] - (dtemp*n[j]/vectorMag(n));
                //Make u a unit vector
                dtemp = vectorMag(u);
                for( j = 0; j<3; j++)
                    u[j] /= dtemp;
                //Make v perpendicular to u and n               
                v[0] = n[1]*u[2] - n[2]*u[1];
                v[1] = n[2]*u[0] - n[0]*u[2];
                v[2] = n[0]*u[1] - u[0]*n[1];
                //Make v a unit vector
                dtemp = vectorMag(v);
                //Make v lie one the most positive N-plane
                if(n[i] < TINY)
                    dtemp *= -1;
                for( j = 0; j<3; j++)
                    v[j] /= dtemp;
                 break;
             }
             else if ( vectorMag(n) < TINY )
             {
                //n will be freshed in the JPanel  - see Component
                n[0] = n[1] = n[2] = 1; 
                i = 0;
             }       
        }                  
        
        for (j = 0; j < 2; j++){
        if( j == 0 && (reDraw || orbit > 2) )
         continue;
        x0 = 0.1;
        y0 = 0;
        z0 = 0;           
        for (i = 0; i < iterations; i++) 
        {
            
            if(method_rk)
            {
                h_f_RungeKutta(x0, y0, z0,k1);
                h_f_RungeKutta(x0 + k1[0]/2, y0 + k1[1]/2, z0 + + k1[2]/2,k2);
                h_f_RungeKutta(x0 + k2[0]/2, y0 + k2[1]/2, z0 + + k2[2]/2,k3);
                h_f_RungeKutta(x0 + k3[0]  , y0 + k3[1]  , z0 + k3[2]    ,k4);
            
                x1 = x0 + k1[0]/6 + k2[0]/3 + k3[0]/3 + k4[0]/6;
                y1 = y0 + k1[1]/6 + k2[1]/3 + k3[1]/3 + k4[1]/6;
                z1 = z0 + k1[2]/6 + k2[2]/3 + k3[2]/3 + k4[2]/6;
            }
            else
            {
                x1 = x0 + h * a * (y0 - x0);
                y1 = y0 + h * (x0 * (b - z0) - y0);
                z1 = z0 + h * (x0 * y0 - c * z0);
            }            
            x0 = x1;
            y0 = y1;
            z0 = z1;
            
            
           if(Math.abs(x1-x0) < TINY && Math.abs(y1-y0) < TINY && Math.abs(z1-z0) < TINY ) 
           {
                if (Math.abs(x1-fixedPoint[0]) < TINY && Math.abs(y1-fixedPoint[1]) < TINY && Math.abs(z1-fixedPoint[2]) < TINY)
                    fixedPointsCount++;
                else
                {
                    fixedPoint[0] = x1; fixedPoint[1] = y1; fixedPoint[2] = z1;
                    fixedPointsCount = 0;
                }
           }
           
           ///////////////////////// Fixed Points and NaNs
           if(Double.isNaN(x0) )
           {    
                gContext.setColor ( Color.magenta );
                gContext.drawString("Lorenz Points have become unbounded",20,20);
                System.out.println("Lorenz Points have become unbounded");
                break;
           } 
           else if( fixedPointsCount > iterations/3)
           {
                gContext.setColor ( Color.magenta );
                gContext.drawString("Lorenze Points have converged to a fixed point - No Chaotic attractor to show",20,20);
                System.out.println("Lorenze Points have converged to a fixed point\n - No Chaotic attractor to show");
                break;
           }
                 
            //Map x0 y0 z0 to the n-plane using u and v            
            nx = (x1 - n[0])*u[0] + (y1 - n[1])*u[1] + (z1 - n[2])*u[2]; 
            ny = (x1 - n[0])*v[0] + (y1 - n[1])*v[1] + (z1 - n[2])*v[2]; 
                
            //Now Rotate!
            dtemp = nx*Math.cos(TWOPI*rotationAngle/360) - ny*Math.sin(TWOPI*rotationAngle/360);
            ny = nx*Math.sin(TWOPI*rotationAngle/360) + ny*Math.cos(TWOPI*rotationAngle/360);
            nx = dtemp;
             
            if(i > iterationLapse )
            {    
                if( nx > xmax )
                    xmax = nx;
                else if( nx < xmin )
                    xmin = nx;
                        
                if( ny > ymax )
                    ymax = ny;
                else if( ny < ymin )
                    ymin = ny;
                lorenzBoundingBox(x1, y1, z1);
            }
            else if( i == iterationLapse )
            {
                xmax = xmin = nx;
                xmax += TINY;
                ymax = ymin = ny;
                ymax += TINY;
                lorenzBoundingBox(x1, y1, z1);
            }
                
            //Cant be updated during paintprocess!
            //Cant change during orbit eaither!
            if(!rangesSet && ((j == 0 && i == iterations - 1) || j == 1) )
            {
                nymax = ymax;
                nymin = ymin;
                nxmax = xmax;
                nxmin = xmin;
                rangesSet = true;
            }
                
             if(orbit == 0 && iterations/4== i )
             {
                for( int m = 0; m < 3; m++)
                {
                    p[m] = (l_max[m] + l_min[m])/2;
                    l_min[m] = l_max[m] = 0.0;
                }
             }    
            
            //Now scale to fit the screen
            pxold = px;
            pyold = py;
                
            py = pHeight*(ny - nymin)/(nymax - nymin);
            py = pHeight - py + pHeight*vShift/10;
            px = pWidth * (nx - nxmin)/(nxmax - nxmin) + pWidth*hShift/10;
            
            //calculate where P is on the screen
            if( orbit > 0 )
            {
                double p_py = 0, p_px = 0, p_ny = 0, p_nx = 0;
                //Map p to the n-plane using u and v
                p_nx = (p[0] - n[0])*u[0] + (p[1] - n[1])*u[1] + (p[2] - n[2])*u[2]; 
                p_ny = (p[0] - n[0])*v[0] + (p[1] - n[1])*v[1] + (p[2] - n[2])*v[2];               
                //Now Rotate!
                dtemp = p_nx*Math.cos(TWOPI*rotationAngle/360) - p_ny*Math.sin(TWOPI*rotationAngle/360);
                p_ny = p_nx*Math.sin(TWOPI*rotationAngle/360) + p_ny*Math.cos(TWOPI*rotationAngle/360);
                p_nx = dtemp;
                //Now scale to fit the screen
                p_py = pHeight*( p_ny - nymin)/(nymax - nymin);
                p_py = pHeight - p_py + pHeight*vShift/10;
                p_px = pWidth*(p_nx - nxmin)/(nxmax - nxmin) + pWidth*hShift/10;
                
                //Now shift the other points going to the screen.
                py += pHeight/2 - p_py;
                px += pWidth/2 - p_px;                
            }
                
            //Zoom Correction
            px = (px-pWidth/2)*(100.0+dZoom)/100.0 + pWidth/2.0;
            py = (py-pHeight/2)*(100.0+dZoom)/100.0 + pHeight/2.0;
                
            if( j == 1 && i > iterationLapse && rangesSet)
            {       
                gContext.draw( new Line2D.Double( pxold, pyold, px, py ));
                switch ( (i/( (int)(iterations/40)))%7)
                {
                    case 0: gContext.setColor ( Color.red );break;
                    case 1: gContext.setColor ( Color.getHSBColor(1,1,1) );break;
                    case 2: gContext.setColor ( Color.black );break;
                    case 3: gContext.setColor ( Color.gray );break;
                    case 4: gContext.setColor ( Color.yellow );break;
                    case 5: gContext.setColor (  Color.blue );break;
                    case 6: gContext.setColor (  Color.magenta );break;
                    
                    default: gContext.setColor ( Color.getHSBColor(3,5,5));break;
                }
            }
            }
        }
    //This stops the nxmax/min nymax.min from being computed again
    reDraw = true;
    }
   
   public void actionPerformed(ActionEvent evt) 
    {
         double dtemp = 0;
         Object source = evt.getSource();
         if ( source == bUp)
         {
            vShift--; repaint();
         }   
         else if (source == bDown)
         {
            vShift++; repaint();
         }
         else if ( source == bLeft )
         {
            hShift--; repaint();
         }
         else if ( source == bRight )
         {
            hShift++; repaint();
         }
         else if ( source == bZoomIn )
         {       
            dZoom += 20; repaint();
         }
         else if ( source == bZoomOut )
         {
            dZoom -= 20; repaint();
         }
         else if( source == bDraw )
         {
           try
           {           
                n[0] = Double.parseDouble( n0Value.getText());
                n[1] = Double.parseDouble( n1Value.getText());
                n[2] = Double.parseDouble( n2Value.getText());
                a    = Double.parseDouble( aValue.getText());
                b    = Double.parseDouble( bValue.getText());
                c    = Double.parseDouble( cValue.getText());                
                rotationAngle = Double.parseDouble( rValue.getText());
                if( (dtemp = Double.parseDouble( iValue.getText())) > 100 )
                    iterations = (long)dtemp;
                else
                    iValue.setText(""+iterations);                                     
                dZoom = vShift = hShift = 0;
                ymin = ymax = xmin = xmax = 0;
                reDraw = false;
                repaint();
                
            }
            catch (NumberFormatException nfe)
            {
                System.out.println("\nCan not redraw with those parameter values");
            }                   
         }
         else if( source == bOrbit )
         {
            try
           {           
                n[0] = Double.parseDouble( n0Value.getText());
                n[1] = Double.parseDouble( n1Value.getText());
                n[2] = Double.parseDouble( n2Value.getText());
                a    = Double.parseDouble( aValue.getText());
                b    = Double.parseDouble( bValue.getText());
                c    = Double.parseDouble( cValue.getText());                
                rotationAngle = Double.parseDouble( rValue.getText());
                if( (dtemp = Double.parseDouble( iValue.getText())) > 100 )
                    iterations = (long)dtemp;
                else
                    iValue.setText(""+iterations);                                    
                reDraw = false;
                p[0] = Double.parseDouble( oP0Value.getText());
                p[1] = Double.parseDouble( oP1Value.getText());
                p[2] = Double.parseDouble( oP2Value.getText());
                v[0] = Double.parseDouble( oA0Value.getText());
                v[1] = Double.parseDouble( oA1Value.getText());
                v[2] = Double.parseDouble( oA2Value.getText());
                orbit = 1;
               
            }
            catch (NumberFormatException nfe)
            {
                System.out.println("\nCan not redraw with those parameter values");
            }                   
            
         }
         else if ( source == bReset)
         {
            
            n[0] = 10; n[1] = -3; n[2] = 1;
            iValue.setText(""+(iterations = 26000));
            rValue.setText("0.0");
            dZoom = vShift = hShift = 0;
            ymin = ymax = xmin = xmax = 0;
            reDraw = false;
            orbit = 0;
            aValue.setText(""+(a = 10.0)); bValue.setText(""+(b=28.0));cValue.setText(""+roundD(8.0/3.0,3));
            repaint();
         }
         else if ( source == bMethod )
         {
            if(method_rk)
            {
                method_rk = false;
                bMethod.setText("Euler");    
            }
            else
            {
                method_rk = true;
                bMethod.setText(" R-K ");
            }
            ymin = ymax = xmin = xmax = 0;
            reDraw = false;
            orbit = 0;
         } 
     }
     
  public void orbitButtonPress()
  {
     orbit++;
     //Make v a unit vector
     double mag =0, frames = 10;
     String myString = "";
     int i = 0;
     for(i = 0, mag = vectorMag(v); i< 1 ; i++)
        v[i]/=mag;
     frames = Math.max(10, 25 - iterations/2000);
     for(i = 0; i< (int)frames ; i++)
     {   
        myString = Double.toString((int)(90*i/frames));
        orbitAngleValue.setText(myString);
        orbitVector((TWOPI/(4*frames)));
        try 
        {
            Thread.sleep(400);
        } 
        catch (InterruptedException e){}
        reDraw = false;
     }
     orbitAngleValue.setText("0.0");
     orbit = 0;
   } 
   
   private void orbitVector(double rA)
   {
     double [] n_to_pRadialVector = new double [3];
     double dtemp, lambda = 0;
     //reduce n to orgin and then rotate (n-p) about v
     for(int i = 0; i< 3; i++)
        n[i] -= p[i];     
     rotateVectorAboutAxis(n, rA, v);
     for(int i = 0; i< 3; i++)
        n[i] += p[i];
     //calculate the radial vector from p-axis n-axis (both axes in v direction)
     lambda = vectorDot(v,n) - vectorDot(p,v);
     for(int i = 0; i< 3; i++)
       n_to_pRadialVector[i] = p[i] - n[i] + lambda*v[i];
     if(vectorMag(n_to_pRadialVector) < TINY)
     {
        System.out.println("Cannot Orbit this point as it lies on the v-axis through n");
        return;
     }
     vectorCrossProduct(v, n_to_pRadialVector, u);
     //Make u a unit vector
     dtemp = vectorMag(u);
     for(int i = 0; i< 3; i++)
        u[i]/= dtemp;
     paintImmediately(0,0,getSize().width, getSize().height);
     print_n_oA_oP();
   }
         
  public void paintComponent ( Graphics gr )
  { 
    //paintComponent already exists inside JPanel library, so 
    //by rewritting its defintion in this class we are overiding 
    //the original version
    //
    //paintComponent is always handed a Graphics object which must 
    //be typecast into Graphics2D form so that wider range of functions
    //can be used on it
  
    Graphics2D gContext = (Graphics2D) gr;
    
    //Improvement of graphical output, see SAM's page 356
    gContext.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
    RenderingHints.VALUE_ANTIALIAS_ON);
    
    //Set Drawing Stroke of the line, see page 365 in Sam's Java
    // arguments are thickness, cap (line end) and join ( style in which lines meet)
    BasicStroke pen = new BasicStroke( 0.3f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND);
    gContext.setStroke( pen );
    
    Dimension dimJPanel = getSize();
    
    //This ensures that the JPanel is completely white
    gContext.setColor ( Color.white  );
    gContext.fillRect( 0, 0, dimJPanel.width, dimJPanel.height);
    calculatePoints( gContext ); 
    if(orbit == 0)
        print_n_oA_oP();
        
  }
} 