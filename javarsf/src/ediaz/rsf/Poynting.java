

package ediaz.rsf;

import ediaz.lib.*;
import java.util.*;
import static java.lang.Math.sqrt;
import java.io.IOException;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import rsf.RSF;
import rsf.Input;
import rsf.Output;


public class Poynting
{
	static { 
	  System.loadLibrary("jrsf");
	}

	public static void main (String[] args){
		RSF par = new RSF(args);

    // I-O files
		Input input = new Input("in");
		Output output = new Output("out");

    /* time axis geometry */
    float o3 = input.getOrigin(3); 
    float d3 = input.getDelta(3);
    int   n3 = input.getN(3);

    /* x axis geometry */
    float o2 = input.getOrigin(2); 
    float d2 = input.getDelta(2);
    int   n2 = input.getN(2);

    /* z axis geometry */
    float o1 = input.getOrigin(1); 
    float d1 = input.getDelta(1);
    int   n1 = input.getN(1);

    output.setN(1,n1); output.setDelta(1,d1);   output.setOrigin(1,o1);
    output.setN(2,n2); output.setDelta(2,d2);   output.setOrigin(2,o2);
    output.setN(3,2) ; output.setDelta(3,1.0f); output.setOrigin(3,0.f);
    output.setN(4,n3); output.setDelta(4,d3);   output.setOrigin(4,o3);

    /* smoothing parameters from user */
    float sigma1 = par.getFloat("sigmax",1.f);
    float sigma3 = par.getFloat("sigmat",2.f);
    int torder = par.getInt("torder",6);
    int xorder = par.getInt("xorder",4);

    /* allocate working arrays */
    float[][][] tder  = zerofloat(n1,n2,n3);
    float[][][] tder2 = zerofloat(n1,n2,n3);
    float[][]   xder  = zerofloat(n1,n2);
    float[][]   zder  = zerofloat(n1,n2);
    float[][]   wfl2d = zerofloat(n1,n2);
    

    ExpSmooth xsm = new ExpSmooth(sigma1);
    ExpSmooth tsm = new ExpSmooth(sigma3);

    input.read(tder);
    copy(tder,tder2);

    System.err.printf("read data: done\n");

    Deriv dt = new Deriv(d3,torder);
    Deriv dx = new Deriv(d2,xorder);
    Deriv dz = new Deriv(d1,xorder);

    dt.apply33(tder,tder);
    System.err.printf("du/dt: done\n");
  

    for (int i3=0; i3<n3; ++i3){
      dz.apply21(tder2[i3],zder);
      dx.apply22(tder2[i3],xder);
      mul(zder,tder[i3],tder2[i3]); // tder2 stores z component pv
      mul(xder,tder[i3],tder[i3]);  // tder stores x component pv
    }
    System.err.printf("du/dt*grad(u): done\n");


    if (sigma1 >0.0f ){
      xsm.apply(tder,tder);
      xsm.apply(tder2,tder2);
    }
    if (sigma3 >0.0f){
      tsm.apply3(tder,tder);
      tsm.apply3(tder2,tder2);
    }
    System.err.printf("smoothing : done\n");

    int perc = n3/10;
    for (int i3=0; i3<n3; ++i3){
      if(i3%perc == 0)
        System.err.printf("%3.1f %s writing done\n",i3/perc*10.0f,"%"); 
      output.write(tder[i3]);
      output.write(tder2[i3]);
    }
    output.close();
    input.close();

  }
}

