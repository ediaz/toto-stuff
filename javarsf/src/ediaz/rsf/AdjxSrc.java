

package ediaz.rsf;

import ediaz.lib.*;
import java.util.*;
import java.io.IOException;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import rsf.RSF;
import rsf.Input;
import rsf.Output;
import edu.mines.jtk.util.Parallel;




public class AdjxSrc{
	static { 
	  System.loadLibrary("jrsf");
	}

  class Coordinates2 {
    /* coordinate structure  for 2d*/
    public float x,z;
  }


	public static void main (String[] args){
		RSF par = new RSF(args);

    // I-O files
		Input  Fswfl = new Input("in");
		Input  Fcoor = new Input("cc");
		Input  Feicx = new Input("eicx");
    Output Fadj  = new Output("out");


    int src = par.getInt("source",1);
    
    if (src !=1 && src!=-1)
      src=1;
 
    System.err.printf("src=%d\n",src);

    /* time axis geometry */
    float o3 = Fswfl.getOrigin(3); 
    float d3 = Fswfl.getDelta(3);
    int   n3 = Fswfl.getN(3);

    /* x axis geometry */
    float o2 = Fswfl.getOrigin(2); 
    float d2 = Fswfl.getDelta(2);
    int   n2 = Fswfl.getN(2);

    /* z axis geometry */
    float o1 = Fswfl.getOrigin(1); 
    float d1 = Fswfl.getDelta(1);
    int   n1 = Fswfl.getN(1);

    /* coordinates point */
    int nr = Fcoor.getN(2);

    /* get extended images lags */
    int   nl = Feicx.getN(1);
    float ol = Feicx.getOrigin(1);
    float dl = Feicx.getDelta(1);


    Fadj.setN(1,n3); Fadj.setDelta(1,d3); Fadj.setOrigin(1,o3);
    Fadj.setN(2,nr); Fadj.setDelta(2,1.0f); Fadj.setOrigin(2,0.0f);
    Fadj.setN(3,1); Fadj.setDelta(3,d1); Fadj.setOrigin(3,o1);



    
    float[][] cc = new float[nr][2];
    Fcoor.read(cc);
    /* read coordinates, and expand into the Coordinates class */


    float[][] adjsrc = new float[nr][n3];
    float[][][] wfl = new float[n3][n2][n1];
    float[] eic = new float[nl];

    float minlx = ol;
    float maxlx = (nl-1)*dl +ol;
    Fswfl.read(wfl);

    int perc = nr/10;
    for (int ir=0; ir<nr; ++ir){

      if(ir%perc == 0)
        System.err.printf("%3.1f %s writing done\n",ir/perc*10.0f,"%"); 

      Feicx.read(eic);
      int ix = (int)(0.5f +(cc[ir][0] -o2)/d2) ;
      int iz = (int)(0.5f +(cc[ir][1] -o1)/d1) ;
      for (int it=0; it<n3; ++it){
        for (int il=0; il<nl; ++il){
          float lx = il*dl+ol;
          float lxx = cc[ir][0] -2*src*lx;
          int ilx = (int)(0.5f+ (lxx-o2)/d2) ;
          if(ilx >=0 && ilx <n2 && il+src >=0 && il+src<nl){
            adjsrc[ir][it] += eic[il]*wfl[it][ilx][iz];
          }
        }      
      }         
    }
    

    
    Fadj.write(adjsrc);








    Fswfl.close();
    Fcoor.close();
    Fadj.close();
  }
                        
  private static void expandCoor(float[] x,Coordinates2[] coor){
    int nr = x.length/2;
    
    for (int ir=0; ir<nr; ir+=2){
      coor[ir].x = x[ir];
      coor[ir].z = x[ir+1];
    }
  }

}

