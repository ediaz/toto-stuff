

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




public class AdjxSrc2{
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
		Output Fcout = new Output("ccout");
		Input  Feicx = new Input("eicx");
    Output Fadj  = new Output("out");


    int src = par.getInt("source",1);
    
    if (src !=1 && src!=-1)
      src=1;
    boolean source = false;
    if (src==1) 
      source = true;

 
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



 
    float[][] cc = new float[nr][2];
    Fcoor.read(cc);
    /* read coordinates, and expand into the Coordinates class */


    float[][] adjsrc = new float[nr][n3];
    float[][][] wfl = new float[n3][n2][n1];
    float[] eic = new float[nl];

    Fswfl.read(wfl);

    int perc = nr/10;
    int nlx = (nl-1)/2;


    // possible adjoint source location values:
    float[][] grid = new float[n2][n1];
    int live =0;
    for (int ir=0; ir<nr ; ++ir){

      int ix = (int)(0.5f +(cc[ir][0] -o2)/d2);
      int iz = (int)(0.5f +(cc[ir][1] -o1)/d1);
      
      int minx1 = Math.min(ix,nlx);
      int minx2 = Math.min(n2-1-ix,nlx);
      int mlx = Math.min(minx1,minx2);

      for (int ilx=-mlx; ilx<mlx+1; ++ilx){
        grid[ix+ilx][iz] =1.f;
      }
    }
    
    live =0;
    for (int i2=0; i2<n2; ++i2){
      for (int i1=0; i1<n1; ++i1){
        if(grid[i2][i1] == 1.f){
          live+=1;
        }
      }
    }

    System.err.printf("live=%d\n",live);
    float[][] ccout = new float[live][2];
    int[]   ccx   = new int[live];
    int[]   ccz   = new int[live];
    live=0;
    for (int i2=0; i2<n2; ++i2){
      for (int i1=0; i1<n1; ++i1){
        if(grid[i2][i1] == 1.f){
          ccout[live][0] = i2*d2+o2;
          ccout[live][1] = i1*d1+o1;
          ccx[live] = i2;
          ccz[live] = i1;
          live+=1;
        }
      }
    }
    
    float[][][] adjsrc2 = new float[n3][n2][n1];

    for (int ir=0; ir<nr; ++ir){
      if(ir%perc == 0)
        System.err.printf("%3.1f %s writing done\n",ir/perc*10.0f,"%"); 

      Feicx.read(eic);
      int ix = (int)(0.5f +(cc[ir][0] -o2)/d2);
      int iz = (int)(0.5f +(cc[ir][1] -o1)/d1);
      int minx1 = Math.min(ix,nlx);
      int minx2 = Math.min(n2-1-ix,nlx);
      int mlx = Math.min(minx1,minx2);

//      parallelLoop(mlx,source,adjsrc2,eic,wfl,ix,iz);

      for (int il=-mlx; il<mlx+1; ++il){
        int imx = ix-il;
        int ipx = ix+il;
        if(source){
          for (int it=0; it<n3; ++it){
            adjsrc2[it][ipx][iz] += eic[nlx+il]*wfl[it][imx][iz];
          }
        }else{
          for (int it=0; it<n3; ++it){
            adjsrc2[it][imx][iz] += eic[nlx+il]*wfl[it][ipx][iz];
          }
        }
      } 
            
    }

    float[][] adjout = new float[live][n3];
    
    for (int i=0; i<live; ++i){
      int ix = ccx[i];
      int iz = ccz[i];
      for (int it=0; it<n3;++it)
        adjout[i][it] = adjsrc2[it][ix][iz]; 
    }
    adjsrc2=null;

    oaxa(Fadj, n3,     d3, o3,   1);
    oaxa(Fadj, live, 1.0f, 0.0f, 2);
    oaxa(Fadj, 1,      d1, o1,   3);

    oaxa(Fcout,2,1.0f,0.0f,1);
    oaxa(Fcout,live,1.0f,0.0f,2);
    oaxa(Fcout,1,1.0f,0.0f,3);

    



    System.err.printf("writing starts\n"); 

    Fcout.write(ccout); 
    Fadj.write(adjout);
    System.err.printf("closing files \n"); 

    Fswfl.close();
    Fcoor.close();
    Fadj.close();
    Fcout.close();
  }
                       


  private static void parallelLoop(final int mlx,final boolean source,
                            final float[][][] adjsrc2,
                            final float[] eic,
                            final float[][][] wfl,
                            final int ix,
                            final int iz){

    final int nlx = (eic.length-1)/2;
    final int n3 = wfl.length;
    


    Parallel.loop(mlx*2+1,new Parallel.LoopInt() {   
    public void compute(int ill) {
      int il = ill-mlx;
//    for (int il=-mlx; il<mlx+1; ++il){
      int imx = ix-il;
      int ipx = ix+il;
      if(source){
        for (int it=0; it<n3; ++it){
          adjsrc2[it][ipx][iz] += eic[nlx+il]*wfl[it][imx][iz];
        }
      }else{
        for (int it=0; it<n3; ++it){
          adjsrc2[it][imx][iz] += eic[nlx+il]*wfl[it][ipx][iz];
        }
      }
    }}); 

  }

  private static void oaxa(Output File,int n, float d, float o, int axis){
    File.setN(axis,n);
    File.setOrigin(axis,o);
    File.setDelta(axis,d);
  }



}

