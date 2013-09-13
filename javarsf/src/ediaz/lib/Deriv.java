package ediaz.lib; 

import edu.mines.jtk.util.Parallel;


/**
 * Digital signal processing.
 * @author Esteban Diaz
 * @version 2012.11.04
 *


*/
public class Deriv{

public static int nthread = Runtime.getRuntime().availableProcessors();

  /*
    Common variables:
  */

  private static float _d;
  private static int _n;

  /*
    Constructor for exponential smoothing
    sigma is an approximation to the standard
    deviation of a gaussian filter.
  */
  public Deriv(float d){
    this(d,4); // 8 smooths is the default
  }
  public Deriv(float d, int order){
    _d = 1.0f/d;
    _n = order;
  }


  public void apply22(final float[][] x, final float[][]y){
    final int n2 = x.length;
    final int n1 = x[0].length;


    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
    
      float[] xx = new float[n2];
      float[] yy = new float[n2];
      for (int i2=0; i2<n2;++i2){
        xx[i2] = x[i2][i1]; 
      }
      sf_deriv(xx,yy,_n);
      for (int i2=0; i2<n2;++i2){
        y[i2][i1] = yy[i2]*_d; 
      }
    }
    });

  }


  public void apply21(final float[][] x,final float[][]y){
    final int n2 = x.length;
    final int n1 = x[0].length;


    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      sf_deriv(x[i2],y[i2],_n);
      for (int i1=0; i1<n1;++i1){
        y[i2][i1] *= _d; 
      }
    }});

  }


  public void apply33(final float[][][] x,final float[][][] y){
    final int n3 = x.length;
    final int n2 = x[0].length;
    final int n1 = x[0][0].length;

    Parallel.loop(n2,new Parallel.LoopInt() {
       public void compute(int i2) {
      float[] xx = new float[n3];
      float[] yy = new float[n3];
      for (int i1=0; i1<n1; ++i1){
        for (int i3=0; i3<n3; ++i3){
          xx[i3] = x[i3][i2][i1];
        }
        sf_deriv(xx,yy,_n);
        for (int i3=0; i3<n3; ++i3){
          y[i3][i2][i1] = yy[i3]*_d;
        }
      }
       }
    });
  }


  /********************************************************************/

  private static void sf_deriv (float[] trace, float[] trace2, int n)
  /*< derivative operator >*/
   {
    int nt = trace.length;
    float[] h = new float[nt];
    int i, it;
    
    for (it=0; it < nt; it++) {
    	h[it] = trace[it];
    } 

    for (i=n; i >= 1; i--) {
	    for (it=1; it < nt-1; it++) {
	      trace2[it] = h[it]-0.5f*(h[it+1]+h[it-1]);
    	}
	    trace2[0] = trace2[1];
  	  trace2[nt-1] = trace2[nt-2];

	    for (it=0; it < nt; it++) {
	      h[it] = trace[it] + trace2[it]*i/(2*i+1);
    	}
    }

    trace2[0] = h[1]-h[0];
    for (it=1; it < nt-1; it++) {
    	trace2[it] = 0.5f*(h[it+1]-h[it-1]);
    }
    trace2[nt-1] = h[nt-1]-h[nt-2];
  }



}
