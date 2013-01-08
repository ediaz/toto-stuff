package lab6;

import edu.mines.jtk.util.Parallel;



/**
 * Digital signal processing.
 * @author Esteban Diaz
 * @version 2012.11.04
 *


*/
public class ExpSmooth{

public static int nthread = Runtime.getRuntime().availableProcessors();

  /*
    Common variables:
  */


  private float _sigma;
  private int _nsmooth ; 
  private float _a ; 

  /*
    Constructor for exponential smoothing
    sigma is an approximation to the standard
    deviation of a gaussian filter.
  */
  public ExpSmooth(float sigma){
    this(sigma,8); // 8 smooths is the default
  }

  public ExpSmooth(float sigma,int nsmooth){
    _sigma = sigma;
    _nsmooth = nsmooth; 

    float epsq = (sigma*sigma)/nsmooth;
    _a= (float)((epsq+1.0f-Math.sqrt(2.0f*epsq+1.0f))/epsq);
  }

  /*
  TODO: create a new constructor that 
        takes a boolean variable to set
        up the boundary condition (zero or zero slope).
  */

  // 1D smoothing
  public void apply(float[] x, float[] y){
    smooth(_a,x,y);
    for (int nr=1 ; nr <_nsmooth ; ++nr) 
      smooth(_a,y,y);
  } 


  // 2D smoothing
  public void apply(float[][] x,float[][] y){
    smooth1P(_a,x,y);
    for (int nr=1 ; nr <_nsmooth ; ++nr) 
      smooth1P(_a,y,y);

    smooth2P(_a,y,y);
    for (int nr=1 ; nr <_nsmooth ; ++nr) 
      smooth2P(_a,y,y);
  }



  ////////////////////////////////////////////////////////////////////////
  // private

  /** 
   * Smooth along second dimension in chunks along first dimension.
   * 
   */

  private static void smooth2_chunk(
    final float a,final float[][] x, final float[][] y,
     final int init, final int chunk){


    int n2 = x.length;
    int n1 = x[0].length;
    int lb = init;
    int ub = lb +chunk -1 ;

    ub = Math.min (ub,n1-1) ;

    float b = 1.0f-a;

    for (int i1=lb;i1<= ub  ; ++i1)
      y[0][i1] = x[0][i1]; //first
 
    for (int i2=1;i2<n2-1; ++i2)
      for (int i1=lb;i1<= ub ; ++i1)
        y[i2][i1] = a*y[i2-1][i1]+b*x[i2][i1]; //forward
   
    for (int i1=lb;i1<= ub ; ++i1)
      y[n2-1][i1] = (a*y[n2-2][i1]+x[n2-1][i1])/(1.0f+a); //last
    
    for (int i2=n2-2;i2>=0; --i2)
      for (int i1=lb;i1<= ub  ; ++i1)
        y[i2][i1] = a*y[i2+1][i1]+b*y[i2][i1]; //reverse 
  }



  /**
   * Smooths a 1D array.
   */
  private static void smooth(float a, float[] x, float[] y) {
    int n = x.length;
    float b = 1.0f-a;
    float yi = y[0] = x[0]; // first
    for (int i=1; i<n-1; ++i) // forward
      y[i] = yi = a*yi+b*x[i];
    y[n-1] = yi = (a*yi+x[n-1])/(1.0f+a); // last
    for (int i=n-2; i>=0; --i) // reverse
      y[i] = yi = a*yi+b*y[i];
  }


  /**
   * Smooths a 2D array along only the 1st dimension.
   * Parallel version.
   */
  private static void smooth1P(
    final float a, final float[][] x, final float[][] y) 
  {
    final int n = x.length;

    Parallel.loop(n,new Parallel.LoopInt() {
       public void compute(int i) {
         smooth(a,x[i],y[i]);
       }
    });
  }


  /**
   * Smooths along the 2nd dimension, with outer loops over the
   * the 2nd dimension and inner loops over the 1st dimension.
   */
  private static void smooth2(float a, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;

    float b = 1.0f-a;

    for (int i1=0;i1<n1; ++i1)
      y[0][i1] = x[0][i1]; //first 

    for (int i2=1;i2<n2-1; ++i2)
      for (int i1=0;i1<n1;++i1)
        y[i2][i1] = a*y[i2-1][i1]+b*x[i2][i1]; //forward

    for (int i1=0;i1<n1; ++i1)
      y[n2-1][i1] = (a*y[n2-2][i1]+x[n2-1][i1])/(1.0f+a); //last
    
    for (int i2=n2-2;i2>=0; --i2)
      for (int i1=0;i1<n1;++i1)
        y[i2][i1] = a*y[i2+1][i1]+b*y[i2][i1]; //reverse
  }


  /**
   * Smooths along the 2nd dimension, with outer loops over the
   * the 2nd dimension and inner loops over the 1st dimension.
   */
  private static void smooth2P(
    final float a,final float[][] x,final float[][] y) {

    final int   n2 = x.length;
    final int   n1 = x[0].length;
    final float  b = 1.0f-a;
    final int chunk_size =(int) n1/nthread ;

    Parallel.loop(0,n1,chunk_size,new Parallel.LoopInt() {
       public void compute(int i1) {
         smooth2_chunk(a,x,y,i1,chunk_size);
       }
     });
    
  } 
  


}
