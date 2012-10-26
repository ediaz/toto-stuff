package lab4;

import java.util.concurrent.atomic.AtomicInteger;
import edu.mines.jtk.util.Parallel;


/**
 * Digital signal processing.
 * @author Esteban Diaz (modified from Dave Hale's code),
 *         Colorado School of Mines
 * @version 2012.10.04
 *




  Fork-join task splitting:
  
  The number of tasks is splitted by the parallel class
  until too many tasks are created, or the remainder
  range for indeces is smaller than chunksize, however
  defaul is chunksize=1, so the first criterion is 
  the most likely to prevail.

  Once splitted, each task has a similar or equal size (
  so the work is as balanced as possible)




*/
public class Dsp {

  /**
   * Number of threads used for parallel processing.
   */
  public static int nthread = Runtime.getRuntime().availableProcessors();

  /**
   * Smooths a 1D array.
   */
  public static void smooth(float a, float[] x, float[] y) {
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
   * Smooths a 2D array.
   */
  public static void smooth(float a, float[][] x, float[][] y) {
    
    smooth1P(a,x,y);
//    smooth2P(a,x,y);

  }


  /**
   * Smooths a 2D array along only the 1st dimension.
   * Serial version.
   */
  public static void smooth1(
    final float a, final float[][] x, final float[][] y) 
  {
    int n = x.length;
    for (int i=0; i<n; ++i)
      smooth(a,x[i],y[i]);
  }

  /**
   * Smooths a 2D array along only the 1st dimension.
   * Parallel version.
   */
  public static void smooth1P(
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
  public static void smooth2(float a, float[][] x, float[][] y) {
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
  public static void smooth2P(final float a,final float[][] x,final float[][] y) {
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
  




  public static float mean(float[] x) {
    int n = x.length;
    float sum = 0.0f;
    for (int i=0; i<n; ++i)
      sum += x[i];
    return sum/n;
  }

  public static float mean(float[][] x) {
    int n = x.length;
    float sum = 0.0f;
    for (int i=0; i<n; ++i)
      sum += mean(x[i]);
    return sum/n;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  /** 
   * Smooth along second dimension in chunks along first dimension.
   * 
   */

  private static void smooth2_chunk(
    final float a,final float[][] x, final float[][] y, final int init, final int chunk){


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
   * Based on the number of threads fills two arrays 
   * of integers with lower and upper bounds.
   * Then, you can use this bounds to split the 
   * work along a given axis.
   */
  private static void getchunks (int n,int []i1, int[] i2) {
    int chunk_size = (int) (n/nthread) ;

    i1[0] = 0;   
    i2[0] = i1[0] + chunk_size;
    for (int i=1 ; i<nthread; ++i){
      i1[i] = i2[i-1] + 1;
      i2[i] = i1[i] + chunk_size;
    }
    i2[nthread-1] = n-1;
  
  }

  /**
   * Parallel transpose
   */

  private static void transposeP(final float[][] x, final float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    final AtomicInteger ai = new AtomicInteger();
    Thread[] threads = new Thread[nthread];

    final int[] lb = new int[nthread]; 
    final int[] ub = new int[nthread];

    // Calculate bounds of the chunk:
    getchunks(n2,lb,ub);


    for (int ithread=0; ithread<nthread; ++ithread) {
      threads[ithread] = new Thread() {
        public void run() {
          for (int i=ai.getAndIncrement(); i<nthread; i=ai.getAndIncrement()) {
            transpose_chunk( x, y,lb[i],ub[i]);
          }
        }
      };
    }
    startAndJoin(threads);
  }


  /** 
   * Serial transpose for one piece of the array x
   */

  private static void transpose_chunk(float[][] x, float[][] y,int lb, int ub) {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=lb ; i2<ub+1; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        y[i1][i2] = x[i2][i1];
      }
    }
  }



  private static void startAndJoin(Thread[] threads) {
    int nthread = threads.length;
    for (int ithread=0; ithread<nthread; ++ithread) {
      threads[ithread].start();
    }
    for (int ithread=0; ithread<nthread; ++ithread) {
      try {
        threads[ithread].join();
      } catch (InterruptedException e) {
        throw new RuntimeException(e);
      }
    }
  }
}
