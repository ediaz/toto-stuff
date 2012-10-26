import edu.mines.jtk.util.Parallel;
/*
*   transpose using the outher loop 
    to be the faster of the destination 
    array

    for i1 
      for i2
        y[i1][i2] = x[i2][i1]

* 
*/ 

/**
 * Array copy and transposing
 * @author Esteban Diaz
 *         Colorado School of Mines
 * @version 2012.10.24
 *
*/
public class Copy {

  /**
   * Number of threads used for parallel processing.
   */
  public static int nthread = Runtime.getRuntime().availableProcessors();

  /**
   * Copy a 1D array.
   */
  public static void copy1S(float[] x, float[] y) {
    int n1 = x.length;

    for (int i1=0 ;i1<n1 ; ++i1)
      y[i1] = x[i1] ; 
  }


  /**
   * Copy a 2D array.
   */
  public static void copy1S(float[][] x, float[][] y) {
    int n2 = x.length;

    for (int i2=0 ;i2<n2 ; ++i2)
      copy1S(x[i2],y[i2]) ; 
  }

  /**
   * Copy a 3D array.
   */
  public static void copy1S(float[][][] x, float[][][] y) {
    int n3 = x.length;

    for (int i3=0 ;i3<n3 ; ++i3)
      copy1S(x[i3],y[i3]) ; 
  }


  /**
   * Copy a 3D array, parallelizing over outher
   * dimension
   */  
  public static void copy1P(final float[][][] x,final float[][][] y) {

    int n3 = x.length ;

    Parallel.loop(n3,new Parallel.LoopInt() {
       public void compute(int i3) {
         copy1S(x[i3],y[i3]);
       }
    });

  }


  /**
   * Copy a 3D array, using "loop unrolling"
   * over the last dimension, with block size=8
   */
  public static void copy1S8(final float [][][] x,final float [][][] y) {
  
    final int n3 = x.length;
    int block_size = 8 ;

    int remainder =  n3%block_size ; 

    for (int i3=0; i3< remainder ; ++i3)
      copy1S(x[i3],y[i3]) ; 
    

    for (int i3=remainder ;i3<n3 ; i3+= block_size)
      copyblock8(x,y,i3); 
  }

  /**
   * Copy a 3D array, using "loop unrolling"
   * over the last dimension, with block size=8,
   * parallelizing over the last dimension
   */

  public static void copy1P8(final float [][][] x,final float [][][] y) {
  
    int n3 = x.length;
    int block_size = 8 ;

    for (int i3=0; i3< n3%block_size ; ++i3){
      copy1S(x[i3],y[i3]) ; 
    }

    Parallel.loop((int) n3/block_size,new Parallel.LoopInt() {
      public void compute(int i3) {
        copyblock8(x,y,i3) ;
      }
    });
    
  }







  /////////////////////////////////////////////////////////////////////////
  // private

  /**
   * Copy a 3D array, using "loop unrolling"
   * over the last dimension, with block size=8
   */
  private static void copyblock8(
    final float [][][] x,final float [][][] y,final int i3) {
    copy1S(x[i3  ],y[i3  ]) ; 
    copy1S(x[i3+1],y[i3+1]) ; 
    copy1S(x[i3+2],y[i3+2]) ; 
    copy1S(x[i3+3],y[i3+3]) ; 
    copy1S(x[i3+4],y[i3+4]) ; 
    copy1S(x[i3+5],y[i3+5]) ; 
    copy1S(x[i3+6],y[i3+6]) ; 
    copy1S(x[i3+7],y[i3+7]) ; 
  }




}
