package lab8;

import edu.mines.jtk.util.Parallel;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.util.*;


public class NearestNeighbor{

  private float _flag;


  public NearestNeighbor(float flag){
    _flag = flag;    
  }


  public void apply(float[] scatter, float[] interp){

    int n1 = scatter.length;
    // get closest point
    ClosestPointTransform cpt = new ClosestPointTransform(_flag);
    int[] i1o = new int[n1];
    float[] dt = new float[n1];
    cpt.apply(scatter,i1o,dt);
    nearest(scatter,i1o,interp);
  }


  public void apply(float[][] scatter, float[][] interp){

    int n2 = scatter.length;
    int n1 = scatter[0].length;

    // get closest point
    ClosestPointTransform cpt = new ClosestPointTransform(_flag);
    int[][] i1o = new int[n2][n1];
    int[][] i2o = new int[n2][n1];
    float[][] dt = new float[n2][n1];

    System.out.printf("calculating Closest point Transform (CPT) %n");
    System.out.printf("n1=%4d n2=%4d %n",n1,n2);
    Stopwatch sw = new Stopwatch();
    sw.start();
    cpt.apply(scatter,i1o,i2o,dt);
    System.out.printf("finished CPT in %g sec%n",sw.time());
    nearest(scatter,i1o,i2o,interp);
    System.out.printf("Nearest Neighbor in %g sec%n",sw.time());
  }


  public void apply(float[][][] scatter, float[][][] interp){

    int n3 = scatter.length;
    int n2 = scatter[0].length;
    int n1 = scatter[0][0].length;

    // get closest point
    ClosestPointTransform cpt = new ClosestPointTransform(_flag);
    int[][][] i1o = new int[n3][n2][n1];
    int[][][] i2o = new int[n3][n2][n1];
    int[][][] i3o = new int[n3][n2][n1];

    float[][][] dt = new float[n3][n2][n1];

    System.out.printf("calculating Closest point Transform (CPT) %n");
    System.out.printf("n1=%4d n2=%4d n3=%4d %n",n1,n2,n3);
    Stopwatch sw = new Stopwatch();
    sw.start();

    cpt.apply(scatter,i1o,i2o,i3o,dt);

    System.out.printf("finished CPT in %g sec%n",sw.time());

    nearest(scatter,i1o,i2o,i3o,interp);
    System.out.printf("Nearest Neighbor in %g sec%n",sw.time());

  }




  /**
   *
   *  Nearest neighbor 3d gridding:
   *
   *  @param sparse : array with sparse samples and null values
   *  @param i1csp : array with i1 index of closest coordinate  
   *  @param i2csp : array with i2 index of closest coordinate 
   *  @param i3csp : array with i3 index of closest coordinate 
   *  @param nninterp : array with nearest neighbor interpolation
   */
  
  private void nearest(
    final float[][][] sparse, final int[][][] i1csp, final int[][][] i2csp, 
    final int[][][] i3csp, final float[][][] nninterp)
  {
    final int n1 = i2csp[0][0].length;
    final int n2 = i2csp[0].length;
    final int n3 = i2csp.length;

    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2 ; ++i2){ 
        for (int i1=0; i1<n1 ; ++i1){
          int i33 = i3csp[i3][i2][i1];
          int i22 = i2csp[i3][i2][i1];
          int i11 = i1csp[i3][i2][i1];
          nninterp[i3][i2][i1] = sparse[i33][i22][i11];
        }
      }
    }
    });
  }


  /**
   *
   *  Nearest neighbor 2d gridding:
   *
   *  @param sparse : array with sparse samples and null values
   *  @param i1csp : array with i1 index of closest coordinate  
   *  @param i2csp : array with i2 index of closest coordinate 
   *  @param nninterp : array with nearest neighbor interpolation
   */
  
  private void nearest(
   final float[][] sparse, final int[][] i1csp, final int[][] i2csp, final float[][] nninterp)
  {
    final int n2 = i2csp.length;
    final int n1 = i2csp[0].length;

    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1 ; ++i1){
        int i22 = i2csp[i2][i1];
        int i11 = i1csp[i2][i1];
        nninterp[i2][i1] = sparse[i22][i11];
      }
    }
    });
  }

  /**
   *
   *  Nearest neighbor 1d gridding:
   *
   *  @param sparse : array with sparse samples and null values
   *  @param i1csp : array with i1 index of closest coordinate  
   *  @param nninterp : array with nearest neighbor interpolation
   */

  private void nearest (
    float[] sparse, int[] i1csp, float[] nninterp)
  {
    int n1 = i1csp.length;
    for (int i1=0; i1<n1 ; ++i1){
      int i11 = i1csp[i1];
      nninterp[i1] = sparse[i11];
    }
  }


}
