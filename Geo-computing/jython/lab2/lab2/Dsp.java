package lab2;

/**
 * Digital signal processing.
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.09.09
 */
public class Dsp {

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
   * Smooths a 2D array along both 1st and 2nd dimensions.
   */
  public static void smooth(float a, float[][] x, float[][] y) {
    smooth1(a,x,y);
    smooth2(a,y,y);
  }

  /**
   * Smooths a 2D array along only the 1st dimension.
   */
  public static void smooth1(float a, float[][] x, float[][] y) {
    int n = x.length;
    for (int i=0; i<n; ++i) {
      smooth(a,x[i],y[i]);
    }
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
   * Smooths along the 2nd dimension, with a simple outer loop over
   * the 1st dimension and inner loops over the 2nd dimension.
   */
  public static void smooth2S(float a, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    float b = 1.0f-a;

    for(int i1=0; i1<n1 ; ++i1){
      float yi = y[0][i1] = x[0][i1]; // first

      for (int i2=1; i2<n2-1; ++i2) 
        y[i2][i1] = yi = a*yi+b*x[i2][i1]; // forward
  
      y[n2-1][i1] = yi = (a*yi+x[n2-1][i1])/(1.0f+a); // last

      for (int i2=n2-2; i2>=0; --i2) // reverse
        y[i2][i1] = yi = a*yi+b*y[i2][i1];     
    }

  }

  /**
   * Smooths along the 2nd dimension by (1) transposing the 1st and
   * 2nd dimensions, (2) smoothing along the 1st dimension, and (3) 
   * again transposing the 1st and 2nd dimensions. 
   */
  public static void smooth2T(float a, float[][] x, float[][] y) {
    int n2=x.length;
    int n1=x[0].length;
    float[][] y1=new float[n1][n2]; // intermediate array to
                                    // keep the transpose
    transpose(x,y1);
    smooth1(a,y1,y1);
    transpose(y1,y);
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

  private static void transpose(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        y[i1][i2] = x[i2][i1];
      }
    }
  }
}
