package lab0;

/**
 * Digital signal processing.
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.08.19
 */
public class Dsp {
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
  public static float mean(float[] x) {
    int n = x.length;
    float sum = 0.0f;
    for (int i=0; i<n; ++i)
      sum += x[i];
    return sum/n;
  }
}
