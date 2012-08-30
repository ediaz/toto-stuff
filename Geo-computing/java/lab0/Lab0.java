package lab0;

/**
 * Benchmark recursive smoothing filter in Java.
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.08.19
 */
public class Lab0 {
  public static void main(String[] args) {
    System.out.println("Java");
    int n = 1001;
    double maxtime = 2.0;
    float a = 0.99f;
    float[] x = new float[n];
    float[] y = new float[n];
    x[0] = x[n/2] = x[n-1] = 1.0f;
    int nsmooth;
    Stopwatch sw = new Stopwatch();
    sw.start();
    for (nsmooth=0; sw.time()<maxtime; ++nsmooth)
      Dsp.smooth(a,x,y);
    sw.stop();
    System.out.println("nsmooth = "+nsmooth);
    System.out.println("   mean = "+Dsp.mean(y));
    System.out.println("   time = "+sw.time());
    System.out.println(" mflops = "+(int)(6.0e-6*n*nsmooth/sw.time()));
  }
}
