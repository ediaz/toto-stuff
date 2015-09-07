package lab8;

import edu.mines.jtk.util.Parallel;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;
import java.util.Vector;
import java.util.Iterator;

public class SibsonInterpolation{

  private float _flag;
  private int nthread = Runtime.getRuntime().availableProcessors();









  public SibsonInterpolation(float flag){
    _flag = flag;
  }


  public void apply(float[] scatter, float[] c ){
    int n1 = scatter.length;

    // get closest point
    ClosestPointTransform cpt = new ClosestPointTransform(_flag);
    int[] i1o = new int[n1];
    float[] dt = new float[n1];


    cpt.apply(scatter,i1o,dt);
    // accumulated hits vector:
    int[] n = new int[n1];

    for (int i1=0; i1<n1 ; ++i1){
      int i1cp = i1o[i1];
      float r = dt[i1]; 
      float f = scatter[i1cp];

      int i1min = max(i1 -(int)r,0); 
      int i1max = min(i1 +(int)r,n1-1); 

      for(int i11=i1min; i11<=i1max ; ++i11 ){
        float d = abs(i11-i1);
        if(d<=r){
          c[i11] += f;      
          n[i11]++;
        }
      }
    }

    for (int i1=0; i1<n1 ; ++i1){
      c[i1] = c[i1]/n[i1];
    }

  }



  // this implements a "scatter" approach for the Sibson's interpolant.



  public void applyP(final float[][] scatter, final float[][] c ){
    final int n1 = scatter[0].length;
    final int n2 = scatter.length;

    // get closest point
    ClosestPointTransform cpt = new ClosestPointTransform(_flag);
    final int[][] i1o = new int[n2][n1];
    final int[][] i2o = new int[n2][n1];
    final float[][] dt = new float[n2][n1];

    System.out.printf("===========================================%n");
    System.out.printf("calculating Closest point Transform (CPT) %n");
    System.out.printf("n1=%4d n2=%4d %n",n1,n2);
    Stopwatch sw = new Stopwatch();
    sw.start();

    cpt.apply(scatter,i1o,i2o,dt);

    System.out.printf("finished CPT in %g sec%n",sw.time());


    // accumulated hits vector:
    final int chunk =(int)n2/nthread;

    final float[][][] e = Parallel.reduce(nthread, new Parallel.ReduceInt<float[][][]>(){
      public float[][][] compute(int ithread){
      float[][][] e = new float[2][n2][n1];

      int lb = ithread*chunk;
      int ub = min((ithread+1)*chunk,n2);
      if (ithread == nthread -1)
        ub =n2;

     // System.out.printf(" ithread =%d  for( int i2=%d;  i2<%d; ++i2)%n",ithread,lb,ub); 

        for (int i2=lb; i2<ub; ++i2){ 
          for (int i1=0; i1<n1 ; ++i1){
            int i1cp = i1o[i2][i1];
            int i2cp = i2o[i2][i1];
            float r = dt[i2][i1]; 
            float f = scatter[i2cp][i1cp];
    
            int i1min = max(i1 -(int)r,0); 
            int i1max = min(i1 +(int)r,n1-1); 
            int i2min = max(i2 -(int)r,0); 
            int i2max = min(i2 +(int)r,n2-1); 
    
            for(int i22=i2min; i22<=i2max ; ++i22 ){
              for(int i11=i1min; i11<=i1max ; ++i11 ){
                float d = sqrt((i11-i1)*(i11-i1)+(i22-i2)*(i22-i2));
                if(d<=r){
                  e[0][i22][i11] += f;      
                  e[1][i22][i11]++;
                }
              }
            }
          }
        }
//         System.out.printf("i2 =%d%n",i2);
        return e;
      }

      public float[][][] combine(float[][][] ea, float[][][] eb){
        return add(ea,eb);
      }
    });
  
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
    
      for (int i1=0; i1<n1 ; ++i1){
        c[i2][i1] = e[0][i2][i1]/e[1][i2][i1];
      }
    }});

    System.out.printf("Sibson in %g sec%n",sw.time());
  }


  public void apply(float[][] scatter, float[][] c ){
    int n1 = scatter[0].length;
    int n2 = scatter.length;

    // get closest point
    ClosestPointTransform cpt = new ClosestPointTransform(_flag);
    int[][] i1o = new int[n2][n1];
    int[][] i2o = new int[n2][n1];
    float[][] dt = new float[n2][n1];

    System.out.printf("===========================================%n");
    System.out.printf("calculating Closest point Transform (CPT) %n");
    System.out.printf("n1=%4d n2=%4d %n",n1,n2);
    Stopwatch sw = new Stopwatch();
    sw.start();

    cpt.apply(scatter,i1o,i2o,dt);

    System.out.printf("finished CPT in %g sec%n",sw.time());


    // accumulated hits vector:
    int[][] n = new int[n2][n1];

    for (int i2=0; i2<n2 ; ++i2){
      for (int i1=0; i1<n1 ; ++i1){
        int i1cp = i1o[i2][i1];
        int i2cp = i2o[i2][i1];
        float r = dt[i2][i1]; 
        float f = scatter[i2cp][i1cp];

        int i1min = max(i1 -(int)r,0); 
        int i1max = min(i1 +(int)r,n1-1); 
        int i2min = max(i2 -(int)r,0); 
        int i2max = min(i2 +(int)r,n2-1); 

        for(int i22=i2min; i22<=i2max ; ++i22 ){
          for(int i11=i1min; i11<=i1max ; ++i11 ){
            float d = sqrt((i11-i1)*(i11-i1)+(i22-i2)*(i22-i2));
            if(d<=r ){
              c[i22][i11] += f;      
              n[i22][i11]++;
            }
          }
        }
      }
    }

    for (int i2=0; i2<n2 ; ++i2){
      for (int i1=0; i1<n1 ; ++i1){
        c[i2][i1] = c[i2][i1]/n[i2][i1];
      }
    }
    System.out.printf("Sibson in %g sec%n",sw.time());
  }
  public void apply(float[][][] scatter, float[][][] c ){
    int n1 = scatter[0][0].length;
    int n2 = scatter[0].length;
    int n3 = scatter.length;

    // get closest point
    ClosestPointTransform cpt = new ClosestPointTransform(_flag);
    int[][][] i1o = new int[n3][n2][n1];
    int[][][] i2o = new int[n3][n2][n1];
    int[][][] i3o = new int[n3][n2][n1];
    float[][][] dt = new float[n3][n2][n1];

    System.out.printf("===========================================%n");
    System.out.printf("calculating Closest point Transform (CPT) %n");
    System.out.printf("n1=%4d n2=%4d %n",n1,n2);
    Stopwatch sw = new Stopwatch();
    sw.start();

    cpt.apply(scatter,i1o,i2o,i3o,dt);

    System.out.printf("finished CPT in %g sec%n",sw.time());




    // accumulated hits vector:
    int[][][] n = new int[n3][n2][n1];

    for (int i3=0; i3<n3 ; ++i3){
      System.out.printf("i3=%d%n",i3);

      for (int i2=0; i2<n2 ; ++i2){
        for (int i1=0; i1<n1 ; ++i1){
          int i1cp = i1o[i3][i2][i1];
          int i2cp = i2o[i3][i2][i1];
          int i3cp = i3o[i3][i2][i1];
          float r = dt[i3][i2][i1]; 
          float f = scatter[i3cp][i2cp][i1cp];
  
          int i1min = max(i1 -(int)r,0); 
          int i1max = min(i1 +(int)r,n1-1); 

          int i2min = max(i2 -(int)r,0); 
          int i2max = min(i2 +(int)r,n2-1); 

          int i3min = max(i3 -(int)r,0); 
          int i3max = min(i3 +(int)r,n3-1); 
  
          for(int i33=i3min; i33<=i3max ; ++i33 ){
            for(int i22=i2min; i22<=i2max ; ++i22 ){
              for(int i11=i1min; i11<=i1max ; ++i11 ){
                float d = sqrt((i33-i3)*(i33-i3)+
                               (i11-i1)*(i11-i1)+
                                (i22-i2)*(i22-i2));
                if(d<=r){
                  c[i33][i22][i11] += f;      
                  n[i33][i22][i11]++;
                }
              } // i11
            } // i22
          } // i33

        } // i1
      } // i2
    } // i3

    for (int i3=0; i3<n3 ; ++i3){
      for (int i2=0; i2<n2 ; ++i2){
        for (int i1=0; i1<n1 ; ++i1){
          c[i3][i2][i1] = c[i3][i2][i1]/n[i3][i2][i1];
        } // i1
      } // i2
    } // i3
    System.out.printf("Sibson in %g sec%n",sw.time());
  }




  public void applyR(float[][][] scatter, float[][][] c ){
    int n1 = scatter[0][0].length;
    int n2 = scatter[0].length;
    int n3 = scatter.length;

    // get closest point
    ClosestPointTransform cpt = new ClosestPointTransform(_flag);
    int[][][] i1o = new int[n3][n2][n1];
    int[][][] i2o = new int[n3][n2][n1];
    int[][][] i3o = new int[n3][n2][n1];
    float[][][] dt = new float[n3][n2][n1];

    System.out.printf("===========================================%n");
    System.out.printf("calculating Closest point Transform (CPT) %n");
    System.out.printf("n1=%4d n2=%4d %n",n1,n2);
    Stopwatch sw = new Stopwatch();
    sw.start();

    cpt.apply(scatter,i1o,i2o,i3o,dt);

    System.out.printf("finished CPT in %g sec%n",sw.time());

    Vector<Integer> toDo = new Vector<Integer>();
    int i=0;
    for (int i3=0; i3<n3; ++i3){
      for (int i2=0; i2<n2; ++i2){
        for (int i1=0; i1<n1; ++i1){
          toDo.add(i3*n2*n1+i2*n1+i1);
        }
      }
    }

    // accumulated hits vector:
    int[][][] n = new int[n3][n2][n1];

    Iterator it = toDo.iterator();
    while(it.hasNext(){
      int i = it.next();
      int i3= i/(n2*n1); 
      int i2= (i-i3*n2*n1)/n1;
      int i1= i-i3*n2*n1 -i2*n1;


      System.out.printf("i3=%d%n",i3);

      int i1cp = i1o[i3][i2][i1];
      int i2cp = i2o[i3][i2][i1];
      int i3cp = i3o[i3][i2][i1];
      float r = dt[i3][i2][i1]; 
      float f = scatter[i3cp][i2cp][i1cp];
  
      int i1min = max(i1 -(int)r,0); 
      int i1max = min(i1 +(int)r,n1-1); 

      int i2min = max(i2 -(int)r,0); 
      int i2max = min(i2 +(int)r,n2-1); 

      int i3min = max(i3 -(int)r,0); 
      int i3max = min(i3 +(int)r,n3-1); 
  
      for(int i33=i3min; i33<=i3max ; ++i33 ){
        for(int i22=i2min; i22<=i2max ; ++i22 ){
          for(int i11=i1min; i11<=i1max ; ++i11 ){
            float d = sqrt((i33-i3)*(i33-i3)+
                           (i11-i1)*(i11-i1)+
                            (i22-i2)*(i22-i2));
            if(d<=r){
              c[i33][i22][i11] += f;      
              n[i33][i22][i11]++;
            }
          } // i11
        } // i22
      } // i33
    }

    for (int i3=0; i3<n3 ; ++i3){
      for (int i2=0; i2<n2 ; ++i2){
        for (int i1=0; i1<n1 ; ++i1){
          c[i3][i2][i1] = c[i3][i2][i1]/n[i3][i2][i1];
        } // i1
      } // i2
    } // i3
    System.out.printf("Sibson in %g sec%n",sw.time());
  }




  public void applyP(final float[][][] scatter, final float[][][] c ){
    final int n1 = scatter[0][0].length;
    final int n2 = scatter[0].length;
    final int n3 = scatter.length;

    // get closest point
    ClosestPointTransform cpt = new ClosestPointTransform(_flag);
    final int[][][] i1o = new int[n3][n2][n1];
    final int[][][] i2o = new int[n3][n2][n1];
    final int[][][] i3o = new int[n3][n2][n1];
    final float[][][] dt = new float[n3][n2][n1];

    System.out.printf("===========================================%n");
    System.out.printf("calculating Closest point Transform (CPT) %n");
    System.out.printf("n1=%4d n2=%4d n3=%4d %n",n1,n2,n3);
    Stopwatch sw = new Stopwatch();
    sw.start();

    cpt.apply(scatter,i1o,i2o,i3o,dt);

    System.out.printf("finished CPT in %g sec%n",sw.time());


    // accumulated hits vector:
    final int chunk =(int)n3/nthread;

    final float[][][][] e = Parallel.reduce(nthread, new Parallel.ReduceInt<float[][][][]>(){
      public float[][][][] compute(int ithread){
      float[][][][] e = new float[2][n3][n2][n1];

      int lb = ithread*chunk;
      int ub = min((ithread+1)*chunk,n3);
      if (ithread == nthread -1) ub = n3;

     // System.out.printf(" ithread =%d  for( int i2=%d;  i2<%d; ++i2)%n",ithread,lb,ub); 

        for (int i3=lb; i3<ub; ++i3){ 
          for (int i2=0; i2<n2 ; ++i2){
            for (int i1=0; i1<n1 ; ++i1){
              int i1cp = i1o[i3][i2][i1];
              int i2cp = i2o[i3][i2][i1];
              int i3cp = i3o[i3][i2][i1];

              float r = dt[i3][i2][i1]; 
              float f = scatter[i3cp][i2cp][i1cp];
      
              int i1min = max(i1 -(int)r,0); 
              int i1max = min(i1 +(int)r,n1-1); 
 
              int i2min = max(i2 -(int)r,0); 
              int i2max = min(i2 +(int)r,n2-1); 

              int i3min = max(i3 -(int)r,0); 
              int i3max = min(i3 +(int)r,n3-1); 
      
                  
              for(int i33=i3min; i33<=i3max ; ++i33 ){
                for(int i22=i2min; i22<=i2max ; ++i22 ){
                  for(int i11=i1min; i11<=i1max ; ++i11 ){
                    float d = sqrt((i11-i1)*(i11-i1)+(i22-i2)*(i22-i2)+(i33-i3)*(i33-i3));
                    if(d<=r){
                      e[0][i33][i22][i11] += f;      
                      e[1][i33][i22][i11]++;
                    }
                  }
                }
              }

            }
          }
        }
//         System.out.printf("i2 =%d%n",i2);
        return e;
      }

      public float[][][][] combine(float[][][][] ea, float[][][][] eb){
        return addT(ea,eb);
      }
    });
  
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2 ; ++i2){
        for (int i1=0; i1<n1 ; ++i1){
          c[i3][i2][i1] = e[0][i3][i2][i1]/e[1][i3][i2][i1];
        }
      }
    }});

    System.out.printf("Sibson in %g sec%n",sw.time());
  }

  
  private float[][][][] addT (final float[][][][] a,final float[][][][]b ){
    int n4 = a.length;
    int n3 = a[0].length;
    int n2 = a[0][0].length;
    int n1 = a[0][0][0].length;
    final float [][][][] c = new float[n4][n3][n2][n1];
  
    for (int i4=0; i4<n4 ; ++i4){
      c[i4] = add(a[i4],b[i4]);
    }

    return c;


  }

  public class Box{
    public int _min1,_max1,_min2,_max2,_min3,_max3;

    public Box(int min1, int max1, int min2,int max2,int min3,int max3){
      _min1=min1;
      _max1=max1;
      _min2=min2;
      _max2=max2;
      _min3=min3;
      _max3=max3;
    }

    public boolean intersect(Box other){
      return _max1 >= other._min1 &&
             _min1 <= other._max1 &&
             _max2 >= other._max2 &&
             _min2 <= other._max2 &&
             _max3 >= other._max3 &&
             _min3 <= other._max3; 
    }
  }


}
