package lab8;

import edu.mines.jtk.util.Parallel;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;


/*
 *
 *  Closest Point Transform (CPT) for Euclidean Distance.
 *  
 *  This class implements CPT in 1, 2 and 3 dimensions
 *  The constructor just needs the flag value (from 
 *  which we can identify null samples).
 *
 *  The output consist on:
 *   1) Coordinates of the closest non-null point
 *   2) Distance to the closest non-null point
 *
 * 
 *
 *  A. MEIJSTER*‚  J.B.T.M.  ROERDINK and W.H. HESSELINK. 
 *  A GENERAL ALGORITHM FOR COMPUTING DISTANCE TRANSFORMS 
 *  IN LINEAR TIME.
 *
 *
 * @author: Esteban D'\{i}az
 *
 */

public class ClosestPointTransform{

  private float _flag, _inf;

  public ClosestPointTransform(float flag){
    _flag = flag;
  }

  /**
   * 1D closest point transform:  
   * 
   * @param sparse: input data with null values= flag 
   * @param i1cpt: closest point coordinate 
   * @param d: distance to closest point coordinate 
   * 
   */

  public void apply(float[] sparse, int[] i1cpt, float[]d){
    phaseI(sparse, d, i1cpt);
  }

  /**
   * 1D closest point transform:  
   * 
   * @param sparse: input data with null values= flag 
   * @param i1cpt: closest point 1 coordinate 
   * @param i2cpt: closest point 2 coordinate 
   * @param d: distance to closest point coordinate 
   * 
   */
  public void apply (
    final float[][] sparse, 
    final int[][] i1cpt, 
    final int[][] i2cpt, 
    final float[][]d)
  {
    
    final int n2 = sparse.length;    
    final int n1 = sparse[0].length;    
    
    _inf = (float)(n1*n1+n2*n2); //greatest posible distance


    // first do the distance transform
    // along i1 direction
  
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        phaseI(sparse[i2],d[i2],i1cpt[i2]);
        }
      });
    
    // after do the distance transform
    // along i2 direction  
    Parallel.loop(n1,new Parallel.LoopInt() {
      public void compute(int i1) {
        float[] g = new float[n2]; 
        float[] daux = new float[n2]; 
        int[] i2aux = new int[n2];
        int[] i1aux = new int[n2];
  
        for (int i2=0 ; i2<n2; ++i2){
          g[i2] = d[i2][i1];
        }
        phaseII(g,daux,i2aux);
  
        for (int i2=0 ; i2<n2 ; ++i2){
          d[i2][i1] = (daux[i2]);
          int i22 = i2aux[i2];
          i2cpt[i2][i1] = i22;
          i1aux[i2] = i1cpt[i22][i1];
        }
        for (int i2=0; i2<n2; ++i2){
          i1cpt[i2][i1]  = i1aux[i2];
        }
      }
      });
  }




  /**
   * 1D closest point transform:  
   * 
   * @param sparse: input data with null values= flag 
   * @param i1cpt: closest point 1 coordinate 
   * @param i2cpt: closest point 2 coordinate 
   * @param i3cpt: closest point 3 coordinate 
   * @param d: distance to closest point coordinate 
   */
  public void apply( 
    final float[][][] sparse, 
    final int[][][] i1cpt, 
    final int[][][] i2cpt,
    final int[][][] i3cpt, 
    final float[][][]d)
  {
    
    final int n3 = sparse.length;    
    final int n2 = sparse[0].length;    
    final int n1 = sparse[0][0].length;    
    
    _inf = (float)(n1*n1+n2*n2+n3*n3); //greatest posible distance

    // first do the distance transform
    // along i1 direction
  
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        float[] f = new float[n1]; 
        for (int i2=0 ; i2<n2; ++i2){
          for (int i1=0; i1<n1; ++i1){
            if(sparse[i3][i2][i1] == _flag){
              f[i1] = _inf;
            }else{ 
              f[i1] = 0.0f; 
            }
          }
          dt(f,d[i3][i2],i1cpt[i3][i2]);
        }
      }
      });


    // along i2 direction
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        float[] f = new float[n2]; 
        float[] d2 = new float[n2];
        int[] i2aux = new int[n2];
        int[] i1aux = new int[n2];

        for (int i1=0 ; i1<n1; ++i1){
          for (int i2=0 ; i2<n2; ++i2){
            f[i2] = d[i3][i2][i1];
          }
          dt(f,d2,i2aux);
          // update squared distance 
          for (int i2=0 ; i2<n2; ++i2){
            d[i3][i2][i1] = d2[i2];
            int i22 = i2aux[i2];
            i1aux[i2] = i1cpt[i3][i22][i1];
          }
          for (int i2=0; i2<n2; ++i2){
            int i11 = i1aux[i2];
            int i22 = i2aux[i2];
            i2cpt[i3][i2][i1] = i22;
            i1cpt[i3][i2][i1] = i11;
          }
        }
        
      }
      });



    // along i3 direction
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        float[] f = new float[n3]; 
        float[] d3 = new float[n3];
        int[] i3aux = new int[n3];
        int[] i2aux = new int[n3];
        int[] i1aux = new int[n3];

        for (int i1=0 ; i1<n1; ++i1){
          for (int i3=0 ; i3<n3; ++i3){
            f[i3] = d[i3][i2][i1];
          }
          dt(f,d3,i3aux);
          // update squared distance 
          for (int i3=0 ; i3<n3; ++i3){
            d[i3][i2][i1] = sqrt(d3[i3]);
            int i33 = i3aux[i3];
            i1aux[i3] = i1cpt[i33][i2][i1];
            i2aux[i3] = i2cpt[i33][i2][i1];
          }
          for (int i3=0 ; i3<n3; ++i3){
            int i33 = i3aux[i3];
            int i22 = i2aux[i3];
            int i11 = i1aux[i3];
            i3cpt[i3][i2][i1] = i33;
            i2cpt[i3][i2][i1] = i22;
            i1cpt[i3][i2][i1] = i11;
            
          }
        }
      }
      });
  }





  /**
   * 3D closest point transform:  
   * 
   * @param sparse: input data with null values= flag 
   * @param i1cpt: closest point 1 coordinate 
   * @param i2cpt: closest point 2 coordinate 
   * @param i3cpt: closest point 3 coordinate 
   * @param d: distance to closest point coordinate 
   */
  public void applyM( 
    final float[][][] sparse, 
    final int[][][] i1cpt, 
    final int[][][] i2cpt,
    final int[][][] i3cpt, 
    final float[][][]d)
  {
    
    final int n3 = sparse.length;    
    final int n2 = sparse[0].length;    
    final int n1 = sparse[0][0].length;    
    
    _inf = (float)(n1*n1+n2*n2+n3*n3); //greatest posible distance

    // first do the distance transform
    // along i1 direction
  
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        float[] f = new float[n1]; 
        for (int i2=0 ; i2<n2; ++i2){
          phaseI(sparse[i3][i2],d[i3][i2],i1cpt[i3][i2]);
        }
      }
      });


    // along i2 direction
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        float[] f = new float[n2]; 
        float[] d2 = new float[n2];
        int[] i2aux = new int[n2];
        int[] i1aux = new int[n2];

        for (int i1=0 ; i1<n1; ++i1){
          for (int i2=0 ; i2<n2; ++i2){
            f[i2] = d[i3][i2][i1];
          }
          phaseII(f,d2,i2aux);
          // update squared distance 
          for (int i2=0 ; i2<n2; ++i2){
            d[i3][i2][i1] = d2[i2];
            int i22 = i2aux[i2];
            i1aux[i2] = i1cpt[i3][i22][i1];
          }
          for (int i2=0; i2<n2; ++i2){
            int i11 = i1aux[i2];
            int i22 = i2aux[i2];
            i2cpt[i3][i2][i1] = i22;
            i1cpt[i3][i2][i1] = i11;
          }
        }
        
      }
      });



    // along i3 direction
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        float[] f = new float[n3]; 
        float[] d3 = new float[n3];
        int[] i3aux = new int[n3];
        int[] i2aux = new int[n3];
        int[] i1aux = new int[n3];

        for (int i1=0 ; i1<n1; ++i1){
          for (int i3=0 ; i3<n3; ++i3){
            f[i3] = d[i3][i2][i1];
          }
          phaseII(f,d3,i3aux);
          // update squared distance 
          for (int i3=0 ; i3<n3; ++i3){
            d[i3][i2][i1] = (d3[i3]);
            int i33 = i3aux[i3];
            i1aux[i3] = i1cpt[i33][i2][i1];
            i2aux[i3] = i2cpt[i33][i2][i1];
          }
          for (int i3=0 ; i3<n3; ++i3){
            int i33 = i3aux[i3];
            int i22 = i2aux[i3];
            int i11 = i1aux[i3];
            i3cpt[i3][i2][i1] = i33;
            i2cpt[i3][i2][i1] = i22;
            i1cpt[i3][i2][i1] = i11;
            
          }
        }
      }
      });
  }



  /**
   *
   * Phase I: 
   * simply does two 1d scans (forward a reverse)
   * to compute the distance along the fastest axis.
   *
   **/


  private void phaseI(float[] f, float[] d, int[] cp){
    int n = f.length;

    if(f[0]== _flag){
      d[0] = _inf;
    }else{
      d[0] = 0.0f; cp[0] = 0;
    }

    for(int i=1; i<n ;++i){
      if(f[i]== _flag){
        d[i] = d[i-1] +1.0f;
        cp[i] = cp[i-1];
      }else{
        d[i] = 0.0f;
        cp[i] =i;
      }
    }

    for(int i=n-2; i>=0 ; --i){
      if(d[i+1]<d[i]){
        d[i] = d[i+1]+1.f;
        cp[i]=cp[i+1];
      }
    }
  }



  /**
   * Phase II:
   * 
   * computes DT by keeping the lower envelope of an
   * stack of hyperbolas
   * 
   * implements phase II algorithm (using euclidean distance metric)  
   *
   *
   *
   */
  public static void phaseII(float[] g, float[] d, int[] cp){

    int n = g.length;
    int q = 0; 
    int[] s = new int[n+1]; 
    float[] t = new float[n+1];
  
    for (int u=1; u<n; ++u){
      while(q>=0 && f_edt(t[q],s[q],g)>f_edt(t[q],u,g)){
        q = q-1;
      }
      if(q<0){
        q = 0; s[0] = u;
      }else{
        float w = 1.0f +sep_edt(s[q],u,g);
        if (w<n) {
          q = q+1 ; s[q] =u ;
          t[q]= (int)w;
        }
      }
    }
    
    for (int u=n-1; u>=0; --u){
      d[u] = sqrt(f_edt(u,s[q],g));
      cp[u] = s[q]; 
      if(u <= t[q]){
        q -=1 ;
      }
    } 
  }




  public static void phaseII_new(float[] g, float[] d, int[] cp){

    int n = g.length;
    int q = 0; 
    int[] s = new int[n+1]; 
    float[] t = new float[n+1];
  
    for (int u=1; u<n; ++u){
      while(q>=0 && f_edt(t[q],s[q],g)>f_edt(t[q],u,g)){
        q = q-1;
      }
      if(q<0){
        q = 0; s[0] = u;
      }else{
        float w = 1.0f +sep_edt(s[q],u,g);
        if (w<n) {
          q = q+1 ; s[q] =u ;
          t[q]= (int)w;
        }
      }
    }
    
    for (int u=n-1; u>=0; --u){
      d[u] = (f_edt(u,s[q],g));
      cp[u] = s[q]; 
//      if(u <= t[q]){
      if(u == t[q]){
        q -=1 ;
      }
    } 
  }





  public static float f_edt_new( float x, int i, float[] g) {
    return (float)((x-i)*(x-i)+g[i]);
  }

  public static float f_edt( float x, int i, float[] g) {
    return (float)((x-i)*(x-i)+g[i]*g[i]);
  }

  public static float sep_edt_new(int i, int u, float[] g){
//    float s = (float)(((u+i)*(u-i)+(g[u]+g[i])*(g[u]-g[i]))/(2*u-2*i));
    //float s = ((float)((u*u +g[u]*g[u]) -(i*i+g[i]*g[i]))/(2*u-2*i));

    float s = (int)(u*u -i*i + g[u] - g[i])/(2*(u - i));

    return s; 
  }

 
  public static float sep_edt(int i, int u, float[] g){
    float s = (float)(((u+i)*(u-i)+(g[u]+g[i])*(g[u]-g[i]))/(2*u-2*i));


    return s; 
  }

  // TODO: create functions F and SEP for other
  //      metrics (Manhattan and Chessboard), as explained
  //      in the paper.


  private void dt(float[] f, float[] d, int[] cp){
    int n = f.length; 
    int[] v = new int[n];
    float[] z = new float[n+1];
    z[0] = -_inf ; z[1] = _inf;
    int k=0;
    for (int q=1; q<n ; ++q){
      float s  = ((f[q]+q*q)-(f[v[k]]+v[k]*v[k]))/(2*q-2*v[k]);
      while (s <= z[k]) {
        k=k-1;
        s  = ((f[q]+q*q)-(f[v[k]]+v[k]*v[k]))/(2*q-2*v[k]);
      }
      ++k;
      v[k] = q;
      z[k] = s;
      z[k+1] = _inf;
    } 
    k=0;
    for (int q=0 ;  q<= n-1; ++q){
      while (z[k+1] <q )
        k+=1;
      d[q] = (q-v[k])*(q-v[k]) +f[v[k]];
      cp[q] = v[k];
    }   
  }


}
