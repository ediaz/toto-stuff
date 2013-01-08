package lab8;

import edu.mines.jtk.util.Parallel;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;




public class Grid{

  private float _flag;


  public Grid(float flag){
     _flag = flag;
  } 






  /**
   *
   * 1D Closest point transform O(NxM) M: # of known samples
   *
   * @param flag: number with the flag for null values
   *
   */

  public void apply(
  float[] sparse, int[] i1csp, float[] distance)
  {
    int n1 = i1csp.length;
    int nk =0; //known number of samples (don't know yet)

    System.out.printf("counting known samples loop %n");  


    for (int i1=0; i1<n1 ; ++i1){ 
      if(sparse[i1] != _flag){ 
        nk +=1;
      } 
    }                                   
    System.out.printf("known samples = %d %n",nk);  
    
    int[] i1known = new int[nk];
    nk=0;

      for (int i1=0; i1<n1 ; ++i1){ 
        // fill output with known samples first
        if(sparse[i1] != _flag){ 
          i1known[nk] = i1;
          i1csp[i1] = i1;
          distance[i1] = 0.0f;
          nk+=1;
        } 
      }                                   

    System.out.printf("pass regridding of known samples loop %n");  
 
    int i11=0;
    int ig = 1;
    int ng = n1-nk;

    for (int i1=0; i1<n1 ; ++i1){ 
      if(sparse[i1] == _flag){ 
        float d2 = (float)(n1*n1);
        // largest posible distance in the array.
        float dmin = d2;                                
        for (int i3=0; i3<nk ; ++i3){
          d2 = (i1known[i3]-i1)*(i1known[i3]-i1); 
          if(d2 <dmin){
            i11 = i1known[i3];
            dmin = d2;
          }            
        }
        i1csp[i1] = i11;
        distance[i1] = sqrt(dmin);
      }
    }                                   
  } // end of csp_naive method.

  /**
   *
   * 2D Closest point transform O(NxM) M: # of known samples
   *
   * @param flag: number with the flag for null values
   *
   */
  public void apply( float[][] sparse, int[][] i1csp,
  int[][] i2csp, float[][] distance)
  {
    int n2 = i2csp.length;
    int n1 = i2csp[0].length;
    int nk =0; //known number of samples (don't know yet)

    System.out.printf("counting known samples loop %n");  


    for (int i2=0; i2<n2 ; ++i2){ 
      for (int i1=0; i1<n1 ; ++i1){ 
        if(sparse[i2][i1] != _flag){ 
          nk +=1;
        } 
      }                                   
    } 
    System.out.printf("known samples = %d %n",nk);  
    
    int[] i1known = new int[nk];
    int[] i2known = new int[nk];
    nk=0;

    for (int i2=0; i2<n2 ; ++i2){ 
      for (int i1=0; i1<n1 ; ++i1){ 
        // fill output with known samples first
        if(sparse[i2][i1] != _flag){ 
          i1known[nk] = i1;
          i2known[nk] = i2;
          i2csp[i2][i1] = i2;
          i1csp[i2][i1] = i1;
          distance[i2][i1] = 0.0f;
          nk+=1;
        } 
      }                                   
    }

    System.out.printf("passed regridding of known samples loop %n");  
 
    int i22 = 0; int i11=0;
    int ig = 1;
    int ng = n1*n2-nk;

    for (int i2=0; i2<n2 ; ++i2){ 
      for (int i1=0; i1<n1 ; ++i1){ 
        if(sparse[i2][i1] == _flag){ 
          float dmin = (float)(n1*n1+n2*n2);
          // largest posible distance in the array.
          for (int ik=0; ik<nk ; ++ik){
             float d2 = (i1known[ik]-i1)*(i1known[ik]-i1)
                 + (i2known[ik]-i2)*(i2known[ik]-i2); 
            if(d2 <dmin){
              i22 = i2known[ik];
              i11 = i1known[ik];
              dmin = d2;
            }            
          }
          i2csp[i2][i1] = i22;
          i1csp[i2][i1] = i11;
          distance[i2][i1] = sqrt(dmin);
        }
      }                                   
    } 
  } // end of csp_naive 2d method.



  /**
   *
   * 3D Closest point transform O(NxM) M: # of known samples
   *
   * This algorithm should work nicely for M<< N
   *
   * @param flag: number with the flag for null values
   *
   */
  public void apply( 
    float[][][] sparse, int[][][] i1csp,
    int[][][] i2csp, int[][][] i3csp, float[][][] distance)
  {
    
    int n3 = i2csp.length;
    int n2 = i2csp[0].length;
    int n1 = i2csp[0][0].length;
    int nk =0; //known number of samples (don't know yet)

    System.out.printf("counting known samples loop %n");  

    for (int i3=0; i3<n3 ; ++i3){
      for (int i2=0; i2<n2 ; ++i2){ 
        for (int i1=0; i1<n1 ; ++i1){ 
          if(sparse[i3][i2][i1] != _flag){ 
            nk +=1;
          } 
        }                                   
      } 
    }
    System.out.printf("known samples = %d %n",nk);  
    
    int[] i1known = new int[nk];
    int[] i2known = new int[nk];
    int[] i3known = new int[nk];
    nk=0;

    for (int i3=0; i3<n3 ; ++i3){
      for (int i2=0; i2<n2 ; ++i2){ 
        for (int i1=0; i1<n1 ; ++i1){ 
          // fill output with known samples first
          if(sparse[i3][i2][i1] != _flag){ 
            i1known[nk] = i1;
            i2known[nk] = i2;
            i3known[nk] = i3;
            i3csp[i3][i2][i1] = i3;
            i2csp[i3][i2][i1] = i2;
            i1csp[i3][i2][i1] = i1;
            distance[i3][i2][i1] = 0.0f;
            nk+=1;
          } 
        }                                   
      }
    }
    System.out.printf("passed regridding of known samples loop %n");  
 
    int i33 = 0; int i22 = 0; int i11=0;
    int ig = 1;
    int ng = n1*n2-nk;

    for (int i3=0; i3<n3 ; ++i3){
      for (int i2=0; i2<n2 ; ++i2){ 
        for (int i1=0; i1<n1 ; ++i1){ 
          if(sparse[i3][i2][i1] == _flag){ 
            float d2 = (float)(n1*n1 + n2*n2 + n3*n3);
            // largest posible distance in the array.
            float dmin = d2;                                
            for (int ik=0; ik<nk ; ++ik){
              d2 = (i1known[ik]-i1)*(i1known[ik]-i1)
                 + (i2known[ik]-i2)*(i2known[ik]-i2) 
                 + (i3known[ik]-i3)*(i3known[ik]-i3); 

              if(d2 <dmin){
                i33 = i3known[ik];
                i22 = i2known[ik];
                i11 = i1known[ik];
                dmin = d2;
              }            
            }
            i3csp[i3][i2][i1] = i33;
            i2csp[i3][i2][i1] = i22;
            i1csp[i3][i2][i1] = i11;
            distance[i3][i2][i1] = sqrt(dmin);
          }
        }                                   
      } 
    }
  } // end of csp_naive 3d method.

}
