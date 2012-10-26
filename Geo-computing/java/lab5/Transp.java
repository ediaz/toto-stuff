import edu.mines.jtk.util.*;
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
public class Transp {



  public static void transp13S1(float[][][] x,float[][][] y ) {

    int n2 = y[0].length;

    for (int i2=0 ; i2 < n2 ; ++i2)
      transp13S1_chunki2(x,y,i2); 
  }


  public static void transp13S8(float[][][] x,float[][][] y ) {

    int n2 = y[0].length;

    for (int i2=0 ; i2 < n2 ; ++i2)
      transp13S8_chunki2(x,y,i2); 
  }




  public static void transp13P1(final float[][][] x, final float[][][] y ) {

    int n2 = y[0].length;

    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        transp13S1_chunki2(x,y,i2); 
      }
    });

  }





  public static void transp13P8(final float[][][] x, final float[][][] y ) {

    int n2 = y[0].length;

    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        transp13S8_chunki2(x,y,i2); 
      }
    });

  }



  

  private static void transp13S1_chunki2( 
    final float[][][] x,final float[][][] y,final int i2 ) {

    int n1 = y[0][0].length;
    int n3 = y.length;
  

    for (int i1=0 ; i1 < n1 ; ++i1) 
      for (int i3=0 ; i3 < n3 ; ++i3) 
          y[i3][i2][i1] = x[i1][i2][i3];

  }









  private static void transp13S8_chunki2( 
    final float[][][] x,final float[][][] y,final int i2 ) {

    int n1 = y[0][0].length;
    int n3 = y.length;
    
    int k3 = n3%8 ;

    
    for (int i3=0 ; i3 < k3 ; ++i3) 
        for (int i1=0 ; i1 < n1 ; ++i1) 
          y[i3][i2][i1] = x[i1][i2][i3];


    for (int i3=k3 ; i3 < n3 ; i3+=8) {
        for (int i1=0 ; i1 < n1 ; ++i1){ 
          y[i3  ][i2][i1] = x[i1][i2][i3  ];
          y[i3+1][i2][i1] = x[i1][i2][i3+1];
          y[i3+2][i2][i1] = x[i1][i2][i3+2];
          y[i3+3][i2][i1] = x[i1][i2][i3+3];
          y[i3+4][i2][i1] = x[i1][i2][i3+4];
          y[i3+5][i2][i1] = x[i1][i2][i3+5];
          y[i3+6][i2][i1] = x[i1][i2][i3+6];
          y[i3+7][i2][i1] = x[i1][i2][i3+7];
        }
     }
  }

}
