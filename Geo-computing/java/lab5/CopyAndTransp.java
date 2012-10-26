import edu.mines.jtk.util.*;
import javax.script.ScriptEngine;



/*
  Homework 5:

  1) To run this program (edu_mines_jtk.jar has to be in 
     ../):

     ./j CopyAndTransp

  2) System report:

      Model Name:	MacBook Pro
      Model Identifier:	MacBookPro8,2
      Processor Name:	Intel Core i7
      Processor Speed:	2.4 GHz
      Number of Processors:	1
      Total Number of Cores:	4
      L2 Cache (per Core):	256 KB
      L3 Cache:	6 MB
      Memory:	8 GB
      Boot ROM Version:	MBP81.0047.B27
      SMC Version (system):	1.69f3
      Serial Number (system):	C02GN02FDY5R
      Hardware UUID:	2DB179F7-48FB-5F04-8283-1EEB1F2FA1DF
      Sudden Motion Sensor:
      State:	Enabled

  3) Printout of this program in my machine:

      ================================================
           Benchmarking array Copy
        (n1,n2,n3)=(101,102,103); rate=9846 floats/s
        (n1,n2,n3)=(201,202,203); rate=3777 floats/s
        (n1,n2,n3)=(301,302,303); rate=3958 floats/s
        (n1,n2,n3)=(401,402,403); rate=4054 floats/s
        (n1,n2,n3)=(501,502,503); rate=3961 floats/s
      ================================================
           Benchmarking transpose13S1 (Serial) 
        (n1,n2,n3)=(101,102,103); rate= 682 floats/s
        (n1,n2,n3)=(201,202,203); rate= 277 floats/s
        (n1,n2,n3)=(301,302,303); rate=  53 floats/s
        (n1,n2,n3)=(401,402,403); rate=  39 floats/s
        (n1,n2,n3)=(501,502,503); rate=  21 floats/s
      ================================================
           Benchmarking transpose13P1 (Parallel) 
        (n1,n2,n3)=(101,102,103); rate=1791 floats/s
        (n1,n2,n3)=(201,202,203); rate= 356 floats/s
        (n1,n2,n3)=(301,302,303); rate= 211 floats/s
        (n1,n2,n3)=(401,402,403); rate= 140 floats/s
        (n1,n2,n3)=(501,502,503); rate= 118 floats/s
      ================================================
           Benchmarking transpose13S8 (Serial) 
        (n1,n2,n3)=(101,102,103); rate= 664 floats/s
        (n1,n2,n3)=(201,202,203); rate= 886 floats/s
        (n1,n2,n3)=(301,302,303); rate= 721 floats/s
        (n1,n2,n3)=(401,402,403); rate= 531 floats/s
        (n1,n2,n3)=(501,502,503); rate= 409 floats/s
      ================================================
           Benchmarking transpose13P8 (Parallel) 
        (n1,n2,n3)=(101,102,103); rate=2114 floats/s
        (n1,n2,n3)=(201,202,203); rate=2053 floats/s
        (n1,n2,n3)=(301,302,303); rate=2166 floats/s
        (n1,n2,n3)=(401,402,403); rate=1471 floats/s
        (n1,n2,n3)=(501,502,503); rate=1248 floats/s


    4) Copying performance did not change with different
       approaches, I guess everything is as optimized
       as possible when copying.
*/ 




/**
 * Array copy and transposing
 * @author Esteban Diaz
 *         Colorado School of Mines
 * @version 2012.10.24
 *
*/
public class CopyAndTransp {



  public static void main(String[] args) {

    int [] n = new int[] {101,201,301,401,501};

    int ns = n.length;
    
    /* Benchmarking Copy (reference fast rate)
       seems to be invariable with the size
       of the array
   */ 
    System.out.printf("================================================\n");
    System.out.printf("     Benchmarking array Copy\n");
    for (int i=0 ; i<ns ; ++i)
      bench_interfaced("copy1P",n[i],n[i]+1,n[i]+2);

    System.out.printf("================================================\n");
    System.out.printf("     Benchmarking transpose13S1 (Serial) \n");
    for (int i=0 ; i<ns ; ++i)
      bench_interfaced("transp13S1",n[i],n[i]+1,n[i]+2);

    System.out.printf("================================================\n");
    System.out.printf("     Benchmarking transpose13P1 (Parallel) \n");
    for (int i=0 ; i<ns ; ++i)
      bench_interfaced("transp13P1",n[i],n[i]+1,n[i]+2);

    System.out.printf("================================================\n");
    System.out.printf("     Benchmarking transpose13S8 (Serial) \n");
    for (int i=0 ; i<ns ; ++i)
      bench_interfaced("transp13S8",n[i],n[i]+1,n[i]+2);

    System.out.printf("================================================\n");
    System.out.printf("     Benchmarking transpose13P8 (Parallel) \n");
    for (int i=0 ; i<ns ; ++i)
      bench_interfaced("transp13P8",n[i],n[i]+1,n[i]+2);


  }


  /* 
    this following interface is quite ugly,
    I tried better but I could not find a way to pass 
    a function name as parameter.

    Still is better than creating many benchmark 
    functions...

    Starting the stopwatch outside the IF statements
    did not affect performance
 
  */
  public static void bench_interfaced(String fname, int n1,int n2, int n3){
    float [][][] x  = new float[n3][n2][n1];
    float [][][] y  = new float[n1][n2][n3];
    float [][][] z  = new float[n3][n2][n1];

    randArray3(x);

    int ncopy=0;
    double maxtime = 2.0;

    int rate=0; 

 

    Stopwatch sw = new Stopwatch();
    sw.start();
    if(fname=="transp13P1"){
      for (ncopy=0; sw.time()<maxtime; ncopy+=2){
        Transp.transp13P1(x,y);
        Transp.transp13P1(y,z);
      }
 
    }else if (fname=="transp13P8"){
       for (ncopy=0; sw.time()<maxtime; ncopy+=2){
         Transp.transp13P8(x,y);
         Transp.transp13P8(y,z);
        }
    }else if (fname == "transp13S8"){
       for (ncopy=0; sw.time()<maxtime; ncopy+=2){
         Transp.transp13S8(x,y);
         Transp.transp13P8(y,z);
       }

    }else if (fname == "transp13S1"){
       for (ncopy=0; sw.time()<maxtime; ncopy+=2){
         Transp.transp13S1(x,y);
         Transp.transp13S1(y,z);
        }

    }else if (fname == "copy1P"){
       for (ncopy=0; sw.time()<maxtime; ncopy+=2){
         Copy.copy1P(x,z);
         Copy.copy1P(z,x);
        }
    }else{
      System.out.printf("function not recognized%n");
      System.exit(2);
    } 

    sw.stop();
    rate = (int)(1e-6*n1*n2*n3*ncopy/sw.time()) ;
  
    if (fname == "copy1P"){
      AssertFloats(x,z); 
    }else{
      AssertFloats(z,x); 

    }

    printrate( n1, n2, n3,rate );

  }





  public static void randArray3(float [][][] x){
  
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;

    RandomFloat r = new RandomFloat(1);

    for (int i3=0 ; i3<n3; ++i3)
      for (int i2=0 ; i2<n2; ++i2)
        for (int i1=0 ; i1<n1; ++i1)
          x[i3][i2][i1] = r.uniform();
    
  }






  public static void AssertFloats(float [][][] x, float [][][] y){

    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;


    for (int i3=0 ; i3<n3; ++i3)
      for (int i2=0 ; i2<n2; ++i2)
        for (int i1=0 ; i1<n1; ++i1)
          assert x[i3][i2][i1] == y[i3][i2][i1]: "x != y ";
  }



  private static void printrate(int n1,int n2, int n3, int rate){

    System.out.printf("  (n1,n2,n3)=(%d,%d,%d); rate=%4d floats/s\n",
            n1,n2,n3,rate);
  }


}
