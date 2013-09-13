

package ediaz.rsf;

import ediaz.lib.*;
import java.util.*;
import java.io.IOException;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import rsf.RSF;
import rsf.Input;
import rsf.Output;
import edu.mines.jtk.util.Parallel;


public class Cicpvector
{
	static { 
	  System.loadLibrary("jrsf");
	}

	public static void main (String[] args){
		RSF par = new RSF(args);

    // I-O files
		Input  Fswfl = new Input("in");
		Input  Frwfl = new Input("rwfl");
		Input  Fspv  = new Input("spv");
		Input  Frpv  = new Input("rpv");
    Output Fout  = new Output("out");
    Output Fimg  = new Output("img");

    float PI = FLT_PI;

    /* time axis geometry */
    float o3 = Fswfl.getOrigin(3); 
    float d3 = Fswfl.getDelta(3);
    int   n3 = Fswfl.getN(3);

    /* x axis geometry */
    float o2 = Fswfl.getOrigin(2); 
    float d2 = Fswfl.getDelta(2);
    int   n2 = Fswfl.getN(2);

    /* z axis geometry */
    float o1 = Fswfl.getOrigin(1); 
    float d1 = Fswfl.getDelta(1);
    int   n1 = Fswfl.getN(1);


    Fout.setN(1,n1); Fout.setDelta(1,d1); Fout.setOrigin(1,o1);
    Fout.setN(2,n2); Fout.setDelta(2,d2); Fout.setOrigin(2,o2);
    Fout.setN(3,1) ; Fout.setDelta(3,d3); Fout.setOrigin(3,o3);

    Fimg.setN(1,n1); Fimg.setDelta(1,d1); Fimg.setOrigin(1,o1);
    Fimg.setN(2,n2); Fimg.setDelta(2,d2); Fimg.setOrigin(2,o2);
    Fimg.setN(3,1) ; Fimg.setDelta(3,d3); Fimg.setOrigin(3,o3);

    /* smoothing parameters from user */
    float a = par.getFloat("minangle",5.f);
    float b = par.getFloat("maxangle",10.f);
    float sigma = par.getFloat("sigma",-1.f);
        


    /* allocate working arrays */
    float[][] swfl  = zerofloat(n1,n2);
    float[][] rwfl  = zerofloat(n1,n2);
    float[][] img   = zerofloat(n1,n2);
    float[][] nfimg = zerofloat(n1,n2);
    float[][] imgsm = zerofloat(n1,n2);
    float[][][] spv = zerofloat(n1,n2,2);
    float[][][] rpv = zerofloat(n1,n2,2);
    
    float theta = 0.0f;
    System.err.printf("a=%g b=%g\n",a,b);
    a = cos(a*PI/180.f);
    b = cos(b*PI/180.f);
    float s = abs(a-b)/3.0f;
    int perc = n3/10;

    for (int i3=0; i3<n3; ++i3){
      Fswfl.read(swfl);  
      Frwfl.read(rwfl);  
      Fspv.read(spv);  
      Frpv.read(rpv);
      if(i3%perc == 0)
        System.err.printf("%3.1f %s done\n",i3/perc*10.0f,"%");

      if(a>b){
        goimage(spv,rpv,swfl,rwfl,img,nfimg,a,s);
      }else{
        goimage2(spv,rpv,swfl,rwfl,img,nfimg,b,s);
      }
    }

    if (sigma >0.0f){
      ExpSmooth xsm = new ExpSmooth(sigma);
      xsm.apply(img,img);
    }
    Fout.write(img);
    Fimg.write(nfimg);

    Fswfl.close();
    Frwfl.close();
    Fspv.close();
    Frpv.close();
    Fout.close();
    Fimg.close();
  }






  /****************** Private methods ***********************/
  private static float theta_w(float a,float s,float theta){
    float w2 = 0.0f;
    if(theta >a){
      w2 = 1.0f;
    }else{
      w2 = exp(-((theta-a)*(theta-a))/(2.0f*s*s));
    }
    if (w2 > 1.0f)
      System.err.printf("something funny here....\n");
    return w2;
  }

  private static float theta_w2(float a,float s,float theta){
    float w2 = 0.0f;
    if(theta <a){
      w2 = 1.0f;
    }else{
      w2 = exp(-((theta-a)*(theta-a))/(2.0f*s*s));
    }
    if (w2 > 1.0f)
      System.err.printf("something funny here....\n");
    return w2;
  }


  private static void goimage(
    final float [][][] spv, final float [][][] rpv,
    final float [][] swfl, final float[][] rwfl, 
    final float [][] img, final float [][] nfimg, 
    final float a, final float s){
    
    final int n2 = swfl.length;
    final int n1 = swfl[0].length;


    Parallel.loop(n2,new Parallel.LoopInt() {   
    public void compute(int i2) {
    //for( int i2=0; i2<n2;++i2){
      for (int i1=0; i1<n1; ++i1){
        float xspv = spv[0][i2][i1];
        float zspv = spv[1][i2][i1];
        float xrpv = rpv[0][i2][i1];
        float zrpv = rpv[1][i2][i1];

        float smag = sqrt(xspv*xspv+zspv*zspv);
        float rmag = sqrt(xrpv*xrpv+zrpv*zrpv);
        float pmag = smag*rmag;
        float num = (xspv*xrpv+zspv*zrpv);
        float costheta = 10.0f;
        if (pmag>0.0f) 
          costheta = num/pmag ;
        img[i2][i1] += theta_w(a,s,costheta)*swfl[i2][i1]*rwfl[i2][i1]; 
        nfimg[i2][i1] += swfl[i2][i1]*rwfl[i2][i1];             
      }
    } 
    });
   }                           

  private static void goimage2(
    final float [][][] spv, final float [][][] rpv,
    final float [][] swfl, final float[][] rwfl, 
    final float [][] img, final float [][] nfimg, 
    final float a, final float s){
    
    final int n2 = swfl.length;
    final int n1 = swfl[0].length;


    Parallel.loop(n2,new Parallel.LoopInt() {   
    public void compute(int i2) {
    //for( int i2=0; i2<n2;++i2){
      for (int i1=0; i1<n1; ++i1){
        float xspv = spv[0][i2][i1];
        float zspv = spv[1][i2][i1];
        float xrpv = rpv[0][i2][i1];
        float zrpv = rpv[1][i2][i1];

        float smag = sqrt(xspv*xspv+zspv*zspv);
        float rmag = sqrt(xrpv*xrpv+zrpv*zrpv);
        float pmag = smag*rmag;
        float num = (xspv*xrpv+zspv*zrpv);
        float costheta = 10.0f;
        if (pmag>0.0f) 
          costheta = num/pmag ;
        img[i2][i1] += theta_w2(a,s,costheta)*swfl[i2][i1]*rwfl[i2][i1]; 
        nfimg[i2][i1] += swfl[i2][i1]*rwfl[i2][i1];             
      }
    } 
    });
   }    



}

