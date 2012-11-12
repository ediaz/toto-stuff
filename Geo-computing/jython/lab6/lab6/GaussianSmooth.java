package lab6;

import edu.mines.jtk.util.Parallel;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Digital signal processing.
 * @author Esteban Diaz
 * @version 2012.11
 *


*/
public class GaussianSmooth{

public static int nthread = Runtime.getRuntime().availableProcessors();

  /*
    Common variables:
  */


  private float _sigma;
  private int _flength, _ngk1,_n1; 
  private int _nfft1, _nfft2, _k1  ; 
  private int _nfft  ;
  private int _ngk ;  
  private float [] _gfilter1; 
  private FftReal _fft; 


  /*
    Constructor for exponential smoothing
    sigma is an approximation to the standard
    deviation of the gaussian.
  */
  public  GaussianSmooth(float sigma){
    this(sigma,(int)(sigma*3.0));

    /*
       default size is related to  3sigma
       which is the distance for which a gaussian
       goes almost to zero. 
     */
  }

  public  GaussianSmooth(float sigma,int flength){
    _sigma = sigma;
    _flength = flength; 
    _nfft1 = FftReal.nfftFast(flength);

    if(flength < (int)(3.0*sigma)){ 
      System.out.printf("Warning: the filter length is too small! %n"); 
      System.out.printf("Recommended length = 3*sigma = %f %n",3.0f*sigma) ;
    }
  }

  public void apply(float[] x, float[] y){
    PrepareFft(x);
    getFilter1();
    filter1d_paux(x,y);
  } 

  public void apply(final float[][] x,final float[][] y){
    final int n2 = x.length; 
    final int n1 = x[0].length;

    PrepareFft(x[0]);
    getFilter1();

    // Filter along first axis
    filter1_2d(x,y);
    float [][] z =zerofloat(n2,n1) ;

    Transpose(y,z);
    PrepareFft(z[0]);
    getFilter1();
  
    // Filter along second axis
    filter1_2d(z,z);
    Transpose(z,y);
  }




  private void filter1_2d(final float [][] x, final float[][] y){
    final int n2 = x.length; 
    final int n1 = x[0].length;

    Parallel.loop(n2, new Parallel.LoopInt(){
      public void compute(int i2){
        filter1d_paux(x[i2],y[i2]);
      }
    });
  }

  //////////////////////////////////////////////////////////////////
  // private

  /** 
   * Smooth along second dimension in chunks along first dimension.
   * 
   */

  private  void getFilter1(){
    _gfilter1 = new float[(int)_nfft1/2 + 1];

    int n1 = _gfilter1.length;
    float scale = 1.f/_nfft1 ;
    float omega = -2.0f*FLT_PI*scale;
    float a1 =(float)(omega*omega*_sigma*_sigma*0.5);

    for(int i1=0; i1< n1 ; ++i1){
      _gfilter1[i1] =exp(-a1*i1*i1);
    }
  }


  private void filter1d_paux(float[] x, float[] y){
    float[] xaux = new float[_nfft1+2];

    copyZ(x,xaux);

    _fft.realToComplex(-1,xaux,xaux);

    for (int i1=0 ; i1<(int)_nfft1/2 +1 ; ++i1){
       int rk = 2*i1;
       xaux[rk  ] *= _gfilter1[i1];
       xaux[rk+1] *= _gfilter1[i1];
    }
    _fft.complexToReal(1, xaux, xaux);
    _fft.scale(_n1, xaux);

    copy(_n1, xaux,y);
  }

  private void copyZ(float[]x, float[] xaux){
    copy(x,xaux);
    int naux= xaux.length;

    for (int i1=_n1; i1<2*_n1;++i1)
      xaux[i1]=x[_n1-1];

    for (int i1=2*_n1; i1<_nfft1;++i1)
      xaux[i1]=x[0];

  }

  private void PrepareFft(float[]x){
    _n1 = x.length;
    _nfft1 = FftReal.nfftFast(3*_n1);
    _fft = new FftReal(_nfft1);
    
    if(_flength > _nfft1-1) 
      _flength = _nfft1;
  }

  private void Transpose(final float[][] x, final float[][] y){
  final int n2 = y.length; 
  final int n1 = y[0].length ;

  int k2=n2%4 ; 

  for (int i2=0;i2<k2 ; ++i2){
    for (int i1=0; i1<n1; ++i1){ 
      y[i2][i1]=x[i1][i2];
    }
  }
  

  Parallel.loop(k2,n2,4, new Parallel.LoopInt(){
    public void compute(int i2){
      for (int i1=0; i1<n1; ++i1){ 
        y[i2  ][i1]=x[i1][i2  ];
        y[i2+1][i1]=x[i1][i2+1];
        y[i2+2][i1]=x[i1][i2+2];
        y[i2+3][i1]=x[i1][i2+3];
      }
    }
    });
  }

}
