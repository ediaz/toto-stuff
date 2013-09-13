package radon; 

import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.dsp.FftReal;

/**
 * Linear and parabolic radon transform:
 * forward and adjoint operators 
 * @author Esteban Diaz
 * @version 2013.06.08
 *


*/
public class Radon{

  /*
    Common variables:
  */

  private static float _dx, _dp, _dt, _dw;
  private static float _ox, _op, _ot, _ow;
  private static int   _nx, _np, _nt, _nw;
  private static int   _nfft1;
  private static FftReal _fft; 
  private static float[][] _wpanelx, _wpanelp;


  public Radon(float ox, float dx, int nx,
               float op, float dp, int np,
               float ot, float dt, int nt){
    // x grid set up
    _ox = ox;
    _dx = dx;
    _nx = nx;

    // p grid set up
    _op = op;
    _dp = dp;
    _np = np;

    // t grid set up
    _ot = ot;
    _dt = dt;
    _nt = nt;

    PrepareFft();
  }

  public static float[][] Forward(float[][] inpanel){
    float[][] outpanel = new float[_nx][_nt]; 
    Forward(inpanel,outpanel);
    return outpanel;
  }

  public static void Forward(float[][] inpanel,float[][] outpanel){
    goFft(inpanel); // this one fills wpanelx    
    
    for (int ix=0; ix<_nx; ++ix){
      for (int i=0; i<_wpanelx[0].length; ++i){
        _wpanelx[ix][i] = 0.0f;
      }
    }
 
    for (int iw=0; iw<_nw; ++iw){
      int ireal = iw*2;
      int iimag = ireal+1;
      float w = iw*_dw;
      for (int ix=0; ix<_nx; ++ix){
        float x = ix*_dx+_ox;
        for(int ip=0; ip<_np; ++ip){
          float p = ip*_dp +_op;
          float phi = w*p*x;  
          _wpanelx[ix][ireal] +=  _wpanelp[ip][ireal]*cos(phi)+ _wpanelp[ip][iimag]*sin(phi);
          _wpanelx[ix][iimag] += -_wpanelp[ip][ireal]*sin(phi)+ _wpanelp[ip][iimag]*cos(phi);
        }
      }
    }
    goiFft(outpanel);
    
  }



  public static float[][] Adjoint(float[][] inpanel){
    float[][] outpanel = new float[_np][_nt]; 
    Adjoint(inpanel,outpanel);
    return outpanel;
  }


  /**
   * This function takes data in (t,x) space and maps it to 
   * (t,p) space
   * @param inpanel: float[nx][nt] (t,x) panel
   * @param outpanel: float[np][nt] (t,p) panel 
   */
  public static void Adjoint(float[][] inpanel, float[][] outpanel){
    goFft(inpanel); // this one fills wpanelx    
    
    for (int ip=0; ip<_np; ++ip){
      for (int i=0; i<_wpanelp[0].length; ++i){
        _wpanelp[ip][i] = 0.0f;
      }
    }
 
    for (int iw=0; iw<_nw; ++iw){
      int ireal = iw*2;
      int iimag = ireal+1;
      float w = iw*_dw;
      for (int ip=0; ip<_np; ++ip){
        float p = ip*_dp +_op;
        for(int ix=0; ix<_nx; ++ix){
          float x = ix*_dx+_ox;
          float phi = w*p*x;
          _wpanelp[ip][ireal] += _wpanelx[ix][ireal]*cos(phi)- _wpanelx[ix][iimag]*sin(phi);
          _wpanelp[ip][iimag] += _wpanelx[ix][ireal]*sin(phi)+ _wpanelx[ix][iimag]*cos(phi);
        }
      }
    }
    goiFft(outpanel);
  }











  private static void goFft(float[][] panel){
    if (panel.length == _nx){
      for (int ix=0; ix<_nx;++ix){
        copy(_nt,panel[ix],_wpanelx[ix]); 
        _fft.realToComplex(-1,_wpanelx[ix],_wpanelx[ix]);
      }   
    }else{
      for (int ip=0; ip<_np;++ip){
        copy(panel[ip],_wpanelp[ip]);    
        _fft.realToComplex(-1,_wpanelp[ip],_wpanelp[ip]);
      }
    }  
  }


  private static void goiFft(float[][] panel){
    if (panel.length == _nx){
      for (int ix=0; ix<_nx;++ix){
        _fft.complexToReal(+1,_wpanelx[ix],_wpanelx[ix]); // ifft
        copy(_nt,_wpanelx[ix],panel[ix]); // copy to output array
        _fft.scale(_nt,panel[ix]);
      }   
    }else{
      for (int ip=0; ip<_np;++ip){
        _fft.complexToReal(+1,_wpanelp[ip],_wpanelp[ip]);
        copy(_nt,_wpanelp[ip],panel[ip]);
        _fft.scale(_nt,panel[ip]);
      }
    }  
  }





  private void PrepareFft(){
    _nfft1 = 2*FftReal.nfftFast(_nt);
    _fft = new FftReal(_nfft1);

    _nw = _nfft1/2 +1;
    _dw = 2.0f*FLT_PI/(_nfft1*_dt);
    _ow = 0.0f;    

    _wpanelx = new float[_nx][_nfft1+2];
    _wpanelp = new float[_np][_nfft1+2];

  }

}
