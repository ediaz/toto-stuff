package fdmod; 

import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Simple 2D muting functions; 
 * @author Esteban Diaz
 * @version 2013.06.08
 *


*/
public class Mute{


  /*
    Common variables:
  */

  private static float _dx, _dt;
  private static float _ox, _ot;
  private static int   _nx, _nt;

  public static float _zr, _otsrc,_ramp;

  private static float[][] _src_coo; // [ns][2] (x,z)

  private static float[][] _rvr_coo; // [nr][2] (x,z)


  /*
  */
  public Mute(float ox, float dx, int nx,
                 float ot, float dt, int nt,
                 float otsrc,
                 float zr, float ramp){

    // x grid set up
    _ox = ox;
    _dx = dx;
    _nx = nx;

    // t grid set up
    _ot = ot;
    _dt = dt;
    _nt = nt;
  
    // receiver depth, same for all
    _zr = zr;

    _otsrc = otsrc;
    _ramp  = ramp;
  }

  public static float[][] linear(float xs, float zs, float vel, float[][] igather){

    float[][] muted = new float[_nx][_nt];
    float zz = abs(_zr-zs);

    for (int ix=0; ix<_nx; ++ix){
      float x = ix*_dx + _ox;
      float xx = abs(x - xs);
      float d = sqrt(xx*xx+zz*zz);
      float tmute = d/vel + _otsrc;

      for (int it=0; it<_nt; ++it){
        float t = it*_dt +_ot;
        if (t>tmute){
          muted[ix][it] = igather[ix][it];
        }
        if ( t >tmute - _ramp && t<tmute){
          float dt = t -(tmute-_ramp);
          muted[ix][it]=igather[ix][it]*dt/_ramp;
        } 
        if (t<tmute - _ramp){
          muted[ix][it] = 0.0f;
        }        
      }
    }
    return muted;
  }

  public static void linear(float xs, float zs, float vel, float[][] igather,float[][] muted){

    float zz = abs(_zr-zs);

    for (int ix=0; ix<_nx; ++ix){
      float x = ix*_dx + _ox;
      float xx = abs(x - xs);
      float d = sqrt(xx*xx+zz*zz);
      float tmute = d/vel + _otsrc;

      for (int it=0; it<_nt; ++it){
        float t = it*_dt +_ot;
        if (t>tmute){
          muted[ix][it] = igather[ix][it];
        }
        if ( t >tmute - _ramp && t<tmute){
          float dt = t -(tmute-_ramp);
          muted[ix][it]=igather[ix][it]*dt/_ramp;
        } 
        if (t<tmute - _ramp){
          muted[ix][it] = 0.0f;
        }        
      }
    }
  }








  public static float[][] hyperbolic(float xs, float zs, float vel, float[][] igather){

    float[][] muted = new float[_nx][_nt];
    float zz = abs(_zr-zs);

    for (int ix=0; ix<_nx; ++ix){
      float x = ix*_dx + _ox;
      float xx = abs(x - xs);
      float d = sqrt(xx*xx+zz*zz);
      float tmute = sqrt( d*d/(vel*vel) +_otsrc*_otsrc);

      for (int it=0; it<_nt; ++it){
        float t = it*_dt +_ot;
        if (t>tmute){
          muted[ix][it] = igather[ix][it];
        }
        if ( t >tmute - _ramp && t<tmute){
          float dt = t -(tmute-_ramp);
          muted[ix][it]= dt/_ramp*igather[ix][it];
        } 
        if (t<tmute - _ramp){
          muted[ix][it] = 0.0f;
        }        
      }
    }
    return muted;
  }


  public static void hyperbolic(float xs, float zs, float vel, float[][] igather,float[][] muted){

    float zz = abs(_zr-zs);

    for (int ix=0; ix<_nx; ++ix){
      float x = ix*_dx + _ox;
      float xx = abs(x - xs);
      float d = sqrt(xx*xx+zz*zz);
      float tmute = sqrt( d*d/(vel*vel) +_otsrc*_otsrc);

      for (int it=0; it<_nt; ++it){
        float t = it*_dt +_ot;
        if (t>tmute){
          muted[ix][it] = igather[ix][it];
        }
        if ( t >tmute - _ramp && t<tmute){
          float dt = t -(tmute-_ramp);
          muted[ix][it]= dt/_ramp*igather[ix][it];
        } 
        if (t<tmute - _ramp){
          muted[ix][it] = 0.0f;
        }        
      }
    }
  }



}
