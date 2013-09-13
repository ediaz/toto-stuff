package ediaz.lib.fdmod; 

import edu.mines.jtk.util.Parallel;


/**
 * Simple 2D finite difference modeling; 
 * @author Esteban Diaz
 * @version 2013.06.08
 *


*/
public class Fdmod2d{

public static int nthread = Runtime.getRuntime().availableProcessors();

  /*
    Common variables:
  */

  private static float _dx, _dz, _dt;
  private static float _ox, _oz, _ot;
  private static int   _nx, _nz, _nt;
  private static int   _ntout;

  private static float[][] _vel, _den;

  private static float[][] _src_coo; // [ns][2] (x,z)
  private static int _ns;
  private static float[][] _wav; // [ns][nt]


  private static float[][] _rvr_coo; // [nr][2] (x,z)
  private static int _nr;
  private static float[][] _dat; // [nr][nt/jsnap]  

  private static int _jsnap;
  private static boolean _wavefield;

  /* finite difference coeff from Muir's paper for staggered grid
     used in variable density
      
    taken from Francesco's sg code, need to check with him....
    */
  private static float c1s = 1225.f/1024.f;
  private static float c2s = -1225/(1024.f*15);
  private static float c3s = 1225/(1024.f*125);
  private static float c4s = -1225/(1024.f*1715);

  private static float _c1sx, _c2sx, _c3sx, _c4sx;
  private static float _c1sz, _c2sz, _c3sz, _c4sz;


  /* 8th order finite difference coeff for constant 
     velocity (centered) */
  private static float c0c = -205.f/72.f;
  private static float c1c = 8.f/315.f;
  private static float c2c = -1.f/5.f;
  private static float c3c = 8.f/315.f;
  private static float c4c = -1.f/560.f;
  
  private static float _c0cx, _c1cx, _c2cx, _c3cx, _c4cx;
  private static float _c0cz, _c1cz, _c2cz, _c3cz, _c4cz;


  private static float[][] _up1, _up0, _um1, _utmp;
  private static float[][] _odat;
  private static float[][][] _owfl;
  private static float[][] _v2dt2;




  /*
    Constructor for exponential smoothing
    sigma is an approximation to the standard
    deviation of a gaussian filter.
  */
  public Fdmod2d(float ox, float dx, int nx,
                 float oz, float dz, int nz,
                 float ot, float dt, int nt,
                 float[][] vel, // velocity array
                 float[][] _src_coo, // source locations (x,z) pairs[ns][2]
                 float[][] wav, // source function [ns][nt]
                 int jsnap, 
                 boolean wavefield){

    // x grid set up
    _ox = ox;
    _dx = dx;
    _nx = nx;

    // z grid set up
    _oz = oz;
    _dz = dz;
    _nz = nz;

    // t grid set up
    _ot = ot;
    _dt = dt;
    _nt = nt;


    // velocity
    _vel = vel;

    // wavelet set up
    _wav = wav;

    // every how many time samples to keep?
    _jsnap = jsnap;  

    _ns = wav.length;
    _nt = wav[0].length;

    // set data length
    _ntout = (int) _nt/jsnap;
    _wavefield = wavefield;

    // set fd coef with dx,dz constant density:
    _c0cx = c0c/_dx;
    _c1cx = c1c/_dx;
    _c2cx = c2c/_dx;
    _c3cx = c3c/_dx;
    _c4cx = c4c/_dx;

    _c0cz = c0c/_dz;
    _c1cz = c1c/_dz;
    _c2cz = c2c/_dz;
    _c3cz = c3c/_dz;
    _c4cz = c4c/_dz;

    
    _up1  = new float[_nx][_nz];
    _up0  = new float[_nx][_nz];
    _um1  = new float[_nx][_nz];
    _utmp = new float[_nx][_nz];
    
    _odat = new float[_nr][_nt];

    if (_wavefield)
      _owfl = new float[_ntout][_nx][_nz];

    _v2dt2 = new float[_nx][_nz];
    for (int ix=0; ix< _nx; ++ix){
      for (int iz=0; iz< _nz; ++iz){
        _v2dt2[ix][iz] = _vel[ix][iz]*_vel[ix][iz]*_dt*_dt;
      }
    }
  }





  public static void forward_all(){

    int itout =0;
    for (int it=0; it<_nt; ++it){
      tstep(it);

      if (it%_jsnap ==0 && _wavefield){
        _owfl[itout] = _up0; 
        itout +=1;
      }
    }
  }









  public static void tstep(int it){
      laplacian8th(_up0,_utmp); // 8th order laplacian operator
      inject_source(_utmp,it); // source injection
      scalev(_v2dt2,_utmp);     // scale by velocity squared 

      // do the time steping:
      for (int ix=0; ix< _nx; ++ix){
        for (int iz=0; iz< _nz; ++iz){
          _up1[ix][iz] = 2.0f*_up0[ix][iz]-_um1[ix][iz]+_utmp[ix][iz];
        }
      }    

      // circulate pointers:
      _um1 = _up0;
      _up0 = _up1;
  }


  private static void laplacian8th(float[][] u, float[][] out){
    float c0 = _c0cx +_c0cz;

    for (int ix=4; ix< _nx-4; ++ix){
      for (int iz=4; iz < _nz-4; ++iz){
        out[ix][iz] = c0*u[ix][iz] +
                   _c1cx*(u[ix-1][iz]+u[ix+1][iz])+
                   _c2cx*(u[ix-2][iz]+u[ix+2][iz])+
                   _c3cx*(u[ix-3][iz]+u[ix+3][iz])+
                   _c4cx*(u[ix-4][iz]+u[ix+4][iz])+
                   _c1cz*(u[ix-1][iz]+u[ix+1][iz])+
                   _c2cz*(u[ix-2][iz]+u[ix+2][iz])+
                   _c3cz*(u[ix-3][iz]+u[ix+3][iz])+
                   _c4cz*(u[ix-4][iz]+u[ix+4][iz]);
      } // iz
    } // ix 
  }

  private static void laplacian(float[][] u, float[][] idensity){

  }


  private static void inject_source(float[][]utmp,int it){

    for (int is=0; is<_ns;++is){
      int ix =  (int)((_src_coo[is][0]-_ox)/_dx +0.5f);
      int iz =  (int)((_src_coo[is][1]-_oz)/_dz +0.5f);
      utmp[ix][iz] += _wav[is][it];
    }  
  }


  private static void scalev(float[][] v2dt2, float[][] utmp){
    for (int ix=0; ix< _nx; ++ix){
      for (int iz=0; iz< _nz; ++iz){
        utmp[ix][iz] *= v2dt2[ix][iz];
      }
    }
  }


}
