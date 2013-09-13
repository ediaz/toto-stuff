package fdmod; 

import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Simple 2D finite difference modeling; 
 * @author Esteban Diaz
 * @version 2013.06.08
 *


*/
public class Fdmod2d{

  /*
    Common variables:
  */

  private static float _dx, _dz, _dt, _dt2;
  private static float _ox, _oz, _ot;
  private static int   _nx, _nz, _nt;
  private static int   _ntout;

  private static float[][] _src_coo; // [ns][2] (x,z)
  private static int _ns;
  private static float[][] _wav; // [ns][nt]

  private static float[][] _rvr_coo; // [nr][2] (x,z)
  private static int _nr;
  private static float[][] _dat; // [nr][nt/jsnap]  

  private static int _jsnap = 10;
  private static boolean _wavefield=false;

  private static int _nb;

  /* 8th order finite difference coeff for constant 
     velocity (centered) */

  private static float _c0c = -205f/72.f;
  private static float _c1c = 8f/5f;
  private static float _c2c = -1.f/5.f;
  private static float _c3c = 8.f/315.f;
  private static float _c4c = -1.f/560f;


 
  private static float _c0cx, _c1cx, _c2cx, _c3cx, _c4cx;
  private static float _c0cz, _c1cz, _c2cz, _c3cz, _c4cz;


  public static float[][] _up1, _up0, _um1, _utmp,_ujsnap;
  private static float[][] _odat;
  private static float[][][] _owfl;
  public static float[][] _v2;

  private static boolean _free;
  private static boolean _record = false;


  private static int[] _abs ;


  /*
  */
  public Fdmod2d(float ox, float dx, int nx,
                 float oz, float dz, int nz,
                 float ot, float dt, int nt,
                 float[][] vel, // velocity array
                 int nb,
                 boolean free){

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

    setCoeficients();
    setBoundary(nb,free);
    setWavefields();
    setv2dt2(vel);
    _dt2 = _dt*_dt;
  }

  /**
   * The object can be recalled with new source locations if 
   * the modeling grid does not change.
  */

  public static void setSources(float[][] src_coo,float[][] wav){
    // wavelet set up
    _wav = wav;
    _src_coo = src_coo;
    _ns = src_coo.length;
  }

  public static void setReceivers(float[][] rvr_coo){
    _rvr_coo = rvr_coo;
    _nr = _rvr_coo.length;
    _odat = new float[_nr][_nt];
    _record = true;
  }

  private static void setBoundary(int nb, boolean free){
    _nb = nb;
    _free = free;
  }
  
  public static void goWavefield(int jsnap){
    _wavefield = true;
    _ntout = (int)(_nt/jsnap)+1;
    _jsnap = jsnap;
    _owfl = new float[_ntout][_nx][_nz];
  }

  private static void setCoeficients(){
    // set fd coef with dx,dz constant density:
    float sx = (1.f/_dx)*(1.f/_dx);
    float sz = (1.f/_dz)*(1.f/_dz);

    _c0cx = _c0c*sx;
    _c1cx = _c1c*sx;
    _c2cx = _c2c*sx;
    _c3cx = _c3c*sx;
    _c4cx = _c4c*sx;

    _c0cz = _c0c*sz;
    _c1cz = _c1c*sz;
    _c2cz = _c2c*sz;
    _c3cz = _c3c*sz;
    _c4cz = _c4c*sz;
  }

  private static void setv2dt2(float[][] vel){
    _v2 = new float[_nx+2*_nb][_nz+2*_nb];
    float[][] tmp = new float[_nx][_nz];
    for (int ix=0; ix< _nx; ++ix){
      for (int iz=0; iz< _nz; ++iz){
        float v = vel[ix][iz];
        tmp[ix][iz] = v*v;
      }
    }
    _v2 = extendModel(tmp);
    tmp = null;
  }
  
  private static void setWavefields(){
    _up1  = new float[_nx+2*_nb][_nz+2*_nb];
    _up0  = new float[_nx+2*_nb][_nz+2*_nb];
    _um1  = new float[_nx+2*_nb][_nz+2*_nb];
    _utmp = new float[_nx+2*_nb][_nz+2*_nb];
    _ujsnap = new float[_nx][_nz];
  }

  public static void setReverse(){
    int n1  = _wav[0].length;
    for (int is=0; is<_ns ;++is){
      float[] tmp = new float[n1];
      copy(_wav[is],tmp);
      for (int it=0; it<n1; ++it){
        _wav[is][it] = tmp[n1-1-it];
      }
    }
  }

  /*==========================================================*/
  public static void forwardAll(){
    setT0();
    int itout = 0;

    for (int it=0; it<_nt; ++it){
      tstep(it);

      if (it%_jsnap ==0 && _wavefield){
        for (int ix=0; ix< _nx; ++ix){
          for (int iz=0; iz< _nz; ++iz){
            _owfl[itout][ix][iz] = _up0[ix+_nb][iz+_nb];
        }}    
        itout +=1;
      }
  
      if(it%(_nt/10) == 0){

        System.out.printf("%4.1f %s done %n",it*100.0f/_nt,"%");
      }
    }
  }

  private static int _it=0,_itout=0;

  public static float[][] forwardJsnap(){
    if(_it ==0) setT0();

    int itout = 0;
    for (int it=_it; it<_nt; ++it){
      tstep(it);
      if (it%_jsnap ==0 && _wavefield){
        _it = it+1;
        break;
      }
    }

    for (int ix=0; ix< _nx; ++ix){
      for (int iz=0; iz< _nz; ++iz){
        _ujsnap[ix][iz] = _up0[ix+_nb][iz+_nb];
    }}  
    return _ujsnap;
  }


  public static float[][][] getWavefield(){
    return _owfl;
  }

  public static float[][] getData(){
    return _odat;
  }

  private static void setT0(){
      Parallel.loop(_nx+2*_nb,new Parallel.LoopInt() {
      public void compute(int ix) {
        for (int iz=0; iz< _nz+2*_nb; ++iz){
          _um1[ix][iz] = 0.0f;
          _up0[ix][iz] = 0.0f;
          _up1[ix][iz] = 0.0f;
        }
      }});
      for (int ir=0; ir<_nr; ++ir){
        for (int it=0; it<_nt; ++it){
        _odat[ir][it]=0.0f;
        }
      }
  }

  /*==========================================================*/


  public static void tstep(int it){
      laplacian8th(); // 8th order laplacian operator
      inject_source(it); // source injection

      // do the time steping:
      Parallel.loop(_nx+2*_nb,new Parallel.LoopInt() {
      public void compute(int ix) {
        for (int iz=0; iz< _nz+2*_nb; ++iz){
          _up1[ix][iz] = 2.0f*_up0[ix][iz]
                        -_um1[ix][iz]
                        +_v2[ix][iz]*_dt2*_utmp[ix][iz];
        }
      }});
      absorb();
      extract_data(it);
      // circulate pointers:
      Parallel.loop(_nx+2*_nb,new Parallel.LoopInt() {
      public void compute(int ix) {
        for (int iz=0; iz< _nz+2*_nb; ++iz){
          _um1[ix][iz] = _up0[ix][iz];
          _up0[ix][iz] = _up1[ix][iz];
      }}});
      if (_free){
      Parallel.loop(_nx+2*_nb,new Parallel.LoopInt() {
      public void compute(int ix) {
        for (int iz=0; iz< _nb; ++iz){
          _up0[ix][iz] = 0.0f;
        }
      }});
      }
  }

  public static void laplacian8th(){
    Parallel.loop(_nx+2*_nb-8,new Parallel.LoopInt() {
      public void compute(int i2) {
      int ix = i2+4;
      for (int iz=4; iz < _nz+2*_nb-4; ++iz){
        _utmp[ix][iz] = (_c0cx +_c0cz)*_up0[ix][iz]+
                         _c1cx*(_up0[ix-1][iz]+_up0[ix+1][iz])+
                         _c2cx*(_up0[ix-2][iz]+_up0[ix+2][iz])+
                         _c3cx*(_up0[ix-3][iz]+_up0[ix+3][iz])+
                         _c4cx*(_up0[ix-4][iz]+_up0[ix+4][iz])+
                         _c1cz*(_up0[ix][iz-1]+_up0[ix][iz+1])+
                         _c2cz*(_up0[ix][iz-2]+_up0[ix][iz+2])+
                         _c3cz*(_up0[ix][iz-3]+_up0[ix][iz+3])+
                         _c4cz*(_up0[ix][iz-4]+_up0[ix][iz+4]);
    }}}); 
  }

 //TODO: implement variable density acoustic wave equation



  public static void inject_source(final int it){
    Parallel.loop(_ns,new Parallel.LoopInt() {
      public void compute(int is) {    
      int ix =  (int)((_src_coo[is][0]-_ox)/_dx +0.5f);
      int iz =  (int)((_src_coo[is][1]-_oz)/_dz +0.5f);
   
      _utmp[ix+_nb][iz+_nb] += _wav[is][it];
    }});  
  }


  public static void extract_data(final int it){
    Parallel.loop(_nr,new Parallel.LoopInt() {
      public void compute(int ir) {    
      int ix =  (int)((_rvr_coo[ir][0]-_ox)/_dx +0.5f);
      int iz =  (int)((_rvr_coo[ir][1]-_oz)/_dz +0.5f);
   
      _odat[ir][it] += _up0[ix+_nb][iz+_nb];
    }});  

  }





  //============================================================


  public static boolean cflCheck(float fmax){
    boolean cfl;
    float maxvel = sqrt(max(_v2));
    cfl = (maxvel*_dt/min(_dx,_dz) < 1.0f);
    System.out.printf("cfl= %g%n",maxvel*_dt/(min(_dx,_dz)));
    return cfl;
  }


  public static boolean dispersionCheck(float fmax, float pmw){
    boolean dispersion;
    float minvel = sqrt(min(_v2));

    dispersion = (minvel/(fmax*pmw) >= max(_dx,_dz));
    if (dispersion==false)
      System.out.printf("I suggest dx = %g%n",minvel/(fmax*pmw));
    return dispersion;
  }


  /**
   * Methods based on Simon Luo's fd code,
   * taken with his permission from 
   * https://github.com/sluo/lss/blob/master/src/lss/dev/BornWavefield.java
   * on June 14th, 2013.
   */


  private static void absorb(){
    Parallel.loop(_nx+2*_nb,new Parallel.LoopInt() {
      public void compute(int ix) {
       absorbSliceX(ix ,_nb, _up1);
       absorbSliceX(ix ,_nb, _up0);

    }});

    absorbClayton(_um1, _up0, _up1);
  }


  private static void absorbSliceX(int ix ,int b, float[][] u) {
    /**
      Sponge type absorbing boundary condition, it dampens
      the wavefield as it enters into the absorbing layer.
      The dampending factor is:
        exp(-(x*x +z*z)*w*w)
        x: horizontal distance to the border of the bc layer
        z: vertical distance to the border of the bc  layer
        w: dampening factor, pretty ad-hoc...
     */

    int nz = u[0].length;
    int nx = u.length;
    int nzmb = nz-b;
    float w = 0.005f;
    float ws = w*w;
    float[] uix = u[ix];
    if (ix<b) {
      float x = (float)(b-ix);
      float xx = x*x;
      for (int iz=0; iz<b; ++iz) {
        float z = (float)(b-iz)*_dz;
        float rs = (xx+z*z)*ws;
        float e = exp(-rs);
        uix[iz     ] *= e;
        uix[nz-1-iz] *= e;
      }
      for (int iz=b; iz<nzmb; ++iz) {
        float r = x*w;
        uix[iz] *= exp(-r*r);
      } 
    } else if (ix<nx-b) {
      for (int iz=0; iz<b; ++iz) {
        float r = (float)(b-iz)*w;
        float e = exp(-r*r);
        uix[iz     ] *= e;
        uix[nz-1-iz] *= e;
      }
    } else {
      float x = 1.0f+(float)(b+ix-nx);
      float xx = x*x;
      for (int iz=0; iz<b; ++iz) {
        float z = (float)(b-iz);
        float rs = (xx+z*z)*ws;
        float e = exp(-rs);
        uix[iz     ] *= e;
        uix[nz-1-iz] *= e;
      }
      for (int iz=b; iz<nzmb; ++iz) {
        float r = x*w;
        uix[iz] *= exp(-r*r);
      }
    }
  }

  private static float[][] extendModel(float[][] c) {
    int nz = c[0].length;
    int nx = c.length;
    float[][] v = new float[nx+2*_nb][nz+2*_nb];

    for (int ix=0; ix<_nx;++ix){
      for (int iz=0; iz<_nz;++iz){
        v[ix+_nb][iz+_nb]=c[ix][iz];
      }
    }
    for (int ix=_nb; ix<nx+_nb; ++ix) {
      for (int iz=0, jz=nz+_nb; iz<_nb; ++iz, ++jz) {
        v[ix][iz] = v[ix][_nb];
        v[ix][jz] = v[ix][nz+_nb-1];
      }
    }
    for (int ix=0, jx=nx+_nb; ix<_nb; ++ix, ++jx) {
      copy(v[_nb],v[ix]);
      copy(v[nx+_nb-1],v[jx]);
    }
    return v;
  }




  private static void absorbClayton(
    float[][] um, float[][] ui, float[][] up)
  {
    // A1: ix = 0
    for (int iz=4, kx=4; iz<_nz+2*_nb-4; ++iz) {
      float a = sqrt(_v2[kx][iz])*(float)(_dt/_dx);
      up[kx][iz] = ui[kx][iz]+a*(ui[kx+1][iz]-ui[kx][iz]);
    }

    // A1: ix = nx-1
    for (int iz=4, kx=_nx+2*_nb-1-4; iz<_nz+2*_nb-4; ++iz) {
      float a = sqrt(_v2[kx][iz])*(float)(_dt/_dx);
      up[kx][iz] = ui[kx][iz]-a*(ui[kx][iz]-ui[kx-1][iz]);
    }

    // A1: iz = 0
    for (int ix=4, kz=4; ix<_nx+2*_nb-4; ++ix) {
      float a = sqrt(_v2[ix][kz])*(float)(_dt/_dz);
      up[ix][kz] = ui[ix][kz]+a*(ui[ix][kz+1]-ui[ix][kz]);
    }

    // A1: iz = nz-1
    for (int ix=4, kz=_nz-1+2*_nb-4; ix<_nx+2*_nb-4; ++ix) {
      float a = sqrt(_v2[ix][kz])*(float)(_dt/_dz);
      up[ix][kz] = ui[ix][kz]-a*(ui[ix][kz]-ui[ix][kz-1]);
    }
  }


}
