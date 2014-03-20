package lab7;

import edu.mines.jtk.util.Parallel;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;



/**
 * 1D ACOUSTIC WAVE EQUATION FINITE DIFFERENCE
 * 
 * 1  d^2      d  1  d
 * -- --- u  - -- -- -- u = s
 * pv dt^2     dx p  dx
 * 
 *
 * p = density (closest look to \rho is ascii characters)
 *  
 * @author Esteban Diaz
 * @version 2012.11.04
 *
 **/


public class Awefd1dp{

  // problem geometry:
  private float _dx, _dt;
  private int   _nx, _nt, _nxx; 
  private float _ox, _ot;

  private int _snap; // snapshot par: skip every _snap
                     // samples of the movie

  private int _ns,_nb;
  private float[] _v,_rho; // velocity and density

  private float[][] _source; // source function
  private float[]  _sourceX; // array with source coordinates
  private int [] _sourceXi; // coordinates converted to closest grid point

  private float[][] _movie; // wavefield movie



  // constructor:
  public Awefd1dp(float ox, float dx, int nx,
                 float ot, float dt, int nt,
                 float[] v, float[] rho,
                 float[][] source,float[] sourceX){
    this(ox,dx,nx,ot,dt,nt,v,rho,source,sourceX,1);
  } 
         
  public Awefd1dp(float ox, float dx, int nx,
                 float ot, float dt, int nt,
                 float[] v, float[] rho,
                 float[][] source,float[] sourceX, int snap){
    this(ox,dx,nx,ot,dt,nt,v,rho,source,sourceX,snap,1);
  } 
                 
  // constructor:
  public Awefd1dp(float ox, float dx, int nx,
                 float ot, float dt, int nt,
                 float[] v, float[] rho,
                 float[][] source,float[] sourceX,
                 int snap,int nb){

    _ox = ox; _dx = dx; _nx = nx; 
    _ot = ot; _dt = dt; _nt = nt; 
    _source = source;
    _sourceX = sourceX;
    _snap = snap;
    _ns = sourceX.length;
    _nb = nb; 
    _nxx = _nx + 2*_nb;

    _v = pad(v) ; 
    _rho = pad(rho);
    get_sourceXgrided();
    _movie = zerofloat(_nx,(int)(_nt));

    Check.argument(rho.length==v.length,"vel and rho must have same dim");
    Check.argument(rho.length==nx,"rho array length != nx");
  } 


  public float[][] apply(){
    fd1d_forward();
    return _movie;
  }



  /*
    Solves the 1d wave acoustic wave equation 
    forward in time.

   Note to myself:

   // solves right hand side:
   *          _              _
   *         |     d  1  d    |  
   * utmp =  | s + -- -- -- u | * pv^2*dt^2
   *         |     dx p  dx   |
   *         |_              _|
   *
   * Note to myself:
   *
   * algorithmic:
   * 1) solve:
   *   _           _
   *  |  d  1  d    |  
   *  |  -- -- -- u | 
   *  |  dx p  dx   |
   *  |_           _|
   *   
   * By:
   * utmp =  -G'(1/p)G u
   *  with G = forward difference operator
   *
   * 2) Inject source:
   *   utmp += source[is][it]   
   *
   * 3) scale:
   *   utmp *= pv^2*dt^2
   *
   * 4) solve forward time-steping:
   *   u[x][t+1] = utmp +2*u[x][t] -u[x][t-1]
   *
   */
  private void fd1d_forward(){

    // Snapshot temporarly arrays:
    float [] um1   = new float[_nxx];
    float [] uo   = new float[_nxx];
    float [] up1   = new float[_nxx];
    float [] utmp = new float[_nxx];

    float [] rhov2 = new float[_nxx]; // precompute v^2*dt^2
    float [] buoyancy = new float[_nxx];  // 
    
    for (int ix=0; ix<_nxx; ++ix){
      rhov2[ix] = _rho[ix]*_v[ix]*_v[ix]*_dt*_dt;
      buoyancy[ix] = 1.0f/(_rho[ix]); //buoyancy*1/dx^2
    }
    float [] up1_aux = zerofloat(_nxx);

    for (int it=0; it<_nt ; ++it){
      // spacial derivatives:
      utmp = div_rho_grad(uo,buoyancy); 
      // source injection:
      source_inject(it, utmp);
      // time step: 
      for (int ix=0; ix<_nxx; ++ix){
        up1[ix] = 2.0f*uo[ix] +rhov2[ix]*utmp[ix]-um1[ix];
      }
      // circulate arrays:
      for (int ix=0; ix<_nxx; ++ix){
        um1[ix] = uo[ix];
        uo[ix]  = up1[ix];
      }
      // send uo to movie:
      for (int ix=0; ix<_nx; ++ix)
        _movie[it][ix] = uo[ix+_nb];
    }
  }





  /*
   * 1) solve:
   *   _           _
   *  |  d  1  d    |  
   *  |  -- -- -- u | 
   *  |  dx p  dx   |
   *  |_           _|

    By simply doing:
    
    G'kGu

    G' is the transpose of G (transpose)
    k is buoyancy
    G is the backward difference operator.
 
  */
  private float[] div_rho_grad (float [] f, float [] k){
    
    int n1 = f.length;
    float [] u = new float[n1];

    float gi=0.0f;
    for (int i1=1; i1<n1;++i1){
      gi  = f[i1  ];  // gather 
      gi -= f[i1-1];  //  = g[i] = f[i] - f[i-1] //backward
      gi *= k[i1  ];    // scale: g[i] *= 1/rho[i] 
      u[i1-1] -= gi;
      u[i1  ]  = gi; // scatter  h[i] = g[i] - g[i+1] //forward
    }
    return mul(-1.0f/(_dx*_dx),u); // proper scaling
  }
  
  private void get_sourceXgrided(){
    // Nearest neighbor gridding
    _sourceXi = new int[_ns];
    for (int is=0; is<_ns; ++is){
      int ix =(int)(Math.round((_sourceX[is]-_ox)/_dx));
      _sourceXi[is]=ix+_nb;
    }
  }


  private void source_inject(int it,float[] ux){
    for (int is=0; is<_ns;++is){
      int ix = _sourceXi[is];
      ux[ix] += _source[is][it];
    }
  }

  private float[] pad(float[] u){
    float [] out = zerofloat(_nxx);

    for (int ix=0; ix<_nb; ++ix)
      out[ix] = u[0];

    for (int ix=_nb; ix<_nb+_nx; ++ix)
      out[ix] = u[ix-_nb];

    for (int ix=_nb+_nx; ix<_nxx; ++ix)
      out[ix] = u[_nx-1];

    return out;
  }

}
