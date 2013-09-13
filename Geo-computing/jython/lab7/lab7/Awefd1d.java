package lab7;

import edu.mines.jtk.util.Parallel;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;



/**
 * 1D ACOUSTIC WAVE EQUATION FINITE DIFFERENCE
 * 
 *   1    d^2      d   1   d
 * -----  --- u  - -- ---  -- u = s
 * rho v  dt^2     dx rho  dx
 * 
 *
 * rho = density 
 *  
 * @author Esteban Diaz
 * @version 2012.11.04
 *
 **/


public class Awefd1d{


  // problem geometry:
  private float _dx, _dt;
  private int   _nx, _nt; 
  private float _ox, _ot;

  private int _snap; // snapshot par: skip every _snap
                     // samples of the movie

  private int _ns;
  private float[] _v,_rho; // velocity and density

  private float[][] _source; // source function
  private float[]  _sourceX; // array with source coordinates
  private int [] _sourceXi; // coordinates converted to closest grid point

  private float[][] _movie; // wavefield movie



  // constructor:
  public Awefd1d(float ox, float dx, int nx,
                 float ot, float dt, int nt,
                 float[] v, float[] rho,
                 float[][] source,float[] sourceX){
    this(ox,dx,nx,ot,dt,nt,v,rho,source,sourceX,1);
  } 
                 
  // constructor:
  public Awefd1d(float ox, float dx, int nx,
                 float ot, float dt, int nt,
                 float[] v, float[] rho,
                 float[][] source,float[] sourceX,
                 int snap){

    _ox = ox; _dx = dx; _nx = nx; 
    _ot = ot; _dt = dt; _nt = nt; 
    _v = v ; 
    _rho = rho;
    _source = source;
    _sourceX = sourceX;
    _snap = snap;
    _ns = sourceX.length;
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
    float [] um1   = new float[_nx];
    float [] uo   = new float[_nx];
    float [] up1   = new float[_nx];
    float [] utmp = new float[_nx];

    float [] rhov2 = new float[_nx]; // precompute v^2*dt^2
    float [] buoyancy = new float[_nx];  // 
    
    for (int ix=0; ix<_nx; ++ix){
      rhov2[ix] = _rho[ix]*_v[ix]*_v[ix]*_dt*_dt;
    }
    for (int ix=0; ix<_nx-1; ++ix){
      buoyancy[ix] = 2.0f/(_rho[ix]+_rho[ix+1]); //buoyancy*1/dx^2
    }
    buoyancy[_nx-1] =1/_rho[_nx-1];


    float [] up1_aux = zerofloat(_nx);

    for (int it=0; it<_nt ; ++it){
      // spacial derivatives:
      utmp = div_rho_grad(uo,buoyancy); 
      // source injection:
      source_inject(it, utmp);
      // time step: 
      up1 = sub(add(mul(2.0f,uo),mul(rhov2,utmp)),um1);

      // circulate arrays:
      um1 = uo;
      uo  = up1; 
      // send uo to movie:
      _movie[it] = uo;       
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
    for (int i1=1; i1<n1-1;++i1){
      gi  = f[i1  ];  // gather 
      gi -= f[i1-1];  //  = g[i] = f[i] - f[i-1] //backward
      gi *= k[i1  ];    // scale: g[i] *= 1/rho[i] 
      u[i1-1] -= gi;
      u[i1  ]  = gi; // scatter  h[i] = g[i] - g[i+1] //forward
    }
    return mul(-1.0f/(_dx*_dx),u); // proper scaling
  }
  



  private float[] div_rho_grad_old (float [] f, float [] k){
    int n1 = f.length;
    float [] u = new float[n1];
  
    d_dx(f,u);
    mul(k,u,u);
    d_dx(u,u); 
    return u; 
  }

  private void d_dx (float[] u, float[] du ){
    float cm2 = -1.0f/(12.0f*_dx);
    float cm1 = +8.0f/(12.0f*_dx);
    float cp1 = -cm1;
    float cp2 = cm2; 

    for (int ix=2 ; ix<_nx-2; ++ix){
      du[ix] = cm2*u[ix-2] + cm1*u[ix-1] + cp1*u[ix+1] + cp2*u[ix+2];
    }

  }





  private void get_sourceXgrided(){

    // Nearest neighbor gridding
    _sourceXi = new int[_ns];
    for (int is=0; is<_ns; ++is){
      int ix =(int)(Math.round((_sourceX[is]-_ox)/_dx));
      _sourceXi[is]=ix;
    }
  }


  private void source_inject(int it,float[] ux){
    for (int ix=0; ix<_nx;++ix){
      ux[ix] += _source[it][ix];
    }
  }

}
