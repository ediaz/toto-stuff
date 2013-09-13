package fdmod; 

import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Parallel;
import edu.mines.jtk.util.Stopwatch;
import edu.mines.jtk.dsp.RecursiveExponentialFilter;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Simple 2D muting functions; 
 * @author Esteban Diaz
 * @version 2013.06.08
 *


*/
public class Rtm{

  /*
    Common variables:
  */

  private static float _dx, _dz, _dt;
  private static float _ox, _oz, _ot;
  private static int   _nx, _nz, _nt;
  private static int   _ntout;

  private static float[][][] _src_coo; // [ns][2] (x,z)
  private static int _ns;
  private static float[] _wav; // [ns][nt]

  private static float[][][] _rvr_coo; // [ns][nr][2] (x,z)
  private static int _nr;
  private static float[][][] _dat; // [ns][nr][nt]

  private static int _nb;

  private static Fdmod2d modeling;

  public static float[][] _image;
  public static float[][] _sillu;
  private static boolean silluflag;

  public static int jsnap = 10;
  public static int jshot = 1, is0=0;
  private static boolean _verb=false;
  private static boolean _setsources=false;

  // extended image space-lag gathers variables:
  private static boolean _extended=false;
  private static int _fg, _jg, _ng, _nh;
  private static float[][][] _ximage;


  /**
   * RTM object constructor, pretty standard stuff. 
   * @param ox origin in the x axis of the domain
   * @param dx distance between grid points along x axis
   * @param nx number of grid points for x axis
   * @param oz origin in the z axis of the domain
   * @param dz distance between gridpoints along z axis
   * @param nz number of grid points along z axis
   * @param vel 2-D array containing the velocity model
   * @param nb number of grid points for the absorving region 43 works good
   *           on the Marmousi model
   * @param dat input data to be migrated
   * @param wav user supplied wavelet to model the source wavefield
   * @param sillu source illumination flag.
   */
  public Rtm(float ox, float dx, int nx,
             float oz, float dz, int nz,
             float ot, float dt, int nt,
             float[][] vel, // velocity array
             int nb,
             float[][][] dat,
             float[] wav, boolean sillu){

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

    _image = new float[nx][nz];
    silluflag = sillu;
    if (sillu)
      _sillu = new float[nx][nz];

    _wav = wav;
    _dat = dat;
    modeling = new Fdmod2d(ox,dx,nx,oz,dx,nz,ot,dt,nt,vel,nb,false);
  }


  /**
   * sets the sources geometry, the outer dimension is the number
   * of shots, the fast dimension =2: [0]: coordinate x, [1] coordinate z
   * @param src_coo source coordinates array: float[ns][1][2]
   */
  public static void setSources(float[][][] src_coo){
    _src_coo = src_coo;
    _setsources = true;
  }

  /**
   * fixed spread acquisition (receiver are the same for all shots)
   * method set up.
   * Sets the receivers geometry, the outer dimension is the number
   * of shots, the fast dimension =2: [0]: coordinate x, [1] coordinate z.
   * @param rvr_coo receiver coordinates array float[nr][2]
   */
  public static void setReceivers(float[][] rvr_coo){
    Check.argument(_setsources,
        "for fixed spread geometry, you must"+
        " first call the setSources() method");   
  
    _ns = _src_coo.length;
    _nr = rvr_coo.length;
    _rvr_coo = new float[_ns][_nr][2];

    for (int is=0; is<_ns; ++is)
      _rvr_coo[is] = rvr_coo;
  }

  /**
   * Variable receivers acquisition (receiver coordinates change for every
   * shot), number of receivers is the same (kind of thinking in streamer 
   * data.
   * sets the receivers geometry, the outer dimension is the number
   * of shots, the fast dimension =2: [0]: coordinate x, [1] coordinate z.
   * @param rvr_coo receiver coordinates array float[nr][2]
   */
  public static void setReceivers(float[][][] rvr_coo){
    _rvr_coo = rvr_coo;
  }

  /**
   * this method sets the jump in the shots for migration
   * @param jshot jump in the shot migration loop
   */
  public static void setShotJump(int jsht){
    jshot = jsht;
  }

  /**
   * this method sets the jump for wavefield snapshots correlation
   * the smaller jsnap the bigger the memory requirements
   * @param  jsnap jump in the output wavefields for imaging
   */
  public static void setJsnap(int j){
    jsnap = j;
  }

  /**
   * this method sets the first shot to migrate, this together with
   * setJsnap can migrate different parts of the data.
   * It can be ventageous to set setJsnap(10), and then run different
   * setFirst shot in different machines, and later combine the 
   * resulting images, thus creating and out of core distributed 
   * parallelization.
   * @param  firstshot  first shot to migrate (default=0)
   */
  public static void setInitialShot(int firstshot){
    is0 = firstshot;
  }

  /**
   * if called the verbosity flag will be turned on
   */
  public static void setVerbosity(){
    _verb = true;
  }




  /**
   * Once the object is initialized, a call to this method
   * will migrate ns/jshot shots
   */
  public static void goImage(){
    Stopwatch sw = new Stopwatch();
    sw.start();
    for (int is= is0; is< _ns; is+=jshot){  
      float[][][] swfl = goSwfl(is);
      float[][][] rwfl = goRwfl(is); 
      goImage(swfl,rwfl);
      if (_extended) // does xgathers
        goImageX(swfl,rwfl);
      if(_verb) 
        System.out.printf("source %d out of %d done in %10.5g sec%n",
                          is+1,_ns/jshot,sw.time());
      sw.restart();
    }
  }

  public static void goImage(int is){
    float[][][] swfl = goSwfl(is);
    float[][][] rwfl = goRwfl(is); 
    goImage(swfl,rwfl);
  }


  /**
   * this method migrates one shot, given by the user input sequential
   * index of the shots array
   * @param is the shot to migrate
   */
  public static void goImage(final float[][][] swfl, 
                             final float[][][] rwfl){
    final int nts = swfl.length;
    Parallel.loop(_nx,new Parallel.LoopInt() {
    public void compute(int ix) {  
      for (int its=0; its<nts; ++its){
        for (int iz=0; iz<_nz; ++iz){
          _image[ix][iz] += rwfl[nts-1-its][ix][iz]*swfl[its][ix][iz];
        }
      }  
    }});
  }
  /**
   * this method returns the migrated image, it must be called
   * after goImage() method, otherwise it will return an array of zeros
   */
  public static float[][] getImage(){
    return _image;
  }


  public static void setExtendedImage(int fg, int jg, int ng, int nh){
    _fg = fg;
    _jg = jg;
    _ng = ng;
    if (fg+ng*jg>_nx)
      _ng = (_nx-fg)/jg;
    _nh = nh;
    _ximage = new float[_ng][_nz][2*_nh+1];
    System.out.printf("n3=%d  o3=%f d3=%f%n",ng,fg*_dx-_ox,_dx*jg);
    System.out.printf("n2=%d  o2=%f d2=%f%n",_nz,_oz,_dz);
    System.out.printf("n1=%d  o1=%f d1=%f%n",1+2*_nh,-nh*2*_dx,2*_dx);
    
    _extended = true;
  }



  public static void goImageX(int is){
    float[][][] swfl = goSwfl(is);
    float[][][] rwfl = goRwfl(is); 
    goImageX(swfl,rwfl);
  }


  public static void goImageX(final float[][][] swfl, 
                              final float[][][] rwfl){
    final int nts = swfl.length;
    for (int ig=0; ig<_ng; ++ig){
      int ixg = _fg+ig*_jg;
      int lg0 = min(ixg,_nh);
      int lg1 = min(_nx-1-ixg,_nh);
      int maxlag = min(lg0,lg1);      
      //System.out.printf("igx=%d  %n",ixg,maxlag,ig);
      for (int its=0; its<nts; ++its){
        for (int iz=0; iz<_nz; ++iz){
          for (int ih=-maxlag; ih<=maxlag; ++ih){
            
            _ximage[ig][iz][ih+_nh] += swfl[its][ixg-ih][iz]*
                                       rwfl[nts-1-its][ixg+ih][iz];
          }
        }
      }  
    }
  }



  /**
   */
  public static float[][][] getImageX(){
    return _ximage;
  }





  // private methods:
  // computes the source wavefield
  private static float [][][] goSwfl(int is){
    float[][] wav_tmp = new float[1][_nt];
    wav_tmp[0] = _wav;

    modeling.setSources(_src_coo[is], wav_tmp);
    modeling.setReceivers(_src_coo[is]);
    modeling.goWavefield(jsnap);
    modeling.forwardAll();
    return modeling.getWavefield();
  }

  // computes the receiver wavefield
  private static float [][][] goRwfl(int is){

    modeling.setSources(_rvr_coo[is],_dat[is]);
    modeling.setReceivers(_src_coo[is]);
    modeling.setReverse();
    modeling.goWavefield(jsnap);
    modeling.forwardAll();
    return modeling.getWavefield();
  }

}
