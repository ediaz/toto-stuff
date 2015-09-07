import java.io.*;
import  java.util.Arrays;
import java.lang.Math;

public class Find{
  static int sampleRate, nsites;
  static Coordinate[] sitesCoor;
  static SiteData[] sitesData;
  static Quake[] gQuakes; 
  static SiteProb[] sitesProb; 
  //
 
  static int hour, oldhour;
  static int nchannels = 3;

  static FileOutputStream fw1 ;  
  static BufferedOutputStream bos ;
  static DataOutputStream dos; 




  public static  void main(String args[]){

    BufferedReader stdin = 
      new BufferedReader(new java.io.InputStreamReader(System.in));


    int sr = readInt(stdin); 
    int nSites = readInt(stdin);
    int sitesLen = readInt(stdin);
    double[] sData = readDouble(stdin,sitesLen);
    
    System.err.printf("entered init\n");
    int ret = init(sr,nSites,sData);
    System.err.printf("passed init\n");

    System.out.println(ret);


    System.err.printf("sample rate=%d\n",sampleRate);
    System.err.printf("number of sites=%d\n",nSites);
    System.err.printf("n elementents= %d\n",sitesLen);

    
    int doTraining = readInt(stdin);
    if (doTraining == 1){
      trainIt(stdin);
    }else{
      int pass = 0;
      while(pass == 0){
        System.err.printf("===========================\n");
        pass = goHour(stdin); 
        System.err.printf("Passed? %d\n",1-pass);
        System.err.printf("===========================\n");
      }
    }
  }


  public static int init(int srate, int ns, double[] sData){

  try{
    fw1 = new FileOutputStream("file1.txt");
    bos=new BufferedOutputStream(fw1);
    dos=new DataOutputStream(bos); 
  }catch(java.io.IOException e) {  
        System.out.println(e);   
  }


    sampleRate = srate;
    nsites = ns;
    sitesCoor = new Coordinate[nsites];
    sitesProb = new SiteProb[nsites];
    sitesData = new SiteData[nsites]; 
    for (int isite=0; isite<nsites ; ++isite){
      double lat = sData[isite*2];
      double lon = sData[isite*2+1];
      sitesCoor[isite] = new Coordinate(lat,lon);
      sitesProb[isite] = new SiteProb(2160);
    }
    return 0;
  }



  public static void trainIt(BufferedReader buff){
    int    gtf_Site, gtf_Hour;
    double  gtf_Latitude, gtf_Longitude, gtf_Magnitude,
                 gtf_DistToQuake;

    gtf_Site = readInt(buff);
    gtf_Hour = readInt(buff);
    gtf_Latitude = readDouble(buff);
    gtf_Longitude = readDouble(buff);
    gtf_Magnitude = readDouble(buff);
    gtf_DistToQuake = readDouble(buff);
  }

   
  public static int goHour(BufferedReader stdin){
    hour = readInt(stdin);
    System.err.printf("processing hour %d\n",hour);    
    
    if (hour==-1  ) return 1;
    int dlen = readInt(stdin);

    int[] data = readInt(stdin,dlen);

    double K = readDouble(stdin);
    System.err.printf("index K=%g\n",K);
    int qlen = readInt(stdin);
    double[] gQuakes = readDouble(stdin,qlen);

    double[] retM = forecast(hour,data,K,gQuakes);  
    System.out.println(retM.length);
    for (int i=0; i<retM.length; i++)
      System.out.println(retM[i]);
    System.out.flush();

    oldhour = hour; 
    return 0;
  }


  public static double[] forecast(int hour, int[] data, double K, double[] globalQuakes){

    initData(data,hour,globalQuakes);
    
    

    double[] retM = new double[2160*nsites];
    for (int isite=0; isite<nsites; ++isite){
      double[] prob = sitesProb[isite].updateProb(
                                  sitesData[isite],
                                  gQuakes,K, 
                                  sitesCoor[isite]); // bad design, each
                                    // site should have it's coordinate
      for(int h=0; h<2160; ++h){
        retM[h*nsites+isite] = prob[h];
      }
    }
    return retM;
  }


  public static double[] int2double(int[] in){

    double[] out = new double[in.length];
    for (int i=0; i<in.length; ++i)
      out[i] = (double)in[i];
    return out;
  }


  public static void initData(int[] Data, int h, double[] gQ){
    hour = h;
    

    int nsamples = sampleRate*3600; 
    for (int isite=0; isite <nsites; ++isite){
      double[] chan1 = int2double(Arrays.copyOfRange(Data,
                                       isite*nsamples*3+0*nsamples,
                                       isite*nsamples*3+0*nsamples+nsamples));
      double[] chan2 = int2double(Arrays.copyOfRange(Data,
                                       isite*nsamples*3+1*nsamples,
                                       isite*nsamples*3+1*nsamples+nsamples));
      double[] chan3 = int2double(Arrays.copyOfRange(Data,
                                       isite*nsamples*3+2*nsamples,
                                       isite*nsamples*3+2*nsamples+nsamples));
      sitesData[isite] = new SiteData(chan1,chan2,chan3,30.0,isite);
    }

    // now, do the quakes:
    int nQ = (int)(gQ.length/5);
    gQuakes = new Quake[nQ];

    for (int iquake=0; iquake<nQ; ++iquake){
      double lat   = gQ[iquake*5+0];
      double lon   = gQ[iquake*5+1];
      double depth = gQ[iquake*5+2];
      double mag   = gQ[iquake*5+3];
      int    time  = (int) gQ[iquake*5+4];
      gQuakes[iquake] = new Quake(time, new Coordinate(lat,lon), 
                                 depth, mag);
    }
  }




  // input buffer readers 
  public static double readDouble(BufferedReader buff){
    double number = -1;
    try{
       number = Double.parseDouble(buff.readLine());
    }catch(java.io.IOException e) {  
        System.out.println(e);   
    }
    return number; 
  }
  public static double[] readDouble(BufferedReader buff,int n){
    double[] array = new double[n];
    for (int i=0; i<n; ++i)
      array[i] = readDouble(buff);
    return array;
  }


  public static int readInt(BufferedReader buff){
    int number = -1;
    try{
       number = Integer.parseInt(buff.readLine());
    }catch(java.io.IOException e) {  
        System.out.println(e);   
    }
    return number; 
  }
  public static int[] readInt(BufferedReader buff,int n){
    int[] array = new int[n];
    for (int i=0; i<n; ++i)
      array[i] = readInt(buff);
    return array;
  }


  public static class SiteProb{
    double[] prior, posterior, prediction, decay;
    int hour;
    double lm, lmold,dec,dprior,dposterior,lq,lqold;

    public SiteProb(int size){
      this(size,(double) size/3.03);
    }

    public SiteProb(int size, double d){
      dec = d;
      prior = new double[size];
      posterior = new double[size];
      prediction = new double[size];
      for (int i=0; i<size; ++i)
        prior[i]=1.0;
      int hour = 0;
      decay = new double[size];
      lmold=1.0;
      dprior=1.0;
      lqold=1.0;
    }
    
    public double[] updateProb(SiteData sD, Quake[] quakes, double K, 
                                  Coordinate coor){
      updateMag(sD, K);
      updateQuake(quakes,coor);

      for (int i=0; i<posterior.length;++i){
        double mydec = dec*dposterior;
        double dd = Math.exp(-(i-hour)*(i-hour)/(2.0*mydec*mydec));
        posterior[i] = prior[i]*dposterior*decay[i];
        prior[i] = posterior[i];
        dprior = dposterior;
        try{
          dos.writeFloat((float)posterior[i]);
        }catch(java.io.IOException e) { 
          System.out.println(e);   
        }
      }
        try{
          dos.flush();
        }catch(java.io.IOException e) { 
          System.out.println(e);   
        }
      hour +=1;
      return posterior;
    }

    private void updateMag(SiteData sD,double K){
      
      lm = sD.likehoodMag()*Math.exp(-K*K/(1.0d)); // Ksd = 2.0
      dposterior = dprior*lm/lmold;
      lmold = lm;
    }


    private void updateQuake(Quake[] quakes,Coordinate coor){
      /* 
        likehood of receiving earthquake signal at this station
        from quake i is given by the time from quake location 
        to station (so it will more likely in the future).
        Also, the  magnitude should make it more likely.
       */
      int size = decay.length;
      int nquakes = quakes.length;
      double lq =0.0;
      double[] qprob = new double[size];
      double quake =0.0;
      for (int iquake=0; iquake<nquakes;++iquake){
        double dist = distFrom(quakes[iquake].coord,coor);
        if(dist <800.0) quake+= quakes[iquake].magnitude;
      }
      lq = 1.0+1.0e-5 -Math.exp(-quake*quake/(2.0));

      dposterior *= lq/lqold;
      lqold = lq;
    }


    public double timeFrom(Coordinate co1, Coordinate co2){

      double d = distFrom(co1,co2);
      System.err.printf("d=%g \n ",d);
      double highSlope = 60./122;
      double lowSlope = 53.3/120;
      double meanSlope = 60 * ( highSlope + lowSlope ) / 2;      
      double earthRadius = 6372.795;
      double earthCircumference = 2 * Math.PI * earthRadius; //%km
      double t= (d/earthCircumference)*360.0*meanSlope;
      System.err.printf("t=%g ms=%g\n",t,meanSlope);
      return t;
    }

    public double distFrom(Coordinate co1, Coordinate co2){
      double lat1 = co1.lat;
      double lat2 = co2.lat;
      double lng1 = co1.lon;
      double lng2 = co2.lon;

      double earthRadius = 6372.795;
      double dLat = Math.toRadians(lat2-lat1);
      double dLng = Math.toRadians(lng2-lng1);
      double sindLat = Math.sin(dLat / 2);
      double sindLng = Math.sin(dLng / 2);
      double a = Math.pow(sindLat, 2) + Math.pow(sindLng, 2)
                *Math.cos(Math.toRadians(lat1)) *
                 Math.cos(Math.toRadians(lat2));

      double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));
      double dist = earthRadius * c;
  
      return dist;
    }


  }


  // aux classes
  public static class Coordinate{
    double lat;
    double lon;

    public Coordinate(double lat1, double lon1) {
      lat = lat1;
      lon = lon1;
    }
  }

  public static class SiteData{
    public double[] magnitude;
    public double mean;
    public double std,stdprior;
    public int isite;


    public SiteData(double[] c1, double[] c2, double[] c3){
      this(c1,c2,c3,50.0d);
    }
    
    public SiteData(double[] c1, double[] c2, double[] c3, double stdvt){
      this(c1,c2,c3,30.0d,0);
    }

    
    public SiteData(double[] c1, double[] c2, double[] c3, double stdvt,int ksite){
      magnitude = new double[c1.length]; 
  
      for (int i=0; i<c1.length; ++i){
        magnitude[i] = Math.sqrt(c1[i]*c1[i] + c2[i]*c2[i]+c3[i]*c3[i]); 
      }
      checkForReset();
      stdprior = stdvt;
      isite = ksite;

    }


    private void checkForReset(){
      mean = meanF(magnitude);
      int resetSamples = 20000;
      double[] check = Arrays.copyOfRange(magnitude,0,resetSamples);
      double mincheck = min(check);
      double maxcheck = max(check);
      int pass =0;
      for (int i=0; i<resetSamples; ++i){
        double val = check[i]-mean;
        if(val >1e5 || val <-1e5) ++pass;
      }  
      if(pass > 5000){
        for (int i=0; i<resetSamples; ++i){
          magnitude[i]= mean;
        }
      }
      mean= meanF(magnitude);
    }

    private double meanF(double[] m) {
      double sum = 0;
      for (int i = 0; i < m.length; i++) {
        sum += m[i];
      }
      return sum / m.length;
    }
    
    private double min(double[]m){
      double x = 1.0e25;
      for (int i = 0; i < m.length; i++) {
        x = (m[i] <x)? m[i]:x; 
      }
      return x;
    }
    private double max(double[]m){
      double x = -1.0e25;
      for (int i = 0; i < m.length; i++) {
        x = (m[i] >x)? m[i]:x; 
      }
      return x;
    }



    public double likehoodMag(){
      double lmag = 0.0;
      double sum = 0.0;
  
      for (int i=0; i<magnitude.length; ++i){
        sum += (magnitude[i]-mean)*(magnitude[i]-mean);
      }

      double st = Math.sqrt(sum/magnitude.length);

      int count =0;
      for (int i=0; i<magnitude.length; ++i){
        if(magnitude[i] > mean+2*st || magnitude[i]<mean-2*st) count++;
      }
      System.err.printf("count site=%d = %d\n",isite,count);
      double p =  1-Math.exp(-(count*count)/(2*2000*2000.0));//stdprior*stdprior));
      if(p < 1e-5) p =1e-5;
      return p;
    }

  }

  public static class Quake{
    int timeSecs;
    double depth;
    double magnitude;
    Coordinate coord;
    public Quake(int tim, Coordinate loc, double dep, double mag) {
      coord = loc;
      depth = dep;
      magnitude = mag;
      timeSecs = tim;
    }
  }
}
