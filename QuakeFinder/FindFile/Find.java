import java.io.*;
import  java.util.Arrays;
import java.lang.Math;

public class Find{
  static int sampleRate, nsites;
  static Coordinate[] sitesCoor;
  static SiteData[] sitesData;
  static Quake[] gQuakes; 
  static ArrayOutputStream dataFile;  
  static ArrayOutputStream sitesFile;
  static ArrayOutputStream quakesFile; 
  static SiteProb[] sitesProb; 
  //
 
  static int hour, oldhour;
  static int nchannels = 3;




  public static  void main(String args[]){

    BufferedReader stdin = 
      new BufferedReader(new java.io.InputStreamReader(System.in));


    int sampleRate = readInt(stdin); 
    int nSites = readInt(stdin);
    int sitesLen = readInt(stdin);
    double[] sitesData = readDouble(stdin,sitesLen);
    
    System.err.printf("entered init\n");
    int ret = init(sampleRate,nSites,sitesData);
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

    sampleRate = srate;
    nsites = ns;
    sitesCoor = new Coordinate[nsites];
    sitesProb = new SiteProb[nsites];
    for (int isite=0; isite<nsites ; ++isite){
      double lat = sData[isite*2];
      double lon = sData[isite*2+1];
      sitesCoor[isite] = new Coordinate(lat,lon);
      sitesProb[isite] = new SiteProb(2160);
    }
    try{
      dataFile = new ArrayOutputStream("data.bin");
      sitesFile = new ArrayOutputStream("sites.bin");
      quakesFile = new ArrayOutputStream("quake.bin");
    }catch(java.io.IOException e) {  
      System.err.println(e);   
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
    System.err.printf("I am in go Hour\n");    
    hour = readInt(stdin);
    System.err.printf("processing hour %d\n",hour);    
    
    if (hour==-1 ) return 1;
    System.err.printf("reading dlen\n");
    int dlen = readInt(stdin);
    System.err.printf("dlen=%d\n reading data\n",dlen);

    double[] data = readDouble(stdin,dlen);
    System.err.printf("data read \n");

    double K = readDouble(stdin);
    System.err.printf("index k=%g\n",K);
    int qlen = readInt(stdin);
    double[] gQuakes = readDouble(stdin,qlen);

    System.err.printf("entering forecast\n");    
    double[] retM = forecast(hour,data,K,gQuakes);  
    System.err.printf("passed forecast\n");
    System.err.printf("processing %d hour\n",hour);
    System.err.printf("retM[0]=%g\n",retM[0]); 
    System.out.println(retM.length);
    for (int i=0; i<retM.length; i++)
      System.out.println(retM[i]);
    System.out.flush();

    oldhour = hour; 
    return 0;
  }


  public static double[] forecast(int hour, double[] data, double K, double[] globalQuakes){

    System.err.printf("entered init Data\n");
    initData(data,hour,globalQuakes);
    System.err.printf("passed init Data nsites=%d\n",nsites);
    
    

    double[] retM = new double[2160*nsites];
    for (int isite=0; isite<nsites; ++isite){
      System.err.printf("computing site probability\n");
      double[] prob = sitesProb[isite].updateProb(
                                  sitesData[isite],
                                  gQuakes,K, 
                                  sitesCoor[isite]); // bad desig, each
                                    // site should have it's coordinate
      System.err.printf("done site probability\n");
      for(int h=0; h<2160; ++h){
        retM[h+2160*isite] = prob[h];
      }
    }
    return retM;
  }





  public static void initData(double[] Data, int h, double[] gQ){
    hour = h;
    

    sitesData = new SiteData[nsites]; 
    int nsamples = sampleRate*3600; 
    System.err.printf("writing nsamp=%d 3 chans nsites=%d \n",nsamples,nsites);
    for (int isite=0; isite <nsites; ++isite){
      double[] chan1 = Arrays.copyOfRange(Data,
                                       isite*nsamples*3+0*nsamples,
                                       isite*nsamples*3+0*nsamples+nsamples);
      double[] chan2 = Arrays.copyOfRange(Data,
                                       isite*nsamples*3+1*nsamples,
                                       isite*nsamples*3+1*nsamples+nsamples);
      double[] chan3 = Arrays.copyOfRange(Data,
                                       isite*nsamples*3+2*nsamples,
                                       isite*nsamples*3+2*nsamples+nsamples);
      sitesData[isite] = new SiteData(chan1,chan2,chan3);
      writeFloat(sitesData[isite].magnitude,dataFile); 
    }
    try{
      dataFile.flush();
    }catch(java.io.IOException e) {  
        System.out.println(e);   
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
      System.err.printf("depth of quake %d is %g\n",iquake,depth);
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
    public static double[] prior, posterior, prediction, decay;
    public static int hour;

    public SiteProb(int size){
      this(size,(double) size/3.33);
    }

    public SiteProb(int size, double dec){
      prior = new double[size];
      posterior = new double[size];
      prediction = new double[size];
      for (int i=0; i<size; ++i)
        prior[i]=1.0;
      int hour = 0;
      decay = new double[size];
      for (int i=0; i<size; ++i)
        decay[i] = Math.exp(-i*i/(2*dec*dec));
    }
    
    public static double[] updateProb(SiteData sD, Quake[] quakes, double K, 
                                  Coordinate coor){
      System.err.printf("updating magnetic likehood\n");
      updateMag(sD, K);
      System.err.printf("updating magnetic likehood DONE\n");
      updateQuake(quakes,coor);
      for (int i=1; i<posterior.length;++i){
        posterior[i] = prior[i-1]*prediction[i];
        prior[i-1] = posterior[i]; // roll prob
      }
      prior[posterior.length-1] = prior[posterior.length-2];
      hour +=1;
      return prior;
    }

    private static void updateMag(SiteData sD,double K){
      double lm = sD.likehoodMag();//*Math.exp(-K*K/(2.0d)); // Ksd = 1.0
      System.err.printf("lm=%g\n",lm);
      // spray likehood along time with a gaussian prob
      for (int i=0; i< decay.length; ++i)
        prediction[i] = decay[i]*lm;
    }


    private static void updateQuake(Quake[] quakes,Coordinate coor){
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

      for (int iquake=0; iquake<nquakes;++iquake){
        double tiq = quakes[iquake].timeSecs 
                    + timeFrom(quakes[iquake].coord,coor);
        tiq = tiq/3600.0-hour-1;//s to h (time where earthq will occur at this
                                //station (relative to current hour)           
        for (int i=0; i<size; ++i){
          qprob[i] += quakes[iquake].magnitude*
                      Math.exp(-(i-tiq)*(i-tiq)/(800.));
        }        
      }
      double sum=0.0; 
      for (int i=0; i<size; ++i){
        sum+= qprob[i]*qprob[i];
      }        
      sum = Math.sqrt(sum);
      for (int i=0; i<size; ++i){
         prediction[i]*=qprob[i]/sum; // normalize so that sum qprob[i] =1;
      }        
    }


    public static double timeFrom(Coordinate co1, Coordinate co2){

      double d = distFrom(co1,co2);
      double highSlope = 60./122;
      double lowSlope = 53.3/120;
      double meanSlope = 60 * ( highSlope + lowSlope ) / 2;      
      double earthRadius = 6372.795;
      double earthCircumference = 2 * Math.PI * earthRadius; //%km
      return (d/earthCircumference)*360.0*meanSlope;
    }

    public static double distFrom(Coordinate co1, Coordinate co2){
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
    public static double[] magnitude;
    public static double mean;
    public static double std,stdprior;



    public SiteData(double[] c1, double[] c2, double[] c3){
      this(c1,c2,c3,4.0d);
    }
    
    
    public SiteData(double[] c1, double[] c2, double[] c3, double stdvt){
      magnitude = new double[c1.length]; 
  
      for (int i=0; i<c1.length; ++i){
        magnitude[i] = Math.sqrt(c1[i]*c1[i] + c2[i]*c2[i]+c3[i]*c3[i]); 
      }
      checkForReset();
      stdprior = stdvt;

    }


    private static void checkForReset(){
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
      for (int i=0; i<magnitude.length; ++i)
         magnitude[i] -= mean;
    }

    private static double meanF(double[] m) {
      double sum = 0;
      for (int i = 0; i < m.length; i++) {
        sum += m[i];
      }
      return sum / m.length;
    }
    
    private static double min(double[]m){
      double x = 1.0e25;
      for (int i = 0; i < m.length; i++) {
        x = (m[i] <x)? m[i]:x; 
      }
      return x;
    }
    private static double max(double[]m){
      double x = -1.0e25;
      for (int i = 0; i < m.length; i++) {
        x = (m[i] >x)? m[i]:x; 
      }
      return x;
    }



    public static double likehoodMag(){
      double lmag = 0.0;
      double sum = 0.0;
      for (int i=0; i<magnitude.length; ++i){
        sum += magnitude[i]*magnitude[i];
      }
      double st = Math.sqrt(sum)/magnitude.length;
      System.err.printf("st=%e\n",st);
      return Math.exp(-(2*stdprior*stdprior)/st);
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


  public static void writeFloat(double[] in, ArrayOutputStream stream){
    float[] out = new float[in.length];
    for (int i=0; i<in.length; ++i)
      out[i] = (float) in[i];

    try{
      stream.writeFloats(out);
    }catch(java.io.IOException e) {  
      System.err.println(e);   
    }
  }


}
