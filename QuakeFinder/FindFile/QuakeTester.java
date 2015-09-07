/*
Change log
----------
2015-07-17 : Testing release
2015-07-20 : Scoring function update
2015-07-26 : BUG FIX: Remove offset when loading other quakes.
2015-07-26 : Add (time in seconds from start) to algorithm for other quakes.
*/


import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.io.*;
import java.util.Arrays;
import java.util.*;
import java.text.*;
import java.nio.ByteBuffer;
import org.w3c.dom.Document;
import org.w3c.dom.*;


import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

public class QuakeTester
{
    public static String execCommand = null;
    public static boolean debug = true;
    public static boolean training = false;
    public static String dataFolder = "data";
    public static int seed = 0;
    public static Process solution;

    public void printMessage(String s) {
        if (debug) {
            System.out.println(s);
        }
    }

    public class Coordinate
    {
        double latitude;
        double longitude;

        public Coordinate(double lat, double lon) {
            latitude = lat;
            longitude = lon;
        }
    }

    public class Quake
    {
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

    public class DataSet
    {
        Coordinate[] sites;
        Quake[] quakes;
        int sampleRate;
        int numOfSites;
        int numOfEMA;
        int numOfQuakes;

        int[] rawData = null;
        byte[] result = null;
        double[] EMA = null;

        String gtfStartTime, gtfEQtime;
        double gtfMagnitude, gtfLatitude, gtfLongitude, gtfDistToEQ, gtfEQSec;
        int gtfEQHour, gtfSite;

        DocumentBuilderFactory docBuilderFactory = null;
        DocumentBuilder docBuilder = null;

        public void loadSiteInfo(String sXmlFile) throws Exception
        {
            // load site information
            Document doc = docBuilder.parse (new File(sXmlFile));
            doc.getDocumentElement ().normalize ();

            NodeList listOfSites = doc.getElementsByTagName("Site");
            numOfSites = listOfSites.getLength();
            printMessage("Total number of sites: " + numOfSites);

            sites = new Coordinate[numOfSites];
            for (int s=0;s<numOfSites;s++) {
                Element siteElement = (Element)listOfSites.item(s);
                if (s==0) {
                    sampleRate = Integer.parseInt( siteElement.getAttribute("sample_rate") );
                    printMessage("Sample Rate = "+sampleRate);
                }
                double lat = Double.parseDouble( siteElement.getAttribute("latitude") );
                double lon = Double.parseDouble( siteElement.getAttribute("longitude") );
                sites[s] = new Coordinate(lat, lon);
                printMessage("site "+s+": "+lat+","+lon);
            }
            // allocate memory for hourly data
            rawData = new int[numOfSites*3*3600*sampleRate];
            result = new byte[3600*sampleRate*4];
        }

        public void loadEarthMagneticActivity(String sXmlFile) throws Exception
        {
            // load earth magnetic activity
            Document doc = docBuilder.parse (new File(sXmlFile));
            doc.getDocumentElement ().normalize ();

            NodeList listOfEMA = doc.getElementsByTagName("kp_hr");
            numOfEMA = listOfEMA.getLength();
            printMessage("Total number of EM activities: " + numOfEMA);

            EMA = new double[numOfEMA];
            for (int i=0;i<numOfEMA;i++) {
                EMA[i] = Double.parseDouble(listOfEMA.item(i).getFirstChild().getNodeValue());
            }
        }

        public void loadOtherQuakes(String sXmlFile) throws Exception
        {
            // load earth magnetic activity
            Document doc = docBuilder.parse (new File(sXmlFile));
            doc.getDocumentElement ().normalize ();

            NodeList listOfQuakes = doc.getElementsByTagName("Quake");
            numOfQuakes = listOfQuakes.getLength();
            printMessage("Total number of other quakes: " + numOfQuakes);

            quakes = new Quake[numOfQuakes];
            for (int i=0;i<numOfQuakes;i++) {
                Element quakeElement = (Element)listOfQuakes.item(i);
                int secs = Integer.parseInt( quakeElement.getAttribute("secs") );
                double lat = Double.parseDouble( quakeElement.getAttribute("latitude") );
                double lon = Double.parseDouble( quakeElement.getAttribute("longitude") );
                double depth = Double.parseDouble( quakeElement.getAttribute("depth") );
                double mag = Double.parseDouble( quakeElement.getAttribute("magnitude") );
                quakes[i] = new Quake(secs, new Coordinate(lat, lon), depth, mag);
            }
        }

        public double[] getOtherQuakes(int hour)
        {
            int hStart = hour*3600;
            int hEnd = (hour+1)*3600;
            int numInHour = 0;
            for (int i=0;i<numOfQuakes;i++) {
                if (quakes[i].timeSecs>=hStart && quakes[i].timeSecs<hEnd) numInHour++;
            }
            double[] oQuake = new double[numInHour*5];
            int q = 0;
            for (int i=0;i<numOfQuakes;i++) {
                if (quakes[i].timeSecs>=hStart && quakes[i].timeSecs<hEnd) {
                    oQuake[q] = quakes[i].coord.latitude;
                    oQuake[q+1] = quakes[i].coord.longitude;
                    oQuake[q+2] = quakes[i].depth;
                    oQuake[q+3] = quakes[i].magnitude;
                    oQuake[q+4] = quakes[i].timeSecs;
                    q+=5;
                }
            }
            return oQuake;
        }

        public DataSet(String sFolder) throws Exception
        {
            docBuilderFactory = DocumentBuilderFactory.newInstance();
            docBuilder = docBuilderFactory.newDocumentBuilder();
            loadSiteInfo(sFolder+"SiteInfo.xml");
            loadEarthMagneticActivity(sFolder+"Kp.xml");
            loadOtherQuakes(sFolder+"Quakes.xml");
        }

        public void readGTF() throws Exception
        {
            BufferedReader br = new BufferedReader(new FileReader("gtf.csv"));
            int numOfCases = Integer.parseInt(br.readLine());
            for (int i=0;i<numOfCases;i++)
            {
                String s = br.readLine();
                String[] token = s.split(",");
                int setID = Integer.parseInt(token[0]);
                if (setID==seed)
                {
                    gtfStartTime = token[1];
                    gtfEQtime = token[2];
                    gtfMagnitude = Double.parseDouble(token[3]);
                    gtfLatitude = Double.parseDouble(token[4]);
                    gtfLongitude = Double.parseDouble(token[5]);
                    gtfSite = Integer.parseInt(token[6]);
                    gtfDistToEQ = Double.parseDouble(token[7]);
                    // Calculate number of hours till EQ
                    SimpleDateFormat ft = new SimpleDateFormat("yyyy-MM-dd hh:mm:ss");
                    Date date1 = ft.parse(gtfStartTime);
                    long startMSec = date1.getTime();
                    Date date2 = ft.parse(gtfEQtime);
                    long eqMSec = date2.getTime();
                    gtfEQSec = (eqMSec-startMSec)/1000;
                    gtfEQHour = (int)(gtfEQSec / (60*60));
                    printMessage("Quake happened at hour = " + gtfEQHour +" and second = " + gtfEQSec+ " at row "+(gtfEQSec-3600.0*gtfEQHour)*sampleRate);
                    break;
                }
            }
            br.close();
        }

        public int[] loadHour(String sFolder, int h) throws Exception
        {
            String fname = sFolder + "test" + seed + "_" + h + ".bin";
            printMessage("loading "+fname);
            int[] dt = new int[numOfSites*3*3600*sampleRate];
            File file = new File(fname);
            InputStream input = new BufferedInputStream(new FileInputStream(file));
            DataInputStream din = new DataInputStream(input);
            int prev = -1;
            int diff;
            for (int i=0;i<dt.length;i++)
            {
                int b1 = (int)din.readByte(); if (b1<0) b1+=256;
                if ((b1&3)==1) {
                    diff = (b1>>2)-(1<<5);
                    dt[i] = prev + diff;
                } else if ((b1&3)==2) {
                    int b2 = (int)din.readByte(); if (b2<0) b2+=256;
                    diff = ((b1+(b2<<8))>>2)-(1<<13);
                    dt[i] = prev + diff;
                } else if ((b1&3)==3) {
                    int b2 = (int)din.readByte(); if (b2<0) b2+=256;
                    int b3 = (int)din.readByte(); if (b3<0) b3+=256;
                    diff = ((b1+(b2<<8)+(b3<<16))>>2)-(1<<21);
                    dt[i] = prev + diff;
                } else
                {
                    int b2 = (int)din.readByte(); if (b2<0) b2+=256;
                    int b3 = (int)din.readByte(); if (b3<0) b3+=256;
                    int b4 = (int)din.readByte(); if (b4<0) b4+=256;
                    dt[i] = ((b1+(b2<<8)+(b3<<16)+(b4<<24))>>2)-1;
                }
                prev = dt[i];
            }
            input.close();
            return dt;
        }

    }



    public double doExec() throws Exception {

        try {
            solution = null;

            try {
                solution = Runtime.getRuntime().exec(execCommand);
            } catch (Exception e) {
                System.err.println("ERROR: Unable to execute your solution using the provided command: "
                        + execCommand + ".");
                return -1;
            }

            BufferedReader reader = new BufferedReader(new InputStreamReader(solution.getInputStream()));
            PrintWriter writer = new PrintWriter(solution.getOutputStream(), true);
            new ErrorStreamRedirector(solution.getErrorStream()).start();

             DataSet dset = new DataSet(dataFolder+seed+"/");
             // load ground truth data
             dset.readGTF();

             // pass site info to algo
             // int init(int sampleRate, int numOfSites, vector<double> sitesData,)
             writer.println(dset.sampleRate);
             writer.println(dset.numOfSites);
             double[] sitesData = new double[dset.numOfSites*2];
             for (int i=0;i<dset.numOfSites;i++) {
                sitesData[i*2] = dset.sites[i].latitude;
                sitesData[i*2+1] = dset.sites[i].longitude;
             }
             writer.println(sitesData.length);
             for (double v : sitesData) {
                writer.println(v);
             }
             writer.flush();
             int initRet = Integer.parseInt(reader.readLine());

             if (training) {
                // training enabled, pass gtf to algorithm
                writer.println(1);
                writer.println(dset.gtfSite);
                writer.println(dset.gtfEQHour);
                writer.println(dset.gtfLatitude);
                writer.println(dset.gtfLongitude);
                writer.println(dset.gtfMagnitude);
                writer.println(dset.gtfDistToEQ);
                writer.flush();
             } else
             {
                writer.println(0);
             }

             // pass hourly data to algorithm
             double score = 0.0;
             int hcnt = 0;
             for (int h=0;h<dset.gtfEQHour;h++) {
                int[] hourlyData = dset.loadHour(dataFolder+seed+"/", h);
                double[] otherQuakes = dset.getOtherQuakes(h);
                // vector<double> forecast(int hour, vector<int> data, double K, vector<double> globalQuakes)
                writer.println(h);
                writer.println(hourlyData.length);
                for (int i=0;i<hourlyData.length;i++) {
                    writer.println(hourlyData[i]);
                }
                writer.println(dset.EMA[h]);
                writer.println(otherQuakes.length);
                for (double v : otherQuakes) {
                    writer.println(v);
                }
                // read forecast matrix
                int numOfRet = Integer.parseInt(reader.readLine());

                if (numOfRet!=dset.numOfSites*(2160)) {
                    System.err.println("ERROR: The number of elements in your return is incorrect. You returned "+numOfRet+" it should be "+(dset.numOfSites*(2160)) );
                    return -1.0;
                }

                double[] ret = new double[numOfRet];
                for (int i=0;i<numOfRet;i++) {
                    ret[i] = Double.parseDouble(reader.readLine());
                }

                // score return
                if (h>=768)
                {
                    int idx = h*dset.numOfSites;
                    // normalize to sum of 1
                    double sum = 0;
                    for (int i=idx;i<numOfRet;i++) {
                        sum += ret[i];
                    }
                    if (Math.abs(sum)>1e-9) {
                        for (int i=idx;i<numOfRet;i++) {
                            ret[i] /= sum;
                        }
                    }
                    // score only from hour 768
                    //F = sizeof(NN) * (2 * G – Sum of squared values in NN) - 1
                    int sizeofNN = dset.numOfSites*(2160 - h);
                    int gtfIndex = (dset.gtfEQHour)*dset.numOfSites + dset.gtfSite;
                    double ssv = 0; // sum of squared values in NN
                    for (int i=idx;i<numOfRet;i++) {
                            ssv += ret[i]*ret[i];
                    }
                    double hourScore = (double)sizeofNN * (2.0 * ret[gtfIndex] - ssv) - 1.0;
                    score += hourScore;
                    hcnt++;
                    printMessage(hourScore + " " + score);
                }
             }
             writer.println(-1);
             writer.flush();

             score /= hcnt;
             return score;

        } catch (Exception e) {
            System.out.println("FAILURE: " + e.getMessage());
            e.printStackTrace();
        }
        return -1.0;
    }


    public static void main(String[] args) throws Exception {


       for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-folder")) {
                dataFolder = args[++i];
            } else if (args[i].equals("-exec")) {
                execCommand = args[++i];
            } else if (args[i].equals("-seed")) {
                seed = Integer.parseInt(args[++i]);
            } else if (args[i].equals("-train")) {
                training = true;
            } else if (args[i].equals("-silent")) {
                debug = false;
            } else {
                System.out.println("WARNING: unknown argument " + args[i] + ".");
            }
        }

        if (execCommand == null) {
            System.err.println("ERROR: You did not provide the command to execute your solution." +
                    " Please use -exec <command> for this.");
            System.exit(1);
        }

        try {
            if (dataFolder != null) {
                double score = new QuakeTester().doExec();
                System.out.println("Score = " + score);
            } else {
                System.out.println("WARNING: nothing to do for this combination of arguments.");
            }
        } catch (Exception e) {
            System.out.println("FAILURE: " + e.getMessage());
            e.printStackTrace();
        }
    }

    class ErrorStreamRedirector extends Thread {
        public BufferedReader reader;

        public ErrorStreamRedirector(InputStream is) {
            reader = new BufferedReader(new InputStreamReader(is));
        }

        public void run() {
            while (true) {
                String s;
                try {
                    s = reader.readLine();
                } catch (Exception e) {
                    // e.printStackTrace();
                    return;
                }
                if (s == null) {
                    break;
                }
                System.out.println(s);
            }
        }
    }
}



