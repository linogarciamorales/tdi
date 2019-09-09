PImage img_in, img_out;

int N, M;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  //img_in = loadImage("DSC05633.png");
  //img_in = loadImage("DSC05634.png");
  //img_in = loadImage("DSC05635.png");
  //img_in = loadImage("DJI_0670.png");
  //img_in = loadImage("DJI_0672.png");
  //img_in = loadImage("DJI_0673.png");
  //img_in = loadImage("calle.png");
  //img_in = loadImage("imagenProyect.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  //img_in.resize(img_in.width/ 4, img_in.height/ 4);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+ 2, M);
}

void setup() {
  image(img_in, 0, 0, N, M);
  Stat(img_in);

  img_out = histoEQ(img_in); 
  image(img_out, N+2, 0, N, M);

  Stat(img_out);

  smooth();
  noLoop();
  updatePixels();
}

void draw() {
  double[] histn_in = new double[256];
  histn_in = histoNormalized(img_in);
  drawHisto(histn_in, 0);
  double[] histn_out = new double[256];
  histn_out = histoNormalized(img_out);
  drawHisto(histn_out, N+2);
  //save("histo_detalle_BW2.png");
}

int[] histo(PImage img) {
  img.loadPixels();
  int[] hist = new int[256];
  // Calculate the histogram
  for(int k = 0; k < img.pixels.length; k++) {
      int bright = int(brightness(img.pixels[k]));
      hist[bright]++; 
  }
  return hist;
}
/*
double[] histoNormalized(PImage img) {
  img.loadPixels();
  double[] histn = new double[256];
  int[] hist = new int[256];
  // Calculate the histogram
  int L = img.pixels.length;
  for(int k = 0; k < L; k++) {
      int bright = int(brightness(img.pixels[k]));
      hist[bright]++; 
  }
  int H = hist.length;
  for(int k = 0; k < H; k++) {
    histn[k] = ((double) hist[k])/L;
  }
  return histn;
}
*/

double[] histoNormalized(PImage img) {
  int L = img.pixels.length;
  int[] hist = histo(img);
  int H = hist.length;
  double[] histn = new double[H];
  for(int k = 0; k < H; k++) {
      histn[k] = ((double) hist[k])/L;
  }
  return histn;
}

void drawHisto(int[] hist, int xdisp) {
  int histMax = max(hist);
  stroke(255);
  // Draw half of the histogram (skip every second value)
  for (int i = 0; i < N; i += 2) {
    // Map i (from 0..img.width) to a location in the histogram (0..255)
    int which = int(map(i, 0, N, 0, 255));
    // Convert the histogram value to a location between 
    // the bottom and the top of the picture
    int y = int(map(hist[which], 0, histMax, M, 0));
    line(i+xdisp, M, i+xdisp, y);
  }
}

double map(double val, double s1, double s2, double t1, double t2) {
  return val * (t2- t1)/(s2 -s1) + t1;
}

void drawHisto(double[] hist, int xdisp) {
  double histMin = getMin(hist);
  double histMax = getMax(hist);
  stroke(255);
  // Draw half of the histogram (skip every second value)
  for (int i = 0; i < N; i += 2) {
    // Map i (from 0..img.width) to a location in the histogram (0..255)
    int which = int(map(i, 0, N, 0, 255));
    // Convert the histogram value to a location between 
    // the bottom and the top of the picture
    Long Y = Math.round(map(hist[which], histMin, histMax, M, 0));
    int y = Integer.valueOf(Y.intValue());
    line(i+xdisp, M, i+xdisp, y);
  }
}

int getMax(PImage img) {
  img.loadPixels();
  int max = Integer.MIN_VALUE; //lowest possible value of an int.
  for(int k = 0; k < img.pixels.length; k++) {
    int bright = int(brightness(img.pixels[k]));
    if(bright > max) {
      max = bright;
    }
  }
  return max;
}

int getMax(int[] hist) {
  for(int k = 255; k >= 0; k--) {
    int bright = hist[k];
    if(bright > 0) {
      return k;
    }
  }
  return 0;
}

double getMax(double[] v) {
  double max = Double.MIN_VALUE; //lowest possible value of an int.
  for(int k = 0; k < v.length; k++) {
    double val = v[k];
    if(val > max) {
      max = val;
    }
  }
  return max;
}

int getMin(PImage img) {
  img.loadPixels();
  int rows = img.height;
  int cols = img.width;
  int MN = rows* cols;
  int min = Integer.MAX_VALUE; //lowest possible value of an int.
  for(int k = 0; k < MN; k++) {
    int bright = int(brightness(img.pixels[k]));
    if(bright < min) {
      min = bright;
    }
  }
  return min;
}

int getMin(int[] hist) {
  for(int k = 0; k < hist.length; k++) {
    int bright = hist[k];
    if(bright > 0) {
      return k;
    }
  }
  return 255;
}

double getMin(double[] v) {
  double min = Double.MAX_VALUE; //lowest possible value of an int.
  for(int k = 0; k < v.length; k++) {
    double val = v[k];
    if(val < min) {
      min = val;
    }
  }
  return min;
}

// Digital Image Processing, Second Edition
//Burger and Burge, p. 50
float getMean(PImage img) {
  float mu = 0;
  img.loadPixels();
  for(int k = 0; k < img.pixels.length; k++) {
      int bright = int(brightness(img.pixels[k]));
      mu = mu + bright; 
  }
  mu = mu/ img.pixels.length;
  return mu;
}

double getMean(double[] v) {
  double mu = 0;
  for(int k = 0; k < v.length; k++) {
      double val = k* v[k];
      mu = mu + val; 
  }
  return mu;
}

float getStd(PImage img) {
  float mu = getMean(img);
  float var = 0;
  img.loadPixels();
  for(int k = 0; k < img.pixels.length; k++) {
      int bright = int(brightness(img.pixels[k]));
      var = var + (bright- mu)* (bright- mu); 
  }
  float std = sqrt(var/ img.pixels.length);
  return std;
}

double getStd(double[] v) {
  double mu = getMean(v);
  double var = 0;
  for(int k = 0; k < v.length; k++) {
      double val = Math.pow(k-mu, 2)* v[k];
      var = var + val; 
  }
  double std = Math.sqrt(var);
  return std;
}

double getSkew(double[] histn) {
  double mu = getMean(histn);
  double Sk = 0;
  for(int k = 0; k < histn.length; k++) {
      double val = Math.pow(k-mu, 3)* histn[k];
      Sk = Sk + val; 
  }
  return Sk/ Math.pow(255, 2);
}

int getEnergy(int[] hist) {
  int e = 0;
  for(int k = 0; k < hist.length; k++) {
      int h = hist[k];
      e = e + h* h; 
  }
  return e;
}

double getEnergy(double[] v) {
  double e = 0;
  for(int k = 0; k < v.length; k++) {
      double val = v[k];
      e = e + Math.pow(val, 2); 
  }
  return e;
}

double log10 (double x) {
  return (Math.log(x) / Math.log(10));
}

double log2 (double x) {
  return (Math.log(x) / Math.log(2));
}

double getEntropy(double[] v) {
  double e = 0;
  for(int k = 0; k < v.length; k++) {
      double val = v[k];
      e = e + val* log2(val+ Double.MIN_VALUE); 
  }
  return -e;
}

int[] cumHisto(int[] ihist) {
  int[] ohist = new int[ihist.length];
  ohist[0] = ihist[0];
  for(int i = 1; i < ihist.length; i++) {
    ohist[i] = ihist[i] + ohist[i-1];
  }
  return ohist;
}

double[] cumHisto(double[] ihist) {
  double[] ohist = new double[ihist.length];
  ohist[0] = ihist[0];
  for(int i = 1; i < ihist.length; i++) {
    ohist[i] = ihist[i] + ohist[i-1];
  }
  return ohist;
}

void normalization(double[] tf, int[] ihist) {
  // normalization
  int P = max(ihist);
  int Q = min(ihist);
  for(int i = 0; i < ihist.length; i++) {
    tf[i] = (double)(ihist[i] - Q)/ (P - Q);
  }  
}

void brightnessEQ(double[][] im, PImage img, double[] tf) {
  int rows = im.length;
  int cols = im[0].length;
  int MN = rows* cols;
  img.loadPixels();
  // equalization
  for(int k = 0, row = 0, col = 0; k < MN; k++, row = k % rows, col = int(k / rows)) {
    int bright = int(brightness(img.pixels[k]));
    im[row][col] = tf[bright];
  }
}

void convert2double(double[][] im, PImage img) {
  int rows = im.length;
  int cols = im[0].length;
  int MN = rows* cols;
  img.loadPixels();
  for(int k = 0, row = 0, col = 0; k < MN; k++, row= k% rows, col= int(k/ rows)) {
    im[row][col] = map(brightness(img.pixels[k]), 0, 255, 0, 1.0);
  }
}

void convert2Pimage(PImage img, double[][] im) {
//  img.loadPixels();
  int rows = im.length;
  int cols = im[0].length;
  int MN = rows* cols;
  for (int k = 0, row = 0, col = 0; k < MN; k++, row= k% rows, col = int(k / rows)) {
    double Y = im[row][col]* 255;
    //img.pixels[k] = color(Y);
    img.set(row, col, color((float)Y));
  }
  //img.updatePixels();
}

void Stat(PImage img) {
/*  
  int[] hist = new int[256];
  hist = histo(img);
  println("min: "+ getMin(hist)+ " max: "+ getMax(hist));
*/
  
  double[] histn = new double[256];
  histn = histoNormalized(img);
  int l = getMin(img);
  int m = getMax(img);
  println("min: "+l+" max: "+ m+" ("+(m-l)+")");
  
  double mu2 = getMean(img);
  double mu = getMean(histn);
  println("media: "+ mu+ " ("+ mu2+")");
  float s2 = getStd(img);  
  double s = getStd(histn);  
  println("std: "+ s+ " ("+ s2+")");
  
  double S = 1- 1/(1+(pow((float)s,2)/pow(255,2)));
  println("smoothness: "+ S);
  
  double Sk = getSkew(histn);
  println("skew: "+ Sk);
  
  double e = getEnergy(histn);
  println("energia: "+ e);
  
  double h = getEntropy(histn);
  println("entropia: "+ h);
  println();
}

PImage histoEQ(PImage imgi) {
  PImage imgo = createImage(imgi.width, imgi.height, GRAY);
 
  double[] hist_in = histoNormalized(imgi);
  double[] hist_cum = cumHisto(hist_in);
  
  double[][] im_in = new double[N][M];
  convert2double(im_in, imgi);
  double[][] im_out = new double[N][M];
  brightnessEQ(im_out, imgi, hist_cum);
  convert2Pimage(imgo, im_out);
  return imgo;
}
