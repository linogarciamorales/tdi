PImage img_in, img_out;

int N, M;
int[] limits = new int[2];
int histMin = 0;
int histMax = 0;
int alpha = 20;
int scale = 2;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  N = img_in.width/ scale; 
  M = img_in.height/ scale;
  size(N, M);
}

void setup() {
  image(img_in, 0, 0, N, M);
  loadPixels();
  int[] hist = histo(img_in);
  int l = getMin(img_in);
  int m = getMax(img_in);
  println("min: "+l+" max: "+ m);
  float mu = getMean(img_in);
  float mu2 = getMean(hist);
  println("media: "+ mu, mu2);
  float var = getVar(img_in);  
  float var2 = getVar(hist);  
  println("varianza: "+ var, var2);
  int e = getEnergy(hist);
  println("energia: "+ e);  
  float h = getEntropy(hist);
  println("entropia: "+ h);

  smooth();
  noLoop();
  updatePixels();
  drawHisto(hist, 0);
}

void draw() {
  //save("histo_detalle_BW_add.png");
}

int[] histo(PImage im_in) {
  int rows = im_in.height;
  int cols = im_in.width;
  int[] hist = new int[256];
  // Calculate the histogram
  //for(int k = 0; k < im_in.pixels.length; k++) {
  for(int k = 0; k < rows* cols; k++) {
      int bright = int(brightness(im_in.pixels[k]));
      hist[bright]++; 
  }
  return hist;
}

double[] histoNormalized(PImage im_in) {
  int rows = im_in.height;
  int cols = im_in.width;
  int[] hist = histo(im_in);
  int L = hist.length;
  double[] histo = new double[L];
  for(int k = 0; k < L; k++) {
      histo[k] = ((double) hist[k])/(rows*cols);
  }
  return histo;
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

int getMax(PImage img) {
  int max = Integer.MIN_VALUE; //lowest possible value of an int.
  for(int k = 0; k < img.pixels.length; k++) {
    int bright = int(brightness(img.pixels[k]));
    if(bright > max) {
      max = bright;
    }
  }
  return max;
}

int getMin(PImage img) {
  int min = Integer.MAX_VALUE; //lowest possible value of an int.
  for(int k = 0; k < img.pixels.length; k++) {
    int bright = int(brightness(img.pixels[k]));
    if(bright < min) {
      min = bright;
    }
  }
  return min;
}

// Digital Image Processing, Second Edition
//Burger and Burge, p. 50
float getMean(PImage img) {
  int rows = img.height;
  int cols = img.width;
  int MN = rows* cols;
  float mu = 0;
  for(int k = 0; k < MN; k++) {
      int bright = int(brightness(img.pixels[k]));
      mu = mu + bright; 
  }
  return mu/ MN;
}

float getMean(int[] histo) {
  int rows = M* scale;
  int cols = N* scale;
  int MN = rows* cols;
  float mu = 0;
  for(int k = 0; k < histo.length; k++) {
      float bright = k* histo[k];
      mu = mu + bright; 
  }
  return mu/ MN;
}

float getVar(PImage img) {
  int rows = img.height;
  int cols = img.width;
  int MN = rows* cols;
  float mu = getMean(img);
  float var = 0;
  for(int k = 0; k < MN; k++) {
      int bright = int(brightness(img.pixels[k]));
      var = var + (bright- mu)* (bright- mu); 
  }
  var = sqrt(var/ MN);
  return var;
}

float getVar(int[] histo) {
  int rows = M* scale;
  int cols = N* scale;
  int MN = rows* cols;
  float mu = getMean(histo);
  float var = 0;
  for(int k = 0; k < histo.length; k++) {
      float bright = (k-mu)*(k- mu)* histo[k];
      var = var + bright; 
  }
  var = sqrt(var/ MN);
  return var;
}

int getEnergy(int[] hist) {
  int e = 0;
  for(int k = 0; k < hist.length; k++) {
      int h = hist[k];
      e = e + h* h; 
  }
  return e;
}

float getEntropy(int[] hist) {
  float e = 0;
  for(int k = 0; k < hist.length; k++) {
      int h = hist[k];
      e = e + h* log(h+ Float.MIN_VALUE); 
  }
  return -e;
}
