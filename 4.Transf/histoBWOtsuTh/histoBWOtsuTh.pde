PImage img_in;

int[] hist_in = new int[256];
int[] hist_out = new int[256];

int N, M;
int[] limits = new int[2];
int histMin = 0, histMax = 0;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+ 2, M);
}

void setup() {
  colorMode(GRAY);
  smooth();
  noLoop();
}

void draw() {
  code();
}

void keyPressed() {
  if (key == 'q' || key == 'Q') {
    exit();
  }
  if (key == 's' || key == 'S') {
    save("histoBWOtsuTh.png");    
  }
}

void code() {
  image(img_in, 0, 0, N, M);
  histo(hist_in, img_in);
  
  PImage img_out = createImage(N, M, GRAY);
  int th = otsu(hist_in);
  println(th);
  brightnessThr(img_out, img_in, th);
  image(img_out, N+2, 0, N, M);
  histo(hist_out, img_out);
  drawHisto(hist_in, 0);
  drawHisto(hist_out, N+2);
}

int otsu(int[] hist) {
  int level = 0;
  // total - number of pixels in the given image.
  double total = (double) sumhist(hist); 
  // OTSU automatic thresholding method
  double sumB = 0;
  double wB = 0;
  double maximum = 0.0;
  int[] v = linspace(0, 255, 256);
  double sum1 = (double) dot(v, hist);
  for (int k = 0; k < hist.length; k++) {
    wB = wB + (double) hist[k];
    double wF = total - wB;
    if (wB == 0 || wF == 0) {
      continue;
    }
    sumB = sumB+ (double) k* hist[k];
    double mF = (sum1- sumB)/ wF;
    double between = wB* wF* ((sumB/ wB)- mF)* ((sumB/ wB)- mF);
    if ( between >= maximum ) {
        level = k;
        maximum = between;
    }
  }
  return level;
}

int sumhist(int[] hist) {
  int sum = 0;
  for(int k = 0; k < hist.length; k++) {
      sum += hist[k];
  }
  return sum;
}

int dot(int[] a, int[] b) {
  int sum = 0;
  int K = a.length;
  for(int k=0; k < K; k++) {
    sum = sum + (a[k]* b[k]);
  }
  return sum;
}

int[] linspace(int min, int max, int points) {  
  int[] d = new int[points];  
  for (int k = 0; k < points; k++){  
    d[k] = min + k * (max - min) / (points - 1);  
  }  
  return d;  
} 

double[] linspace(double min, double max, int points) {  
  double[] d = new double[points];  
  for (int k = 0; k < points; k++){  
    d[k] = min + k * (max - min) / (points - 1);  
  }  
  return d;  
} 

void brightnessThr(PImage out, PImage in, int th) {
  // thresholding
  out.loadPixels();
  for (int k = 0; k < in.pixels.length; k++) {    
    if(red(in.pixels[k]) <= th) {
      out.pixels[k] = color(0);
    } else { 
      out.pixels[k] = color(255);
    }
  }
  out.updatePixels();
}
