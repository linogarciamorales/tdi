PImage img_in;
double[][] im_in;

int[] hist_in = new int[256];
int[] hist_out = new int[256];
int[] hist_cum = new int[256];
double[] tfunc = new double[256];

int N, M;
int[] limits = new int[2];
int histMin = 0;
int histMax = 0;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+ 2, M);
}

void setup() {
  colorMode(RGB);
  image(img_in, 0, 0, N, M);
  loadPixels();
  histo(hist_in, img_in);
  cumHisto(hist_cum, hist_in);
  normalization(tfunc, hist_cum);
  
  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in);
  double[][] im_out = new double[N][M];
  brightnessEQ(im_out, img_in, tfunc);
  PImage img_out = createImage(N, M, GRAY);
  convert2Pimage(img_out, im_out);
  image(img_out, N+2, 0, N, M);
  histo(hist_out, img_out);

  smooth();
  noLoop();
}

void draw() {
  drawHisto(hist_in, 0);
  //drawHisto(hist_cum, N+2);
  drawHisto(hist_out, N+2);
  //save("histo_detalle_BW_nor.png");
}

void histo(int[] hist, PImage img) {
  // Calculate the histogram
  for(int k = 0; k < img.pixels.length; k++) {
      int bright = int(brightness(img.pixels[k]));
      hist[bright]++;
  }
}

void drawHisto(int[] hist, int xdisp) {
  double histMax = max(hist);
  stroke(255);
  // Draw half of the histogram (skip every second value)
  for (int i = 0; i < N; i += 2) {
    // Map i (from 0..img.width) to a location in the histogram (0..255)
    int which = int(map(i, 0, N, 0, 255));
    // Convert the histogram value to a location between 
    // the bottom and the top of the picture
    int y = int(map((float)hist[which], 0, (float)histMax, M, 0));
    line(i+xdisp, M, i+xdisp, y);
  }
}

void cumHisto(int[] ohist, int[] ihist) {
  ohist[0] = ihist[0];
  for(int i = 1; i < ihist.length; i++) {
    ohist[i] = ihist[i] + ohist[i-1];
  }
}

void normalization(double[] tf, int[] ihist) {
  // normalization
  int P = max(ihist);
  int Q = min(ihist);
  for(int i = 0; i < ihist.length; i++) {
    tf[i] = (double)(ihist[i] - Q)/ (P - Q);
  }  
}

int maxIndex(int[] hist) {
  int pos = -1;
  int max = Integer.MIN_VALUE; //lowest possible value of an int.
  for(int i=0; i < hist.length; i++) {
    if(hist[i] > max) {
      pos = i;
      max = hist[i];
    }
  }
  return pos;
}

int minIndex(int[] hist) {
  int pos = -1;
  int min = Integer.MAX_VALUE; //lowest possible value of an int.
  for(int i = 0; i < hist.length; i++) {
    if(hist[i] < min) {
      pos = i;
      min = hist[i];
    }
  }
  return pos;
}

void brightnessEQ(double[][] im, PImage img, double[] tf) {
  // equalization
  for(int k = 0, x = 0, y = 0; k < img.pixels.length; k++, x = k % im.length, y = int(k / im.length)) {
    int bright = int(brightness(img.pixels[k]));
    im[x][y] = tf[bright];
  }
}

void constrainLimits(double[][] im) {
  for(int k=0,x=0,y=0; k < im.length*im[0].length; k++, x=k % im.length, y = int(k/im.length)) {
      if(im[x][y] > 1.0) {
        im[x][y] = 1.0;
      }
      if(im[x][y] < 0.0) {
        im[x][y] = 0.0;
      }
  }
}

void convert2double(double[][] im, PImage img) {
  int rows = im.length;
  for(int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col= int(k/ rows)) {
    im[row][col] = map(red(img.pixels[k]), 0, 255, 0, 1.0);
  }
}

void convert2Pimage(PImage img, double[][] im) {
  img.loadPixels();
  int rows = im.length;
  for (int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col = int(k / rows)) {
    double Y = im[row][col]* 255;
    //img.pixels[k] = color(Y);
    img.set(row, col, color((float)Y));
  }
  img.updatePixels();
}
