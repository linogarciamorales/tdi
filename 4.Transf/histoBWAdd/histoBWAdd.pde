PImage img_in;
int[] hist_in = new int[256];
int[] hist_out = new int[256];

int N, M;
int[] limits = new int[2];
int histMin = 0;
int histMax = 0;
int alpha = 20;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+ 2, M);
}

void setup() {
  img_in.loadPixels();
  image(img_in, 0, 0, N, M);
  hist_in = histo(img_in);
  PImage img_out = createImage(N, M, GRAY);
  brightnessAdd(img_out, img_in, alpha);
  image(img_out, N+2, 0, N, M);
  hist_out = histo(img_out);

  smooth();
  noLoop();
}

void draw() {
  drawHisto(hist_in, 0);
  drawHisto(hist_out, N+2);
  //save("histo_detalle_BW_add.png");
}

int[] histo(PImage im_in) {
  int[] hist = new int[256];

  // Calculate the histogram
  for(int k = 0; k < im_in.pixels.length; k++) {
      int bright = int(brightness(im_in.pixels[k]));
      hist[bright]++; 
  }
  return hist;
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

void brightnessAdd(PImage img_out, PImage img_in, int alpha) {
  for (int k = 0, x = 0, y = 0; k < img_in.pixels.length; x = k % N, y = int(k / N), k++) {
    img_out.set(x, y, color(brightness(img_in.pixels[k])+ alpha));
    //img_out.set(x, y, color(brightness(img_in.pixels[k])));
  }
}
