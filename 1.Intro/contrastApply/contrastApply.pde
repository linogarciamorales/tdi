PImage img_in;

int[] hist_in = new int[256];
int[] hist_out = new int[256];

int N, M;
int[] limits = new int[2];

void settings() {
  // Make a new instance of a PImage by loading an image file
  //img_in = loadImage("../../img/patt3.png");
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+ 2, M);
}

void setup() {
  colorMode(GRAY);
  image(img_in, 0, 0, N, M);
  histo(hist_in, img_in);
  
  float[] range = {0, 250};
  PImage img_out = createImage(N, M, GRAY);
  setLinearHisto(img_out, img_in, range);
  
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

float getMax(PImage in) {
  float Imax = 0.0;
  for (int k = 0; k < in.pixels.length; k++) {
    float i = brightness(in.pixels[k]);
    if(i > Imax) {
      Imax = i;
    }
  }
  return Imax;
}

float getMin(PImage in) {
  float Imin = 255;
  for (int k = 0; k < in.pixels.length; k++) {
    float i = brightness(in.pixels[k]);
    if(i < Imin) {
      Imin = i;
    }
  }
  return Imin;
}

void setLinearHisto(PImage out, PImage in, float[] range) {
  loadPixels();
  float s1 = getMin(in), s2 = getMax(in);
  float t1 = range[0], t2 = range[1];
  println(s1, s2, t1, t2);
  float a = (t2 - t1)/(s2 - s1);
  //float b = t1 - s1* a;
  for (int k = 0, row = 0, col = 0; k < in.pixels.length; k++, row= k% N, col = int(k / N)) {
    float i = brightness(in.pixels[k]);
    float Y = 0;
    if(i < s1) { Y = t1; }
    if(i > s2) { Y = t2; }
    if (i >= s1 && i <= s2) {
      //Y = a* i+ b;
      Y = (i-s1)* a + t1;
    }
    out.set(row, col, color(Y));
  }
}
