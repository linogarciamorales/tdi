PImage img_in;
double[][] im_in;

double[][] hsobel = {{ 0.25, 0, -0.25 }, 
                    { 0.5, 0, -0.5 }, 
                    { 0.25, 0, -0.25 }};
double[][] vsobel = {{ -0.25, -0.5, -0.25 }, 
                    { 0, 0, 0 }, 
                    { 0.25, 0.5, 0.25 }};

int N, M;
int W;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(3*N+ 4, M);
}

void setup() {
  background(255);
  
  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in);
  
  double[][] im_h = new double[N][M];
  filter(im_h, im_in, hsobel);
  PImage img_out1 = createImage(N, M, GRAY);
  convert2Pimage(img_out1, im_h);
  image(img_out1, 0, 0, N, M);

  double[][] im_v = new double[N][M];
  filter(im_v, im_in, vsobel);
  PImage img_out2 = createImage(N, M, GRAY);
  convert2Pimage(img_out2, im_v);
  image(img_out2, N+2, 0, N, M);

  double[][] mag = new double[N][M];
  magnitude(mag, im_h, im_v);
  PImage img_out3 = createImage(N, M, GRAY);
  convert2Pimage(img_out3, mag);
  image(img_out3, 2*N+4, 0, N, M);

  smooth();
  noLoop();
} 

void draw() {
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("LP_BW_Sobel.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

void histo(PImage img, int[] hist) {
  // Calculate the histogram
  for(int k = 0; k < img.pixels.length; k++) {
      int bright = int(brightness(img.pixels[k]));
      hist[bright]++;
  }
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

void filter(double[][] im_out, double[][] im_in, double[][] kernel) {
  // Loop through every pixel in the image
  int q = floor(kernel.length/2.0);
  for (int x = q; x < N-q; x++) {  // Skip left and right edges
    for (int y = q; y < M-q; y++) {   // Skip top and bottom edges

      double sum = 0; // Kernel sum for this pixel
      for (int v = -q; v <= q; v++) {
        for (int u = -q; u <= q; u++) {
          // Calculate the adjacent pixel for this kernel point
          double val = im_in[x+u][y+v];
          // Multiply adjacent pixels based on the kernel values
          sum += kernel[q-u][q-v] * val;
        }
      }
      // For this pixel in the new image, set the gray value
      // based on the sum from the kernel
      im_out[x][y] = sum;
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
  //img.loadPixels();
  int rows = im.length;
  for (int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col = int(k / rows)) {
    double Y = im[row][col]* 255;
    //img.pixels[k] = color(Y);
    img.set(row, col, color((float)Y));
  }
  //img.updatePixels();
}

void magnitude(double[][] mag, double[][] im_h, double[][] im_v) {
  for (int x=0; x < N; x++) {
    for (int y=0; y < M; y++) { 
      mag[x][y] = sqrt((float)(im_h[x][y]*im_h[x][y] + im_v[x][y]*im_v[x][y]));
    }
  }
}
