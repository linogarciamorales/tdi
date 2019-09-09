PImage img_in;
double[][] im_in;

int[] hist_in = new int[256];
int[] hist_out = new int[256];

int N, M;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  //img_in.resize(img_in.width/ 2, img_in.height/ 2);
  img_in.resize(img_in.width/ 4, img_in.height/ 4);
  N = img_in.width; 
  M = img_in.height;
  size(3*N+ 4, 2*M+2);
}

void setup() {
  background(255);
  image(img_in, 0, M+2, N, M);
  //histo(img_in, hist_in);
  
  im_in = new double[N][M];
  convert2double(im_in, img_in);
  double[][] im_out = new double[N][M];

  double[][] kernel1 = {{ 1.0/9.0, 1.0/9.0, 1.0/9.0 }, 
                      {  1.0/9.0, 1.0/9.0, 1.0/9.0 }, 
                      {  1.0/9.0, 1.0/9.0, 1.0/9.0 }};
  double[][] kernel2 = {{ 1.0/10.0, 1.0/10.0, 1.0/10.0 }, 
                      {  1.0/10.0, 2.0/10.0, 1.0/10.0 }, 
                      {  1.0/10.0, 1.0/10.0, 1.0/10.0 }};
  double[][] kernel3 = {{ 1.0/16.0, 2.0/16.0, 1.0/16.0 }, 
                      {  2.0/16.0, 4.0/16.0, 2.0/16.0 }, 
                      {  1.0/16.0, 2.0/16.0, 1.0/16.0 }};
  double[][] kernel4 = {{ 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0 }, 
                        { 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0 }, 
                        { 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0 }, 
                        { 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0 }, 
                        { 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0 }};
  int W = 11;
  double[][] kernel5 = new double[W][W];
  for (int x = 0; x < W; x++) {
    for (int y = 0; y < W; y++) {
      kernel5[x][y] = 1.0/((double) W*W);
    }
  }
                      
  filter(im_out, im_in, kernel1);
  PImage img_out1 = createImage(N, M, GRAY);
  convert2Pimage(img_out1, im_out);
  image(img_out1, 0, 0, N, M);
  //histo(img_out1, hist_out);
  
//  filter(im_out, im_in, kernel2);
  filter(im_out, im_in, kernel4);
  PImage img_out2 = createImage(N, M, GRAY);
  convert2Pimage(img_out2, im_out);
  image(img_out2, N+2, 0, N, M);

//  filter(im_out, im_in, kernel3);
  filter(im_out, im_in, kernel5);
  PImage img_out3 = createImage(N, M, GRAY);
  convert2Pimage(img_out3, im_out);
  image(img_out3, 2*N+4, 0, N, M);

  smooth();
  noLoop();
} 

void draw() {
  //drawHisto(hist_in, 0);
  //drawHisto(hist_out, N+2);
  //save("LP_BW.png");
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
  int W = kernel.length;
  int W2 = (W-1)/2;
//  for (int x = 1; x < N-1; x++) {  // Skip left and right edges
//    for (int y = 1; y < M-1; y++) {   // Skip top and bottom edges
  for (int x = W2; x < N-W2; x++) {  // Skip left and right edges
    for (int y = W2; y < M-W2; y++) {   // Skip top and bottom edges

      double sum = 0; // Kernel sum for this pixel
//      for (int v = -1; v <= 1; v++) {
//        for (int u = -1; u <= 1; u++) {
      for (int v = -W2; v <= W2; v++) {
        for (int u = -W2; u <= W2; u++) {
          // Calculate the adjacent pixel for this kernel point
          double val = im_in[x+u][y+v];
          // Multiply adjacent pixels based on the kernel values
//          sum += kernel[u+1][v+1] * val;
          sum += kernel[u+W2][v+W2] * val;
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
