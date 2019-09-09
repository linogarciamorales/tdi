PImage img_in;
double[][] im_in;

double[][] laplace = {{ 0,  1, 0 }, 
                     { 1, -4, 1 }, 
                     { 0,  1, 0 }};

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
  image(img_in, 0, 0, N, M);
  
  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in);
  
  double[][] im_gr = new double[N][M];
  filter(im_gr, im_in, laplace);
  PImage img_out1 = createImage(N, M, GRAY);
  convert2Pimage(img_out1, im_gr);
  image(img_out1, N+2, 0, N, M);

  double[][] im_out = new double[N][M];
  minus(im_out, im_in, im_gr);
  PImage img_out3 = createImage(N, M, GRAY);
  convert2Pimage(img_out3, im_out);
  image(img_out3, 2*N+4, 0, N, M);

  smooth();
  noLoop();
} 

void draw() {
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("LP_BW_Enh.png");    
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

void minus(double[][] im_out, double[][] im_in1, double[][] im_in2) {
  //out = in1-in2
  // Loop through every pixel in the image
  for (int x = 0; x < N; x++) {  
    for (int y = 0; y < M; y++) {
      im_out[x][y] = im_in1[x][y] - im_in2[x][y];
    }
  }
}
