PImage img_in;
double[][] im_in;

int N, M;
int W;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+ 2, M);
}

void setup() {
  background(255);
  image(img_in, 0, 0, N, M);
  
  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in);

  double sigma = 7.0/6.0;
  W = 5;
  double[][] kernel = new double[W][W];
  gaussianKernel2D(kernel, sigma);

  double[][] mu = new double[N][M];
  filter(mu, im_in, kernel);

  double[][] dev = new double[N][M];
  local_deviation(dev, im_in, mu, kernel);

  double[][] im_out = new double[N][M];
  MSCN(im_out, mu, dev, im_in);

  PImage img_out = createImage(N, M, GRAY);
  convert2Pimage(img_out, im_out);
  image(img_out, N+2, 0, N, M);

  smooth();
  noLoop();
} 

void draw() {
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("LP_BW_Gauss.png");
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
  img.loadPixels();
  int rows = im.length;
  for(int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col= int(k/ rows)) {
    im[row][col] = map(brightness(img.pixels[k]), 0, 255, 0, 1.0);
  }
}

void convert2Pimage(PImage img, double[][] im) {
  int rows = im.length;
  for (int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col = int(k / rows)) {
    double Y = im[row][col]* 255;
    //img.pixels[k] = color(Y);
    img.set(row, col, color((float)Y));
  }
  img.updatePixels();
}

void gaussianKernel2D(double[][] kernel, double sigma) {
  int W = kernel.length;
  int q = floor(kernel.length/2.0);
  double sum = 0;
  for (int v=-q, y=0; v <= q; v++, y=v+q) {
    for (int u=-q, x=0; u <= q; u++, x=u+q) {
      double A = 1.0/(2.0* PI* sigma*sigma);
      //double A = 1.0;
      kernel[x][y] = A* exp((float)(-(u*u+v*v)/(2*sigma*sigma)));
      sum += kernel[x][y];
    }
  }
  
  for(int x = 0; x < W; x++) {
    for(int y = 0; y < W; y++) {
      kernel[x][y] /= sum;
    }
  }
}

void local_mean(double[][] im_out, double[][] im_in, double[][] kernel) {
  filter(im_out, im_in, kernel);
}

void local_deviation(double[][] im_out, double[][] im_in, double[][] local_mean, double[][] kernel) {
  int rows = im_in.length;
  int cols = im_in[0].length;
  double[][] sigma = new double[rows][cols];
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      sigma[row][col] = im_in[row][col]* im_in[row][col];
    }
  }
  double[][] local_sigma = new double[rows][cols];
  filter(local_sigma, sigma, kernel);
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      im_out[row][col] = sqrt(abs((float) (local_mean[row][col]* local_mean[row][col]- sigma[row][col])));
    }
  }
}

void MSCN(double[][] im_out, double[][] mu, double[][] dev, double[][] im_in){
  int rows = im_in.length;
  int cols = im_in[0].length;
  double C = Double.MIN_VALUE;
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      im_out[row][col] = (im_in[row][col]- mu[row][col])/ (dev[row][col]+ C);
    }
  }
}
