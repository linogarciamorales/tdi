import peasy.*;
boolean record;
PeasyCam cam;

int N, M;
double[][] kernel1, kernel2, kernel3;
double[][] laplace = {{ 0,  1, 0 }, 
                     { 1, -4, 1 }, 
                     { 0,  1, 0 }};
double[][] gradx = {{ 0, -1.0, 0 }, 
                   { 0, 1.0, 0 }, 
                   { 0, 0, 0 }};

int meshSize = 5;
int resX, resY;
double[][] val; 
double sigma = 1.0;
int W = 55;

void settings() {
  // Make a new instance of a PImage by loading an image file
  N = W; 
  M = W;
  size(512, 512, P3D);
  resX = W;
  resY = W;
  val = new double[resX][resY];
  kernel1 = new double[W][W];
  kgauss(kernel1, sigma);
  //printarray(kernel1);

  kernel2 = new double[W][W];
  gaussianKernel(kernel2, sigma);
  //printarray(kernel2);
  
  kernel3 = new double[W][W];
  LoG(kernel3, sigma);
  
  //filter(kernel2, kernel2, gradx);  // first derivative
  //filter(kernel2, kernel2, laplace);  // second derivative
}

void setup() {
  smooth();
  background(255);
  cam = new PeasyCam(this, 128);
}

void draw() {
  background(255);
  //show(kernel1, 64);
  //show(kernel2, 2048);
  show(kernel3, 64);
}

void kgauss(double[][] kernel, double sigma) {
  int q = floor(kernel.length/2.0);
  for (int v=-q, y=0; v <= q; v++, y=v+q) {
    for (int u=-q, x=0; u <= q; u++, x=u+q) {
      double A = 1.0/(2.0* PI* sigma*sigma);
      //double A = 1.0;
      kernel[x][y] = A* exp((float)(-(u*u+v*v)/(2*sigma*sigma)));
    }
  }
}

void gaussianKernel(double[][] kernel, double sigma) {
  int W = kernel.length;
  int q = floor(kernel.length/2.0);
  double sum = 0;
  for (int v=-q, y=0; v <= q; v++, y=v+q) {
    for (int u=-q, x=0; u <= q; u++, x=u+q) {
      //double A = 1.0/(2.0* PI* sigma*sigma);
      double A = 1.0;
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

void filter(double[][] im_out, double[][] im_in, double[][] win) {
  // Loop through every pixel in the image
  int q = floor(win.length/2.0);

  for (int x = q; x < N-q; x++) {  // Skip left and right edges
    for (int y = q; y < M-q; y++) {   // Skip top and bottom edges

      double sum = 0; // Kernel sum for this pixel
      for (int v = -q; v <= q; v++) {
        for (int u = -q; u <= q; u++) {
          // Calculate the adjacent pixel for this kernel point
          double val = im_in[x+u][y+v];
          // Multiply adjacent pixels based on the kernel values
          sum += win[u+q][v+q] * val;
        }
      }
      // For this pixel in the new image, set the gray value
      // based on the sum from the kernel
      im_out[x][y] = sum;
    }
  }
}

void show(double[][] win, double scale) {
  translate(-resX/2*meshSize,-resY/2*meshSize);

  for(int x =0; x<resX; x++){
    for(int y =0; y<resY; y++){
      val[x][y] = scale* win[x][y]*mouseX/width;
    }
  }

  for(int x =0; x < resX-1; x++){
    for(int y =0; y < resY-1; y++){
      beginShape();
      colorMode(HSB, 255);
      fill((float)val[x][y], 255, 255);
      vertex(x*meshSize, y*meshSize, (float)val[x][y] );
      vertex((x+1)*meshSize,y*meshSize,(float)val[x+1][y] );
      vertex((x+1)*meshSize,(y+1)*meshSize,(float)val[x+1][y+1] );
      vertex(x*meshSize,(y+1)*meshSize,(float)val[x][y+1] );
      endShape(CLOSE);
    } 
  }
}

void printarray(double[][] a) {
  for(int x = 0; x < a.length; x++) {
    for(int y = 0; y < a[0].length; y++) {
      print(a[x][y]+ " ");
    }
    println();
  }
  println();
}

void LoG(double[][] kernel, double sigma) {
  int q = floor(kernel.length/2.0);
  double sigma2 = sigma*sigma;
  double sum = 0;
  for (int v=-q, y=0; v <= q; v++, y=v+q) {
    for (int u=-q, x=0; u <= q; u++, x=u+q) {
      double r2 = u*u + v*v;
      double r = sqrt((float)r2);
      kernel[x][y] = -(1.0-r/sigma2)* exp((float)(-(r2)/(2*sigma2)));
      sum += kernel[x][y];
    }
  }

  int W = kernel.length;
  for(int x = 0; x < W; x++) {
    for(int y = 0; y < W; y++) {
      kernel[x][y] /= sum;
    }
  }
}
