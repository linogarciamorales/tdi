PImage img_in;

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
  size(2*N+ 2, M);
}

void setup() {
  background(255);
  smooth();
  noLoop();
} 

void draw() {
  image(img_in, 0, 0, N, M);
  
  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in);
  
  double[][] im_gr = new double[N][M];
  filter(im_gr, im_in, laplace);
  
  PImage img_out = createImage(N, M, GRAY);
  convert2Pimage(img_out, im_gr);
  image(img_out, N+2, 0, N, M);
  
  double v = variance_of_laplacian(im_in);
  println("blurry: "+ v);

  double Q = 9.0;
  double[][] kernel1 = {{ 1.0/Q, 1.0/Q, 1.0/Q }, 
                        {  1.0/Q, 1.0/Q, 1.0/Q }, 
                        {  1.0/Q, 1.0/Q, 1.0/Q }};
  Q = 25;                      
  double[][] kernel2 = {{ 1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q }, 
                        {  1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q }, 
                        {  1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q },
                        {  1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q },
                        {  1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q }};
                      
  Q = 49;            
  double[][] kernel3 = {{ 1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q}, 
                        {  1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q}, 
                        {  1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q},
                        {  1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q},
                        {  1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q},
                        {  1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q},
                        {  1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q}};
                    
  double[][] im = new double[N][M];
  filter(im, im_in, kernel1);
  v = variance_of_laplacian(im);
  println("Q: 3, blurry: "+ v);

  filter(im, im_in, kernel2);
  v = variance_of_laplacian(im);
  println("Q: 5, blurry: "+ v);

  filter(im, im_in, kernel3);
  v = variance_of_laplacian(im);
  println("Q: 7, blurry: "+ v);
}

void keyPressed() {
  if (key == 'q' || key == 'Q') {
    exit();    
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

double mean(double[][] a) {
  int rows = a.length;
  int cols = a[0].length;
  double mu = 0;
  for(int row=0; row < rows; row++) {
    for(int col=0; col < cols; col++) {
      mu += a[row][col];
    }
  }
  mu /= (rows* cols);
  return mu;
}

double var(double[][] a) {
  int rows = a.length;
  int cols = a[0].length;
  double mu = mean(a); 
  double d2 = 0;
  for(int row=0; row < rows; row++) {
    for(int col=0; col < cols; col++) {
      double xy = a[row][col]- mu;
      d2 += xy* xy;
    }
  }
  d2 /= (rows* cols);
  return d2;
} 

int[] size(double[][] im) {
  int[] p = new int[2];
  p[0] = im.length;
  p[1] = im[0].length;
  return p;
}

double[][] zeros(int[] p) {
  double[][] out = new double[p[0]][p[1]];
  for(int row = 0; row < p[0]; row++) {
    for(int col = 0; col < p[1]; col++) {
      out[row][col] = 0;
    }
  }
  return out;
}

double variance_of_laplacian(double[][] im) {
  double[][] im_gr = zeros(size(im));
  filter(im_gr, im, laplace);

  double v = var(im_gr);
  return v;
}
