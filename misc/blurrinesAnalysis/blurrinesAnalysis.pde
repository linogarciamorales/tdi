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
  //convert2Pimage(img_out, im_gr);
  //image(img_out, N+2, 0, N, M);
  
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
 
  Q = 3;
  double[][] kernel4 = {{     0, -1.0/Q, 0 }, 
                        { -1.0/Q, 1.0/Q, 0 }, 
                        {      0,     0, 0 }};
                        
  double[][] N1 = new double[N][M];
  distance(N1, im_in, kernel4, 1);
  double[][] N2 = new double[N][M];
  distance(N2, im_in, kernel4, 2);
  double[][] R = new double[N][M];
  bluriness(R, N2, N1);
  convert2Pimage(img_out, R);
  image(img_out, N+2, 0, N, M);
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

double[][] imcopy(double[][] in) {
  int rows = in.length;
  int cols = in[0].length;
  double[][] out = new double[rows][cols];
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      out[row][col] = in[row][col];
    }
  }
  return out;
}

void imcopy(double[][] out, double[][] in) {
  int rows = out.length;
  int cols = out[0].length;
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      out[row][col] = in[row][col];
    }
  }
}

void imcopy(double[][] out, double[][] in, int start_row, int start_col, int end_row, int end_col) {
  for(int row = 0; row < end_row; row++) {
    for(int col = 0; col < end_col; col++) {
      out[row][col] = in[row+start_row][col+start_col];
    }
  }
}

void imcopy(double[][] out, int start_row_out, int start_col_out, 
            double[][] in, int start_row_in, int start_col_in, int end_row, int end_col) {
  for(int row = 0; row < end_row; row++) {
    for(int col = 0; col < end_col; col++) {
      out[row+start_row_out][col+start_col_out] = in[row+start_row_in][col+start_col_in];
    }
  }
}

void distance(double[][] im_out, double[][] im_in, double[][] kernel, int L) {
  // Loop through every pixel in the image
  int p = kernel.length;
  int q = floor(p/2.0);
  double[][] im_aux = new double[im_in.length+2*(p-1)][im_in[0].length+2*(p-1)];
  double[][] im_aux_out = new double[im_in.length+(p-1)][im_in[0].length+(p-1)];
  
  imcopy(im_aux, p-1, p-1, im_in, 0, 0, im_in.length, im_in[0].length);

  for (int row = q; row < im_aux.length-q; row++) {  // Skip left and right edges
    for (int col = q; col < im_aux[0].length-q; col++) {   // Skip top and bottom edges

      double sum = 0; // Kernel sum for this pixel
      for (int u = -q; u <= q; u++) {
        for (int v = -q; v <= q; v++) {
          // println(x+u, y+v);
          // Calculate the adjacent pixel for this kernel point
          double val = im_aux[row+u][col+v];
          // Multiply adjacent pixels based on the kernel values
          if (L == 1) {
            sum += (kernel[q-u][q-v] + val);
          }
          if (L == 2) {
            sum += (kernel[q-u][q-v] + val)* (kernel[q-u][q-v] + val);
          }
        }
      }
      // For this pixel in the new image, set the gray value
      // based on the sum from the kernel
      im_aux_out[row-q][col-q] = sum;
    }
  }
  imcopy(im_out, im_aux_out, 0, 0, im_out.length, im_out[0].length);
}

/*
 RATIO OF L1-L2 DISTANCE
Let p(x,y) denotes your current pixel.
calculate the average L1-distance of adjacent pixels:
N1=1/(2*N_pixel) * sum( abs(p(x,y)-p(x-1,y)) + abs(p(x,y)-p(x,y-1)) )
then the average L2 distance:
N2= 1/(2*N_pixel) * sum( (p(x,y)-p(x-1,y))^2 + (p(x,y)-p(x,y-1))^2 )
then take ratio R which is measure of blurriness
R = N2 / (N1*N1)
This is for gray scale images, for color you do this for each channel separately.
*/
void bluriness(double[][] R,double[][] N2,double[][] N1) {
  int rows = R.length;
  int cols = R[0].length;
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      R[row][col] = N2[row][col]/ (N1[row][col]* N1[row][col]);
    }
  }
}
