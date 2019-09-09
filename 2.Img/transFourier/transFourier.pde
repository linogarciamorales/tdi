
PImage img_in;
double[][] im_in;
int N, M;

void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW_128_128.png");
  N = img_in.width; 
  M = img_in.height;
  size(3*N+4, M);
  println(N, M);
}

void setup() {
  background(255);
  smooth();
  noLoop();
} 

void draw() {
  //test();
  image(img_in, 0, 0, N, M);

  double[][] im_re = new double[N][M];
  convert2double(im_re, img_in);
  int rows = im_re.length;
  int cols = im_re[0].length;
  double[][] im_im = zeros(rows, cols);

  fft2(im_re, im_im);
  double[][] im_abs = new double[N][M];
  abs(im_abs, im_re, im_im); 
  //normalize(im_abs_full); 
  PImage img_abs = createImage(im_abs.length, im_abs[0].length, GRAY);
  convert2Pimage(img_abs, im_abs);
  image(img_abs, N+2, 0, N, M);
  
  ifft2(im_re, im_im);
  PImage img_rec = createImage(im_re.length, im_re[0].length, GRAY);
  convert2Pimage(img_rec, im_re);
  image(img_rec, 2*N+4, 0, N, M);

}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_DFT.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

double[][] a = {{1, 2, 3, 4, 5},
               {4, 5, 6, 5, 8},
               {7, 8, 9, 10, 11},
               {12, 13, 14, 15, 16},
               {17, 18, 19, 20, 21}};

double[][] x = {{1, 2, 3, 4, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0}};

double[] y = {1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

double[][] z = {{12, 13, 14, 15, 16, 15, 10, 5},
               { 0,  0,  0,  0,  0,  0,  0, 0}};

void test() {
/*
  // FFT 1D test
  int n = x[0].length;

  double[] re = new double[n];
  arraycopy(re, y, 0, n);
  double[] im = new double[n];
  arraycopy(im, y, n, 2*n);
  fft(re, im, true);
  printarray(re);
  printarray(im);
*/
  // FFT 2D test

  int rows = a.length;
  int cols = a[0].length;
  double[][] b = zeros(rows, cols);
  double[][] re = reflect(a);
  double[][] im = reflect(b);
  int arows = re.length;
  int acols = re[0].length;

  fft2(re, im);
  ifft2(re, im);
  
  printarray(a);
  println();
  printarray(re);
  println();
  double[][] rec_a = new double[rows][cols];
  imcopy(rec_a, re, 0, 0, rows, cols);
  printarray(rec_a);

}
