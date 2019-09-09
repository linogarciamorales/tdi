
PImage img_in, img_mask;
int N, M;

void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW_128_128.png");
  img_mask = loadImage("../../img/LP.png");
  N = img_in.width; 
  M = img_in.height;
  size(4*N+6, M);
}

void setup() {
  background(255);
  smooth();
  noLoop();
} 

void draw() {
  //test();

  image(img_in, 0, 0, N, M);

  //double[][] mask_re = ones(N, M);
  double[][] mask_re = new double[N][M];
  convert2double(mask_re, img_mask);
  double[][] mask_im = zeros(N, M);

  PImage img_mask = createImage(N, M, GRAY);
  convert2Pimage(img_mask, mask_re);
  image(img_mask, 3*N+6, 0, N, M);

  double[][] im_re = new double[N][M];
  convert2double(im_re, img_in);
  double[][] im_im = zeros(N, M);

  double[][] im_noised = new double[N][M];
  //imcopy(im_noised, im_re);

  double mu_in = mean(im_re);
  double std_in = sqrt((float)var(im_re, mu_in));
  float SNR = 0;
  double g = std_in* pow(10,-SNR/10);

  gaussianNoiseGen(im_noised, im_re, g);
  PImage img_noised = createImage(N, M, GRAY);
  convert2Pimage(img_noised, im_noised);
  image(img_noised, N+2, 0, N, M);

  fft2(im_noised, im_im);
  fft2(mask_re, mask_im);

  multiply(im_noised, im_im, mask_re, mask_im); 
  
  ifft2(im_noised, im_im);
  normalize(im_noised);

  //double[][] im_abs = new double[N][M];
  //abs(im_abs, mask_re, mask_im); 
  double[][] im_scr = new double[N][M];
  im_scr = shift(im_noised);

  PImage img_rec = createImage(N, M, GRAY);
  convert2Pimage(img_rec, im_scr);
  image(img_rec, 2*N+4, 0, N, M);

}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_DFT_LP.png");
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
  printarray(a);
  println();
  printarray(a, 1, 2, 3, 3);
*/
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

  printarray(re);
  println();
  fft2(re, im);
  printarray(re);
  println();
  ifft2(re, im);
  
  double[][] rec_a = new double[rows][cols];
  imcopy(rec_a, re, 0, 0, rows, cols);
  printarray(rec_a);

}
