PImage img_in;
double[][] im_in;
int N, M;

double Q = 16;
double[][] kernel = {{ 1/Q, 2/Q, 1/Q }, 
                    { 2/Q, 4/Q, 2/Q }, 
                    { 1/Q, 2/Q, 1/Q }};

void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(3*N+4, M);
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

  double mu_in = mean(im_in);
  double std_in = sqrt((float)(var(im_in, mu_in)));
  println(std_in);
  float SNR = 30;
  double g = std_in* pow(10,-SNR/10);

  double[][] im_noised = new double[N][M];
  //addPoissonNoise(im_noised, im_in, g);
  salypimienta(im_noised, im_in, 0.15, 0.45);

  PImage img_noised = createImage(im_noised.length, im_noised[0].length, GRAY);
  convert2Pimage(img_noised, im_noised);
  image(img_noised, N+2, 0, N, M);

  double[][] im_out = new double[N][M];
  medianFilter(im_out, im_noised, 3);
  PImage img_out = createImage(N, M, GRAY);
  convert2Pimage(img_out, im_out);
  image(img_out, 2*N+4, 0, N, M);
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    //save("BW_PoissonNoise.png");
    save("BW_SalYPimientaNoise.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}
void addPoissonNoise(float[][] im_o, float[][] im_i, float gain) {
  for (int x = 0; x < im_o.length; x++) {
    for (int y = 0; y < im_o[0].length; y++) {
      // Creates additive poisson noise
      im_o[x][y] = im_i[x][y]+ gain* poissonNoiseGen(im_i[x][y]);
    }
  }
}
   
int poissonNoiseGen(float pixVal) {
  double L = Math.exp(-(pixVal));
  int k = 0;
  float p = 1;
  do {
    k++;
    // Generate uniform random number u in [0,1] and let p ← p × u.
    p *= randomGaussian();
  } while (p >= L);
  return (k - 1);
}

void salypimienta(double[][] im_o, double[][] im_i, double a, double b) {
  // a + b < 1
  for (int x = 0; x < im_o.length; x++) {
    for (int y = 0; y < im_o[0].length; y++) {
      // Creates additive poisson noise
      double R = 0.5;
      double X = randomGaussian();
      if (X <= a) {
        R = 0;
      }
      if(X > a && X < (a+ b)) {
        R = 1;
      }
      im_o[x][y] = im_i[x][y]+ R;
    }
  }
}

void medianFilter(double[][] im_out, double[][] im_in, int p) {
  // Loop through every pixel in the image
  int q = floor(p/2.0);
  double[][] im_aux = new double[im_in.length+2*(p-1)][im_in[0].length+2*(p-1)];
  imcopy(im_aux, p-1, p-1, im_in, 0, 0, im_in.length, im_in[0].length);
  double[][] im_aux_out = new double[im_in.length+(p-1)][im_in[0].length+(p-1)];
  
  for (int x = q; x < im_aux.length-q; x++) {  // Skip left and right edges
    for (int y = q; y < im_aux[0].length-q; y++) {   // Skip top and bottom edges
      // get win
      float[][] win = new float[p][p];
      for (int u = -q; u <= q; u++) {
        for (int v = -q; v <= q; v++) {
          // println(x+u, y+v);
          // Calculate the adjacent pixel for this kernel point
          win[u+q][v+q] = (float)im_aux[x+u][y+v];
        }
      }
      float[] win_sort = new float[p*p];
      vectorize(win_sort, win);
      win_sort = sort(win_sort);

      im_aux_out[x-q][y-q] = win_sort[(int)ceil(p/2.0)];
    }
  }
  imcopy(im_out, im_aux_out, 0, 0, im_out.length, im_out[0].length);
}

void vectorize(float[] v, float[][] m) {
  int p = m.length;
  for (int x = 0; x < p; x++) {
    for (int y = 0; y < p; y++) {
      v[y*p + x] = m[x][y];   
    }
  }
}
