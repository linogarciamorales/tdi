
PImage img_in;
int M, N;

double[][] m1 = {{0, 1, 0},
                 {1, 1, 1},
                 {0, 1, 0}};

double[][] m2 = {{0, 0, 1, 0, 0},
                 {0, 1, 1, 1, 0},
                 {1, 1, 1, 1, 1},
                 {0, 1, 1, 1, 0},
                 {0, 0, 1, 0, 0}};

double[][] m3 = {{0, 0, 0, 1, 0, 0, 0},
                 {0, 0, 1, 1, 1, 0, 0},
                 {0, 1, 1, 1, 1, 1, 0},
                 {1, 1, 1, 1, 1, 1, 1},
                 {0, 1, 1, 1, 1, 1, 0},
                 {0, 0, 1, 1, 1, 0, 0},
                 {0, 0, 0, 1, 0, 0, 0}};

double[][] m4 = {{1, 1, 0},
                 {1, 0, 0},
                 {0, 0, 0}};

double[][] m5 = {{0, 0, 0},
                 {0, 0, 1},
                 {0, 1, 1}};

void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+2, M);
}

void setup() {
  background(255);
  smooth();
  noLoop();
} 

void draw() {
  double[][] im1 = new double[N][M];
  convert2double(im1, img_in);

  double mu = mean(im1);
  double std = sqrt((float)var(im1, mu));
  println(std);
  float SNR = 5;
  double g = std* pow(10,-SNR/10);

  double[][] im1_noised = new double[N][M];
  gaussianNoiseGen(im1_noised, im1, g);
  
  double[][] im2 = new double[N][M];
  float alpha = 0.5;
  brightnessThr(im2, im1_noised, alpha);
  
  PImage img1 = createImage(N, M, GRAY);
  //convert2Pimage(img1, invert(im2));
  convert2Pimage(img1, im2);
  image(img1, 0, 0, N, M);

  double[][] im3 = new double[N][M];
  opening(im3, invert(im2), m1);
  //opening(im3, im2, m1);
  double[][] im4 = new double[N][M];
  closing(im4, im3, m1);

  PImage img5 = createImage(N, M, GRAY);
  convert2Pimage(img5, invert(im4));
  image(img5, N+2, 0, N, M);
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_morphFilter.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}
