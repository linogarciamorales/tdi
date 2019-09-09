
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
                 {0, 1, 0},
                 {1, 0, 0}};

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
  //test();
  code();
}

void code() {
  double[][] im1 = new double[N][M];
  convert2double(im1, img_in);

  double[][] im2 = new double[N][M];
  float alpha = 0.5;
  brightnessThr(im2, im1, alpha);
  cmp(im2);

  PImage img1 = createImage(N, M, GRAY);
  convert2Pimage(img1, im2);
  image(img1, 0, 0, N, M);

  double[][] m = m4;
  double[][] im3 = new double[N][M];
  opening(im3, im2, m);

  PImage img3 = createImage(N, M, GRAY);
  convert2Pimage(img3, im3);
  image(img3, N+2, 0, N, M);

  double[][] im4 = new double[N][M];
  closing(im4, im2, m);

  PImage img4 = createImage(N, M, GRAY);
  convert2Pimage(img4, im4);
  image(img4, 2*N+4, 0, N, M);
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_MorphOpeningClosing.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

double[][] A1 = {{0, 0, 0, 0, 0, 0, 0, 0},
                {1, 1, 1, 1, 1, 1, 1, 0},
                {0, 0, 0, 1, 1, 1, 1, 0},
                {0, 0, 0, 1, 1, 1, 1, 0},
                {0, 0, 1, 1, 1, 1, 1, 0},
                {0, 0, 0, 1, 1, 1, 1, 0},
                {0, 0, 1, 1, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0}};

double[][] a1 = {{1, 1, 1},
                {1, 1, 1},
                {1, 1, 1}};

double[][] A2 = {{0, 0, 0, 0, 0, 0, 0, 0},
                 {0, 0, 0, 0, 1, 0, 0, 0},
                 {0, 0, 1, 0, 1, 1, 0, 0},
                 {0, 1, 1, 1, 1, 1, 1, 0},
                 {0, 0, 1, 1, 1, 1, 1, 0},
                 {0, 0, 0, 1, 1, 1, 1, 0},
                 {0, 0, 0, 0, 1, 1, 0, 0},
                 {0, 0, 0, 0, 0, 1, 0, 0}};

double[][] a2 = {{1, 1, 0},
                 {0, 1, 0},
                 {1, 0, 0}};


void test() {
  double[][] q = zeros(8, 8);
  double[][] A = A2;
  double[][] a = a2;
  closing(q, A, a);
  printarraybin(q);
  opening(q, A, a);
  printarraybin(q);
}
