
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
  //test();
  code();
}

void code() {
  double[][] im1 = new double[N][M];
  convert2double(im1, img_in);

  double[][] im2 = new double[N][M];
  float alpha = 0.5;
  brightnessThr(im2, im1, alpha);
  //cmp(im2);

  PImage img1 = createImage(N, M, GRAY);
  convert2Pimage(img1, im2);
  image(img1, 0, 0, N, M);

  double[][] im3 = new double[N][M];
  hitandmiss(im3, im2, m4);
  //cmp(im3);

  PImage img3 = createImage(N, M, GRAY);
  convert2Pimage(img3, im3);
  image(img3, N+2, 0, N, M);
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_MorphHit_Miss.png");
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
                 {1, 0, 0},
                 {0, 0, 0}};


void test() {
  double[][] q = zeros(8, 8);
  double[][] A = A2;
  double[][] a = a2;
  hitandmiss(q, A, a);
  printarraybin(q);
}
