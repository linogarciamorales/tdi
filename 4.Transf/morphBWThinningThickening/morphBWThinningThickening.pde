
PImage img_in;
int M, N;
/*
double[][] m = {{0, 0, 0},
                {0, 1, 0},
                {1, 1, 1}};
*/
double[][] m = {{0, 0, 0},
                {1, 1, 0},
                {0, 0, 0}};

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

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_MorphThinningThickening.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
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

  double[][] im3 = new double[N][M];
  thinning(im3, im2, m);

  PImage img3 = createImage(N, M, GRAY);
  convert2Pimage(img3, im3);
  image(img3, N+2, 0, N, M);

  double[][] im4 = new double[N][M];
  thickening(im4, im2, m);
  
  PImage img4 = createImage(N, M, GRAY);
  convert2Pimage(img4, im4);
  image(img4, 2*N+4, 0, N, M);
}

double[][] A1 = {{0, 0, 1, 1, 1, 0, 0, 0},
                 {0, 0, 0, 1, 1, 1, 0, 0},
                 {0, 0, 0, 0, 1, 1, 1, 0},
                 {0, 0, 0, 0, 0, 1, 1, 1},
                 {0, 0, 0, 0, 1, 1, 0, 0},
                 {0, 0, 0, 0, 1, 0, 0, 0},
                 {0, 0, 1, 1, 0, 0, 0, 0},
                 {0, 1, 1, 0, 0, 0, 0, 0}};

double[][] a1 = {{1, 0, 0},
                 {0, 0, 1},
                 {0, 0, 0}};


void test() {
  double[][] q = zeros(8, 8);
  double[][] A = A1;
  double[][] a = a1;
  thinning(q, A, a);
  printarraybin(q);
}
