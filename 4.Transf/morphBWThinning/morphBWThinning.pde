
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
  img_in = loadImage("../../img/patt5.png");
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

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_MorphThinning.png");
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
  
  PImage img1 = createImage(N, M, GRAY);
  convert2Pimage(img1, im2);
  image(img1, 0, 0, N, M);

  double[][] im3 = new double[N][M];
  thinning(im3, im2, m);

  PImage img3 = createImage(N, M, GRAY);
  convert2Pimage(img3, im3);
  image(img3, N+2, 0, N, M);
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
