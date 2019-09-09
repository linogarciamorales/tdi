
PImage img_in, img1;
int M, N;
double[][] im1, im2, im3;

double[][] m = {{0, 1, 0},
                {1, 1, 1},
                {0, 1, 0}};

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
  img1 = createImage(N, M, GRAY);
  im1 = new double[N][M];
  convert2double(im1, img_in);

  im2 = new double[N][M];
  float alpha = 0.5;
  brightnessThr(im2, im1, alpha);
  cmp(im2);
  im3 = imcopy(im2);
  //test();
} 

void draw() {
  convert2Pimage(img1, im2);
  image(img1, 0, 0, N, M);
  PImage img2 = createImage(N, M, GRAY);
  convert2Pimage(img2, im3);
  image(img2, N+2, 0, N, M);
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_MorphRegionFilling.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

void mousePressed() {
  int[] pos = {mouseX, mouseY};
  im3 = regionfilling(im3, m, pos, 20);
}

double[][] A1 = {{0, 0, 0, 0, 0, 0, 0, 0},
                 {0, 0, 0, 0, 0, 0, 0, 0},
                 {0, 0, 1, 1, 1, 1, 1, 0},
                 {0, 0, 1, 0, 0, 0, 1, 0},
                 {0, 0, 1, 0, 0, 1, 1, 0},
                 {0, 0, 0, 1, 0, 0, 1, 0},
                 {0, 0, 0, 0, 1, 1, 1, 0},
                 {0, 0, 0, 0, 0, 0, 0, 0}};
                 
double[][] a1 = {{0, 1, 0},
                 {1, 1, 1},
                 {0, 1, 0}};

void test() {
  double[][] A = A1;
  double[][] a = a1;
  printarraybin(A);
  int[] pos = {3, 3};
  double[][] q = regionfilling(A, a, pos, 5);
  printarraybin(q);
}
