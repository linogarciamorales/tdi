
PImage img_in;
int M, N;

double[][] m = {{0, 0, 0},
                {0, 1, 0},
                {1, 1, 1}};

void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/silueta.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  //img_in = loadImage("../../img/patt2.png");
  N = img_in.width; 
  M = img_in.height;
  size(3*N+ 4, M);
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
    save("BW_MorphSkeletons.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

void code() {
  double[][] im1 = new double[N][M];
  convert2double(im1, img_in);

  //double[][] im2 = imcopy(im1);
  double[][] im2 = new double[N][M];
  float alpha = 0.5;
  brightnessThr(im2, im1, alpha);
  
  PImage img1 = createImage(N, M, GRAY);
  convert2Pimage(img1, im2);
  image(img1, 0, 0, N, M);
  cmp(im2);
 
  double[][] im3 = MAT(im2);

  PImage img3 = createImage(N, M, GRAY);
  normalize(im3);
  convert2Pimage(img3, im3);
  image(img3, N+2, 0, N, M);
  
  double[][] im4 = skel(im2);

  PImage img4 = createImage(N, M, GRAY);
  convert2Pimage(img4, im4);
  image(img4, 2*N+4, 0, N, M);
}

double[][] A1 = {{0, 0, 0, 0, 0, 0, 0, 0},
                 {0, 1, 1, 1, 1, 1, 1, 0},
                 {0, 1, 1, 1, 1, 1, 1, 0},
                 {0, 1, 1, 1, 1, 1, 1, 0},
                 {0, 1, 1, 1, 1, 1, 1, 0},
                 {0, 1, 1, 1, 1, 1, 1, 0},
                 {0, 1, 1, 1, 1, 1, 1, 0},
                 {0, 0, 0, 0, 0, 0, 0, 0}};

void test() {
  double[][] A = A1;
  double[][] q = skel(A);
  printarray(q);
}
