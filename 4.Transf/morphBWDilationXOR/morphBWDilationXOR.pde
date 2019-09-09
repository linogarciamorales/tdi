
PImage img_in;
int M, N;

double[][] m = {{0, 1, 0},
                {1, 1, 1},
                {0, 1, 0}};

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
  PImage img1 = createImage(N, M, GRAY);
  convert2Pimage(img1, im1);
  image(img1, 0, 0, N, M);

  double[][] im2 = new double[N][M];
  float alpha = 0.5;
  brightnessThr(im2, im1, alpha);
  cmp(im2);

  double[][] im3 = new double[N][M];
  dilation(im3, im2, m);

  double[][] im4 = xor(im2, im3);
  PImage img4 = createImage(N, M, GRAY);
  convert2Pimage(img4, invert(im4));
  image(img4, N+2, 0, N, M);
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_MorphDilationXOR.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}
