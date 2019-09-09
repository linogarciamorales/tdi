
PImage img_in;
int M, N;

double[][] m = {{1, 1, 1},
                {1, 1, 1},
                {1, 1, 1}};

/*
double[][] m = {{0, 0, 1, 0, 0},
                {0, 1, 1, 1, 0},
                {1, 1, 1, 1, 1},
                {0, 1, 1, 1, 0},
                {0, 0, 1, 0, 0}};

double[][] m = {{0, 1, 0},
                {0, 1, 1},
                {1, 1, 0}};

double[][] m = {{1, 1, 0},
                {0, 1, 0},
                {1, 0, 0}};
*/

void settings() { 
  // Make a new instance of a PImage by loading an image file
  //img_in = loadImage("../../img/detalle_BW.png");
  //img_in.resize(img_in.width/ 2, img_in.height/ 2);
  img_in = loadImage("../../img/patt4.png");
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

  PImage img1 = createImage(N, M, GRAY);
  convert2Pimage(img1, im1);
  image(img1, 0, 0, N, M);

  double[][] im2 = new double[N][M];
  float alpha = 0.5;
  //brightnessThrInv(im2, im1, alpha);
  brightnessThr(im2, im1, alpha);

  PImage img2 = createImage(N, M, GRAY);
  convert2Pimage(img2, im2);
  image(img2, N+2, 0, N, M);

  double[][] im3 = new double[N][M];
  erosion(im3, im2, m);
  //erosion(im3, im3, m);
  //erosion(im3, im3, m);

  PImage img3 = createImage(N, M, GRAY);
  convert2Pimage(img3, im3);
  image(img3, 2*N+4, 0, N, M);
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_MorphErotion.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

double[][] A = {{0, 0, 0, 0, 0, 0, 0, 0},
                {1, 1, 1, 1, 1, 1, 1, 0},
                {0, 0, 0, 1, 1, 1, 1, 0},
                {0, 0, 0, 1, 1, 1, 1, 0},
                {0, 0, 1, 1, 1, 1, 1, 0},
                {0, 0, 0, 1, 1, 1, 1, 0},
                {0, 0, 1, 1, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 1, 0, 0}};

double[][] a = {{1, 1, 1},
                {1, 1, 1},
                {1, 1, 1}};

double[][] b = {{1, 1, 1, 1, 1},
                {1, 1, 1, 1, 1},
                {1, 1, 1, 1, 1},
                {1, 1, 1, 1, 1},
                {1, 1, 1, 1, 1}};
                
void test() {
  double[][] q = zeros(8, 8);
  erosion(q, A, b);
  printarray(q);
  //printarray(reflect(invert(a)));

}
