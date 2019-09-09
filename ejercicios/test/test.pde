PImage black, white;
int N = 100, M = 100;

void settings() {
  black = loadImage("black.png");  
  white = loadImage("white.png");
  size(N, M);
}

void setup() {
  smooth();
  noLoop();
}

void draw() {
  test();
}

void create() {
  double[][] b = zeros(N, M);
  PImage black = createImage(N, M, GRAY);
  convert2Pimage(black, b);
  image(black, 0, 0, N, M);
  save("black.png");
  double[][] w = ones(100, 100);
  PImage white = createImage(w.length, w[0].length, GRAY);
  convert2Pimage(white, w);
  image(white, 0, 0, N, M);
  save("white.png");
}

void test() {
  //PImage img = setColorImg(white, black, black);
  //PImage img = setColorImg(white, black, white);
  PImage img = setColorImg(black, white, white);
  image(img, 0, 0, N, M);
}


PImage setColorImg(PImage chn1, PImage chn2, PImage chn3) {
  double[][] r = new double[N][M];
  convert2double(r, chn1);
  double[][] g = new double[N][M];
  convert2double(g, chn2);
  double[][] b = new double[N][M];
  convert2double(b, chn3);

  PImage img_RGB = createImage(N, M, RGB);
  convert2Pimage(img_RGB, r, g, b);
  return img_RGB;
}
