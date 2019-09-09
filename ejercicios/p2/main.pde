PImage setColorImg(PImage chn1, PImage chn2, PImage chn3) {
  double[][] r = new double[N][M];
  convert2double(r, chn1);
  double[][] g = new double[N][M];
  convert2double(g, chn2);
  double[][] b = new double[N][M];
  convert2double(b, chn3);

  PImage img = createImage(N, M, RGB);
  convert2Pimage(img, r, g, b);
  return img;
}

void testGamma(double k) {
  println("0.25, "+gamma(0.25, k));
  println("0.5, "+gamma(0.5, k));
  println("0.75, "+gamma(0.75, k));
  println("1.0, "+gamma(1.0, k));
  println("1.5, "+gamma(1.5, k));
}

double gamma(double x, double g) {
  return Math.pow(x, g);
}

double[][] gammaCorrect(double[][] imi, double g) {
  int rows = imi.length;
  int cols = imi[0].length;
  double[][] imo = new double[rows][cols];
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
     double val = imi[row][col];
      imo[row][col] = gamma(val, g); 
    }
  }
  return imo;
}

PImage setColorImgGammaCorrect(PImage chn, double k) {
  double[][] imi = new double[N][M];
  convert2double(imi, chn);

  double[][] imo = gammaCorrect(imi, k);
  PImage img = createImage(N, M, GRAY);
  convert2Pimage(img, imo);
  return img;
}

void testSigmoidContrast(double k) {
  println("0.25, "+sigmoidContrast(0.25, k));
  println("0.5, "+sigmoidContrast(0.5, k));
  println("0.75, "+sigmoidContrast(0.75, k));
  println("1.0, "+sigmoidContrast(1.0, k));
  println("1.5, "+sigmoidContrast(1.5, k));
}

double sigmoid(double x) {
  return 1/(1+ Math.exp(-x));
}

double sigmoidContrast(double x, double k) {
  return x+ x* k* sigmoid(x); 
}

double[][] sigmoidContrast(double[][] imi, double k) {
  int rows = imi.length;
  int cols = imi[0].length;
  double[][] imo = new double[rows][cols];
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
     double val = imi[row][col];
      imo[row][col] = sigmoidContrast(val, k); 
    }
  }
  return imo;
}

PImage setColorImgSigmoidCorrect(PImage chn1, PImage chn2, PImage chn3, double k) {
  double[][] r = new double[N][M];
  convert2double(r, chn1);
  double[][] g = new double[N][M];
  convert2double(g, chn2);
  double[][] b = new double[N][M];
  convert2double(b, chn3);

  double[][] re = sigmoidContrast(r, k);
  double[][] ge = sigmoidContrast(g, k);
  double[][] be = sigmoidContrast(b, k);
  PImage img = createImage(N, M, RGB);
  convert2Pimage(img, re, ge, be);
  return img;
}
