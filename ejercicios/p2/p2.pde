int N, M, N2, M2;
PImage b1, b2, b3, b4, b5, b6, b7, b8;
//String filename = "band";
String filename = "B";

void settings() {
  // Make a new instance of a PImage by loading an image file
  b1 = loadImage(filename+"1.png");
  b2 = loadImage(filename+"2.png");
  b3 = loadImage(filename+"3.png");
  b4 = loadImage(filename+"4.png");
  b5 = loadImage(filename+"5.png");
  b6 = loadImage(filename+"6.png");
  b7 = loadImage(filename+"7.png");
  b8 = loadImage(filename+"8.png");

  N = b1.width; 
  M = b1.height;
  N2 = N/2; 
  M2 = M/2;
  size(N, M);
  //size(2*N, 2*M);
}

void setup() {
  smooth();
  noLoop();
  updatePixels();
}

void draw() {
  //e1();
  e2();
  //e3();
  //e4();
  //e5();
  //e6();
  //e7d();
  //e7e();
  //e8();
  //e9();
  //save("RGBimgGamma.png");
}

void e1() {
  double[][] za = zeros(N, M);
  double[][] oa = ones(N, M);
  PImage z = createImage(N, M, GRAY);
  convert2Pimage(z, za);
  PImage o = createImage(N, M, GRAY);
  convert2Pimage(o, oa);
  PImage img1 = setColorImg(o, z, z);
  //img1.save("img1.png");
  PImage img2 = setColorImg(z, o, z);
  PImage img3 = setColorImg(z, z, o);
  PImage img4 = setColorImg(o, o, o);
  image(img1, 0, 0, N2, M2);
  image(img2, N2, 0, N2, M2);
  image(img3, 0, M2, N2, M2);
  image(img4, N2, M2, N2, M2);  
}

void e2() {
  PImage img1 = setColorImg(b3, b2, b1);
  //img1.save("img1.png");
  PImage img2 = setColorImg(b4, b3, b2);
  PImage img3 = setColorImg(b5, b4, b3);
  PImage img4 = setColorImg(b4, b5, b3);
  image(img1, 0, 0, N2, M2);
  image(img2, N2, 0, N2, M2);
  image(img3, 0, M2, N2, M2);
  image(img4, N2, M2, N2, M2);
}

void e3() {
  double alfa = 1.0;
  testSigmoidContrast(alfa);
}

void e4() {
  double alpha = 1.0;
  PImage img1 = setColorImgSigmoidCorrect(b3, b2, b1, alpha);
  PImage img2 = setColorImgSigmoidCorrect(b4, b3, b2, alpha);
  PImage img3 = setColorImgSigmoidCorrect(b5, b4, b3, alpha);
  PImage img4 = setColorImgSigmoidCorrect(b4, b5, b3, alpha);
  image(img1, 0, 0, N2, M2);
  image(img2, N2, 0, N2, M2);
  image(img3, 0, M2, N2, M2);
  image(img4, N2, M2, N2, M2);
}

void e5() {
  double gamma = 1.03;
  testGamma(gamma);
}

void e6() {
  double alpha = 1.0;
  double gamma_red = 1.03;
  double gamma_blue = 0.925;
  PImage img1 = setColorImgSigmoidCorrect(setColorImgGammaCorrect(b3, gamma_red), 
                                          b2, 
                                          setColorImgGammaCorrect(b1, gamma_blue), alpha);
  PImage img2 = setColorImgSigmoidCorrect(setColorImgGammaCorrect(b4, gamma_red), 
                                          b3, 
                                          setColorImgGammaCorrect(b4, gamma_blue), alpha);
  PImage img3 = setColorImgSigmoidCorrect(setColorImgGammaCorrect(b5, gamma_red), 
                                          b4, 
                                          setColorImgGammaCorrect(b3, gamma_blue), alpha);
  PImage img4 = setColorImgSigmoidCorrect(setColorImgGammaCorrect(b4, gamma_red), 
                                          b5, 
                                          setColorImgGammaCorrect(b3, gamma_blue), alpha);
  image(img1, 0, 0, N2, M2);
  //img1.save("img1.png");
  image(img2, N2, 0, N2, M2);
  image(img3, 0, M2, N2, M2);
  image(img4, N2, M2, N2, M2);
}

void e7d() {
  int scl = 2;
  PImage rgb = loadImage("img1.png");
  double[][] h = new double[N][M];
  convert2double(h, rgb, "hue");
  double[][] s = new double[N][M];
  convert2double(s, rgb, "saturation");
  double[][] b = new double[N][M];
  convert2double(b, rgb, "brightness");

  double[][] m = scale(h, scl, scl);
  double[][] H = new double[scl* N][scl* M];
  //transform(H, h, m, bilinear);
  //transform(H, h, m, biquadratic);
  transform(H, h, m, bicubic);
  double[][] S = new double[scl* N][scl* M];
  //transform(S, s, m, bilinear);
  //transform(S, s, m, biquadratic);
  transform(S, s, m, bicubic);
  
  //PImage B8 = createImage(N, M, GRAY);
  //convert2Pimage(B8, b); 
  //matchHisto(b8, B8);
  
  double[][] V = new double[scl*N][scl*M];
  convert2double(V, b8);

  double[][] R = new double[scl*N][scl*M];
  double[][] G = new double[scl*N][scl*M];
  double[][] B = new double[scl*N][scl*M];
  HSV2RGB(R, G, B, H, S, V);
  PImage img = createImage(scl*N, scl*M, RGB);
  convert2Pimage(img, R, G, B);
  //img.save("img3.png");
  image(img, 0, 0, scl*N, scl*M);
}

void e7e() {
  int scl = 2;
  PImage rgb = loadImage("img1.png");
  double[][] h = new double[N][M];
  convert2double(h, rgb, "hue");
  double[][] s = new double[N][M];
  convert2double(s, rgb, "saturation");
  double[][] b = new double[N][M];
  convert2double(b, rgb, "brightness");

  double[][] m = scale(h, scl, scl);
  double[][] H = new double[scl* N][scl* M];
  //transform(H, h, m, bilinear);
  //transform(H, h, m, biquadratic);
  transform(H, h, m, bicubic);
  double[][] S = new double[scl* N][scl* M];
  //transform(S, s, m, bilinear);
  //transform(S, s, m, biquadratic);
  transform(S, s, m, bicubic);
  
  //PImage B8 = createImage(N, M, GRAY);
  //convert2Pimage(B8, b); 
  //matchHisto(b8, B8);
  
  double[][] V = new double[scl*N][scl*M];
  PImage b8eq = histoEQ(b8);
  convert2double(V, b8eq);

  double[][] R = new double[scl*N][scl*M];
  double[][] G = new double[scl*N][scl*M];
  double[][] B = new double[scl*N][scl*M];
  HSV2RGB(R, G, B, H, S, V);
  PImage img = createImage(scl*N, scl*M, RGB);
  convert2Pimage(img, R, G, B);
  //img.save("img3.png");
  image(img, 0, 0, scl*N, scl*M);
}

void e8() {
  double alpha = 1.0;
  double gamma = 1.25;
  PImage img4 = setColorImgSigmoidCorrect(b4, setColorImgGammaCorrect(b3, gamma), 
                                          setColorImgGammaCorrect(b2, gamma), alpha);
  image(img4, 0, 0, N, M);
  //img2.save("b432sig.png");
}

double[] avrc = {5/255.0, 50/255.0, 50/255.0};

double[][] m5 = {{1, 1, 1},
                 {1, 1, 1},
                 {1, 1, 1}};
                 
void e9() {
  double beta = 0.15;
  PImage img_in = loadImage("b432sig.png");
  double[][] R = new double[N][M];
  convert2double(R, img_in, red);

  double[][] G = new double[N][M];
  convert2double(G, img_in, green);

  double[][] B = new double[N][M];
  convert2double(B, img_in, blue);

  double[][] D = new double[N][M];
  D = dist(avrc, R, G, B);
  double[][] S = new double[N][M];
  brightnessThr(S, D, beta);
  //PImage D_out = createImage(N, M, GRAY);
  //convert2Pimage(D_out, D);
  //D_out.save("compmatrix.png");

  PImage img_out = createImage(N, M, RGB);
  convert2Pimage(img_out, S);
  //img_out.save("seg.png");
  image(img_out, 0, 0, N, M);
}
