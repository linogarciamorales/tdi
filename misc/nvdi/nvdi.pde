int N, M;
PImage b1, b2, b3, b4, b5, b6, b7, b8;
//String filename = "band";
String filename = "../../img/data/B";

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
//  N = b1.width; 
//  M = b1.height;

  N = b1.width; 
  M = b1.height;
  //size(N, M);
  size(2*N, M);
  //size(2*N, 2*M);
}

void setup() {
  smooth();
  noLoop();
  updatePixels();
}

void draw() {
  //code1();
  //code2();
  //code3();
  //code4();
  //code5();
  //code6();
  //save("RGBimgGamma.png");
  code7();
}

void code1() {
  PImage img1 = setColorImg(b3, b2, b1);
  //img1.save("img1.png");
  PImage img2 = setColorImg(b4, b3, b2);
  PImage img3 = setColorImg(b5, b4, b3);
  PImage img4 = setColorImg(b4, b5, b3);
  image(img1, 0, 0, N, M);
  image(img2, N, 0, N, M);
  image(img3, 0, M, N, M);
  image(img4, N, M, N, M);
}

void code2() {
  double k = 1.0;
  //testSigmoidContrast(k);

  PImage img1 = setColorImgSigmoidCorrect(b3, b2, b1, k);
  PImage img2 = setColorImgSigmoidCorrect(b4, b3, b2, k);
  PImage img3 = setColorImgSigmoidCorrect(b5, b4, b3, k);
  PImage img4 = setColorImgSigmoidCorrect(b4, b5, b3, k);
  image(img1, 0, 0, N, M);
  image(img2, N, 0, N, M);
  image(img3, 0, M, N, M);
  image(img4, N, M, N, M);
}

void code3() {
  double k = 1.03;
  //testGamma(k);
  PImage img1 = setColorImgSigmoidCorrect(setColorImgGammaCorrect(b3, 1.03), 
                                          b2, 
                                          setColorImgGammaCorrect(b1, 0.925), k);
  PImage img2 = setColorImgSigmoidCorrect(setColorImgGammaCorrect(b4, 1.03), 
                                          b3, 
                                          setColorImgGammaCorrect(b4, 0.925), k);
  PImage img3 = setColorImgSigmoidCorrect(setColorImgGammaCorrect(b5, 1.03), 
                                          b4, 
                                          setColorImgGammaCorrect(b3, 0.925), k);
  PImage img4 = setColorImgSigmoidCorrect(setColorImgGammaCorrect(b4, 1.03), 
                                          b5, 
                                          setColorImgGammaCorrect(b3, 0.925), k);
  image(img1, 0, 0, N, M);
  //img1.save("img1.png");
  image(img2, N, 0, N, M);
  image(img3, 0, M, N, M);
  image(img4, N, M, N, M);
}

void code4() {
  double k = 1.03;
  PImage img4 = setColorImgSigmoidCorrect(b4, setColorImgGammaCorrect(b5, 1.25), 
                                          setColorImgGammaCorrect(b3, 1.25), k);
  PImage img2 = setColorImgSigmoidCorrect(b4, setColorImgGammaCorrect(b3, 1.25), 
                                          setColorImgGammaCorrect(b2, 1.25), k);
  image(img4, 0, 0, N, M);
  image(img2, N, 0, N, M);
  //img2.save("b432sig.png");
}

void code5() {
  int scl = 2;
  PImage rgb = loadImage("img1.png");
  double[][] h = new double[N][M];
  convert2double(h, rgb, "hue");
  double[][] m = scale(h, scl, scl);
  double[][] H = new double[scl* N][scl* M];
  //transform(H, h, m, bilinear);
  //transform(H, h, m, biquadratic);
  transform(H, h, m, bicubic);
  double[][] s = new double[N][M];
  convert2double(s, rgb, "saturation");
  double[][] S = new double[scl* N][scl* M];
  //transform(S, s, m, bilinear);
  //transform(S, s, m, biquadratic);
  transform(S, s, m, bicubic);
  
  double[][] b = new double[N][M];
  convert2double(b, rgb, "brightness");

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

double[] avrc = {5/255.0, 50/255.0, 50/255.0};

double[][] m5 = {{1, 1, 1},
                 {1, 1, 1},
                 {1, 1, 1}};
                 
void code6() {
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
  brightnessThr(S, D, 0.15);
  //PImage D_out = createImage(N, M, GRAY);
  //convert2Pimage(D_out, D);
  //D_out.save("compmatrix.png");

  PImage img_out = createImage(N, M, RGB);
  convert2Pimage(img_out, S);
  //img_out.save("seg.png");
  image(img_out, 0, 0, N, M);
}

void code7() {
  PImage img1 = setColorImg(b4, b3, b2);
  image(img1, 0, 0, N, M);
  PImage img2 = sum(b4, b3);        // NIR+R
  PImage img3 = sub(b4, b3);        // NIR-R
  PImage img4 = div(img3, img2);    // NDVI = (NIR-R)/(NIR+R)
  PImage img5 = brightnessThr(img4, 0.25);  // NDVI > 0, vegetacion
  image(img5, N, 0, N, M);
}
