
PImage img_in;
int M, N;
int K = 512;

String bilinear = new String("bilinear");
String biquadratic = new String("biquadratic");
String bicubic = new String("bicubic");

void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW_128_128.png");
  N = img_in.width; 
  M = img_in.height;
  size((int)(3.5*N)+K+4, K);
}

void setup() {
  background(255);
  smooth();
  noLoop();
} 

void draw() {
  double[][] im = new double[N][M];
  convert2double(im, img_in);
  double[][] im1 = ones(K, K);
  imcopy(im1, 0, 0, im, 0, 0, N, M);
  PImage img1 = createImage(K, K, GRAY);
  convert2Pimage(img1, im1);
  image(img1, 0, 0, K, K);

  double[][] m;
  m = scale(im, 2.5, 2.5);
  printarray(m);
  double[][] imx2 = new double[(int)(2.5* N)][(int)(2.5* M)];
  transform(imx2, im, m, bilinear);
  //transform(imx2, im, m, biquadratic);
  //transform(imx2, im, m, bicubic);
  double[][] im2 = ones(K, K);
  imcopy(im2, 0, 0, imx2, 0, 0, imx2.length, imx2[0].length);
  PImage img2 = createImage(K, K, GRAY);
  convert2Pimage(img2, im2);
  image(img2, N+2, 0, K, K);

  m = scale(im, 4, 4);
  printarray(m);
  double[][] imx3 = new double[4*N][4*M];
  transform(imx3, im, m, bilinear);
  //transform(imx3, im, m, biquadratic);
  transform(imx3, im, m, bicubic);
  double[][] im3 = ones(K, K);
  imcopy(im3, 0, 0, imx3, 0, 0, imx3.length, imx3[0].length);
  PImage img3 = createImage(K, K, GRAY);
  convert2Pimage(img3, im3);
  image(img3, 2.5*N+N+4, 0, K, K); 
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_GeoScale.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

double[][] scale(double[][] in, double scaleX, double scaleY) {
  int rows_in = in.length;
  int cols_in = in[0].length;
  int rows = (int) (rows_in * scaleX);
  int cols = (int) (cols_in * scaleY);
  double[][] m = zeros(3, 3);
  m[0][0] = ((double)rows_in - 1)/rows;
  m[1][1] = ((double)cols_in - 1)/cols;
  m[2][2] = 1;
  return m;
}
