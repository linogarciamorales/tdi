
PImage img_in;
int M, N;
int K = 512;
double[][] m = zeros(3, 3);

String bilinear = new String("bilinear");
String biquadratic = new String("biquadratic");
String bicubic = new String("bicubic");

void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW_128_128.png");
  N = img_in.width; 
  M = img_in.height;
  size(3*K+4, K);
}

void setup() {
  background(255);
  smooth();
  noLoop();
} 

void draw() {
  //test();

  double[][] im = new double[N][M];
  convert2double(im, img_in);
  double[][] im1 = scale(im, 4, 4);
  transform(im1, im, bilinear);
  PImage img1 = createImage(K, K, GRAY);
  convert2Pimage(img1, im1);
  image(img1, 0, 0, K, K);

  double[][] im2 = scale(im, 4, 4);
  transform(im2, im, biquadratic);
  PImage img2 = createImage(K, K, GRAY);
  convert2Pimage(img2, im2);
  image(img2, K+2, 0, K, K);

  double[][] im3 = scale(im, 4, 4);
  transform(im3, im, bicubic);
  PImage img3 = createImage(K, K, GRAY);
  convert2Pimage(img3, im3);
  image(img3, 2*K+4, 0, K, K); 
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_GeoScaleComp.png");
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
  double[][] out = new double[rows][cols];
  m[0][0] = ((double)rows_in - 1)/rows;
  m[1][1] = ((double)cols_in - 1)/cols;
  m[2][2] = 1;
  return out;
}

void test(){
  double scaleX = 2.5, scaleY = 2.5;
  m[0][0] = scaleX;
  m[1][1] = scaleY;
  m[2][2] = 1;
  printarray(m);
  double[] v = {1, 1, 1};
  double[] o = matrix_x_vector(m, v);
  printarray(o);
  
}
