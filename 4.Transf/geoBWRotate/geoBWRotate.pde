
PImage img_in;
int M, N;

String bilinear = new String("bilinear");
String biquadratic = new String("biquadratic");
String bicubic = new String("bicubic");

void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW_128_128.png");
  //img_in.resize(img_in.width/ 4, img_in.height/ 4);
  N = img_in.width; 
  M = img_in.height;
  size(4*N+6, M);
}

void setup() {
  background(255);
  smooth();
  noLoop();
} 

void draw() {
  double[][] im1 = new double[N][M];
  convert2double(im1, img_in);
  image(img_in, 0, 0, N, M);

  double[][] m;
  m = rotate(im1, PI/2, 0, 1);
  printarray(m);
  double[][] im2 = new double[N][M];
  transform(im2, im1, m, bilinear);
  //transform(im2, im1, m, biquadratic);
 // transform(im2, im1, m, bicubic);
  PImage img2 = createImage(N, M, GRAY);
  convert2Pimage(img2, im2);
  image(img2, N+2, 0, N, M);
  
  m = rotate(im1, PI, 1, 1);
  printarray(m);
  //double[][] m = rotate(im1, -PI/2, 1, 0);
  double[][] im3 = new double[N][M];
  transform(im3, im1, m, bilinear);
  //transform(im3, im1, m, biquadratic);
  //transform(im3, im1, m, bicubic);
  PImage img3 = createImage(N, M, GRAY);
  convert2Pimage(img3, im3);
  image(img3, 2*N+4, 0, N, M);

  m = rotate(im1, PI/4, 0, 0.5);
  printarray(m);
  double[][] im4 = new double[N][M];
  transform(im4, im1, m, bilinear);
  //transform(im4, im1, m, biquadratic);
  //transform(im4, im1, m, bicubic);
  PImage img4 = createImage(N, M, GRAY);
  convert2Pimage(img4, im4);
  image(img4, 3*N+6, 0, N, M);

}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_GeoRotate.png");
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

double[][] rotate(double[][] in, double angle, double x0, double y0) {
  int rows_in = in.length;
  int cols_in = in[0].length;
  int rows = rows_in;
  int cols = cols_in;
  double[][] m = zeros(3, 3);
  m[0][0] = cos((float)angle);
  m[0][1] = sin((float)angle);
  m[1][0] = -sin((float)angle);
  m[1][1] = cos((float)angle);
  m[0][2] = x0 * (rows - 1);     // point to rotate about;
  m[1][2] = y0 * (cols - 1);
  m[2][2] = 1;
  return m;
}
