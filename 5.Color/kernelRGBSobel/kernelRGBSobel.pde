PImage img_in;

int N, M, W;
double[][] hsobel = {{ 0.25, 0, -0.25 }, 
                    { 0.5, 0, -0.5 }, 
                    { 0.25, 0, -0.25 }};
double[][] vsobel = {{ -0.25, -0.5, -0.25 }, 
                    { 0, 0, 0 }, 
                    { 0.25, 0.5, 0.25 }};  

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(3*N+ 4, M);
}

void setup() {
  smooth();
  noLoop();
}

void draw() {
  code();
}

void keyPressed() {
  if (key == 'q' || key == 'Q') {
    exit();
  }
  if (key == 's' || key == 'S') {
    save("RGB_Sobel.png");    
  }
}

void code() {
  //image(img_in, 0, 0, N, M);

  double[][] R = new double[N][M];
  convert2double(R, img_in, red);
  double[][] Rx = new double[N][M];
  conv2D(Rx, R, hsobel);

  double[][] G = new double[N][M];
  convert2double(G, img_in, green);
  double[][] Gx = new double[N][M];
  conv2D(Gx, G, hsobel);

  double[][] B = new double[N][M];
  convert2double(B, img_in, blue);
  double[][] Bx = new double[N][M];
  conv2D(Bx, B, hsobel);
  double[][] gxx = grad(Rx, Gx, Bx); 

  PImage img_gxx = createImage(N, M, GRAY);
  convert2Pimage(img_gxx, gxx);
  image(img_gxx, 0, 0, N, M);  

  double[][] Ry = new double[N][M];
  conv2D(Ry, R, vsobel);
  double[][] Gy = new double[N][M];
  conv2D(Gy, G, vsobel);
  double[][] By = new double[N][M];
  conv2D(By, B, vsobel);
  double[][] gyy = grad(Ry, Gy, By); 
  double[][] gxy = grad(Rx, Gx, Bx, Ry, Gy, By); 

  PImage img_gyy = createImage(N, M, GRAY);
  convert2Pimage(img_gyy, gyy);
  image(img_gyy, N+2, 0, N, M);
  
  double[][] theta = ang(gxx, gyy, gxy);
  double[][] g = grad(gxx, gyy, gxy, theta);

  PImage img_g = createImage(N, M, GRAY);
  convert2Pimage(img_g, g);
  image(img_g, 2*N+4, 0, N, M);
}

double[][] grad(double[][] R, double[][] G, double[][] B) {
  int rows = R.length;
  int cols = R[0].length;
  double[][] out = new double[rows][cols];
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      double red = R[row][col];
      double green = G[row][col];
      double blue = B[row][col];
      out[row][col] = (red* red) + (green* green)+ (blue* blue);
    }
  }
  return out;  
}

double[][] grad(double[][] R1, double[][] G1, double[][] B1, double[][] R2, double[][] G2, double[][] B2) {
  int rows = R1.length;
  int cols = R1[0].length;
  double[][] out = new double[rows][cols];
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      double red1 = R1[row][col];
      double green1 = G1[row][col];
      double blue1 = B1[row][col];
      double red2 = R2[row][col];
      double green2 = G2[row][col];
      double blue2 = B2[row][col];
      out[row][col] = (red1* red2) + (green1* green2)+ (blue1* blue2);
    }
  }
  return out;  
}

double[][] ang(double[][] Gxx, double[][] Gyy, double[][] Gxy) {
  int rows = Gxx.length;
  int cols = Gxx[0].length;
  double[][] out = new double[rows][cols];
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      float gxx = (float) Gxx[row][col];
      float gxy = (float) Gxy[row][col];
      float gyy = (float) Gyy[row][col];
      out[row][col] = 0.5 * atan(2*gxy/(gxx-gyy));
    }
  }
  return out;  
}

double[][] grad(double[][] Gxx, double[][] Gyy, double[][] Gxy, double[][] Theta) {
  int rows = Gxx.length;
  int cols = Gxx[0].length;
  double[][] out = new double[rows][cols];
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      float gxx = (float) Gxx[row][col];
      float gxy = (float) Gxy[row][col];
      float gyy = (float) Gyy[row][col];
      float theta = (float) Theta[row][col];
      out[row][col] = sqrt(0.5*((gxx+ gyy) + (gxx- gxy)* cos(2* theta) + 2* gxy* sin(2* theta)));
    }
  }
  return out;  
}
