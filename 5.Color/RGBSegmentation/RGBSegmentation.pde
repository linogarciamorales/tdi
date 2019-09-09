PImage img_in;

int N, M;
double[] avrc = {150/255.0, 110/255.0, 80/255.0};

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+ 2, M);
}

void setup() {
  colorMode(RGB);
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
    save("RGB_Seg.png");    
  }
}

void code() {
  image(img_in, 0, 0, N, M);
  double[][] R = new double[N][M];
  convert2double(R, img_in, red);

  double[][] G = new double[N][M];
  convert2double(G, img_in, green);

  double[][] B = new double[N][M];
  convert2double(B, img_in, blue);

  double[][] D = new double[N][M];
  D = dist(avrc, R, G, B);
  double[][] S = new double[N][M];
  brightnessThr(S, D, 0.2);

  PImage img_out = createImage(N, M, RGB);
  convert2Pimage(img_out, invert(S));
  image(img_out, N+2, 0, N, M);
}

double[][] dist(double[] avr, double[][] R, double[][] G, double[][] B) {
  int rows = R.length;
  int cols = R[0].length;
  double[][] d = new double[rows][cols];
  float ar = (float) avr[0];
  float ag = (float) avr[1];
  float ab = (float) avr[2];
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      float xr = (float) R[row][col];
      float xg = (float) G[row][col];
      float xb = (float) B[row][col];
      d[row][col] = sqrt((xr- ar)* (xr- ar) + (xg- ag)* (xg- ag) + (xb- ab)* (xb- ab));
    }
  }
  return d;
}
