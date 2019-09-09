int N, M;

void setup() {
  size(512, 512);
  N = width;
  M = height;
  colorMode(HSB, 360, 100, 100);
  smooth();
  noLoop();
}

void draw() {
  float[][] H = new float [N][M];
  float[][] S = new float [N][M];
  float[][] B = new float [N][M];
  
  setValue(H, 0, 359);
  //setValue(H, 300);
  setValue(S, 99);
  setValue(B, 99);
  PImage img_out = createImage(N, M, HSB);
  setColor(img_out, H, S, B);
  image(img_out, 0, 0);
}

void setColor(PImage img, float[][] h, float[][] s, float[][] v) {
  int rows = img.width;
  int cols = img.height;
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      color c = color(h[row][col], s[row][col], v[row][col]);
      img.set(row, col, c);
    }
  }
}

void setValue(float[][] m, int v) {
  int rows = m.length;
  int cols = m[0].length;
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      m[row][col] = v;
    }
  }
}

void setValue(float[][] m, float v1, float v2) {
  int rows = m.length;
  int cols = m[0].length;
  for(int col = 0; col < cols; col++) {
    float val = map(col, 0, cols-1, 0, 359);
    for(int row = 0; row < rows; row++) {
      m[col][row] = val;
    }
  }
}
