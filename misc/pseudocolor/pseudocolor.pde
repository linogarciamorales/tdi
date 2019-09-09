PImage img_in;
int N, M;

void settings() {
  // Make a new instance of a PImage by loading an image file
  //img_in = loadImage("../../img/detalle_BW.png");
  img_in = loadImage("../../img/detalle_BWEQ.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+ 2, M);
}

void setup() {
  colorMode(HSB, 360, 100, 100);
  smooth();
  noLoop();
  image(img_in, 0, 0, N, M);
}

void draw() {
  float[][] H = new float [N][M];
  float[][] S = new float [N][M];
  float[][] V = new float [N][M];
  
  setValue(H, img_in, 300, 0);
  setValue(S, 99);
  setValue(V, 99);
/*  
  float[][] R = new float[N][M];
  float[][] G = new float[N][M];
  float[][] B = new float[N][M];

  HSV2RGB(R, G, B, H, S, V);

  PImage img_out = createImage(N, M, RGB);
  convert2Pimage(img_out, R, G, B);
*/
  PImage pimg = createImage(N, M, HSB);
  setColor(pimg, H, S, V);
  image(pimg, N+2, 0, N, M);
}

void setColor(PImage img, float[][] h, float[][] s, float[][] v) {
  int rows = h.length;
  int cols = h[0].length;
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      color c = color(h[row][col], s[row][col], v[row][col]);
      //color c = color(200, 99, 99);
      img.set(row, col, c);
    }
  }
  img.updatePixels();
}

void setValue(float[][] m, PImage pimg, float v1, float v2) {
  int rows = m.length;
  int cols = m[0].length;
//  for(int k = 0, row = 0, col = 0; k < pimg.pixels.length; k++, row= k% rows, col= int(k/ rows)) {
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      float Y = brightness(pimg.get(row, col));
      //println(Y);
      //float Y = brightness(pimg.pixels[k]);
      float h = map(Y, 0, 255, v1, v2);
      m[row][col] = h;
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
  int cols = m[0].length;
  int rows = m.length;
  for(int col = 0; col < cols; col++) {
    float val = map(col, 0, cols-1, v1, v2);
    for(int row = 0; row < rows; row++) {
      m[row][col] = val;
    }
  }
}
