PImage img_in;

int N, M;
int[] limits = new int[2];
int histMin = 0;
int histMax = 0;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+ 2, M);
}

void setup() {
  smooth();
  noLoop();
}

void draw() {
  //test();
  code();
}

void keyPressed() {
  if (key == 'q' || key == 'Q') {
    exit();
  }
  if (key == 's' || key == 'S') {
    save("HSV_EQHisto.png");    
  }
}

void test() {
  //double[] hsv = {62, 64, 99};
  double[] hsv = {300/360.0, 20/100.0, 99/100.0};
  double[] rgb = HSV2RGB(hsv);
  printarray(rgb);
}

void code() {
  image(img_in, 0, 0, N, M);
  
  double[][] H = new double[N][M];
  convert2double(H, img_in, hue);
 
  double[][] S = new double[N][M];
  convert2double(S, img_in, saturation);

  double[][] V = new double[N][M];
  convert2double(V, img_in, brightness);

  int[] hist_in = new int[256];
  int[] hist_cum = new int[256];
  double[] tfunc = new double[256];
  
  histo(hist_in, V);
  cumHisto(hist_cum, hist_in);
  normalization(tfunc, hist_cum);
  
  brightnessEQ(V, tfunc);

  double[][] R = new double[N][M];
  double[][] G = new double[N][M];
  double[][] B = new double[N][M];

  HSV2RGB(R, G, B, H, S, V);
  PImage img_out = createImage(N, M, RGB);
  convert2Pimage(img_out, R, G, B);
  image(img_out, N+2, 0, N, M);
}
