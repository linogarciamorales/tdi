PImage img_in;

int N, M;
int[] hist_in = new int[256];
int[] hist_cum = new int[256];
double[] tfunc = new double[256];
  

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
    save("RGB_HSV_EQHisto.png");    
  }
}

void code() {
  image(img_in, 0, 0, N, M);

  double[][] R = new double[N][M];
  convert2double(R, img_in, red);

  histo(hist_in, R);
  cumHisto(hist_cum, hist_in);
  normalization(tfunc, hist_cum);
  
  brightnessEQ(R, tfunc);
  
  double[][] G = new double[N][M];
  convert2double(G, img_in, green);

  histo(hist_in, G);
  cumHisto(hist_cum, hist_in);
  normalization(tfunc, hist_cum);
  
  brightnessEQ(G, tfunc);

  double[][] B = new double[N][M];
  convert2double(B, img_in, blue);

  histo(hist_in, B);
  cumHisto(hist_cum, hist_in);
  normalization(tfunc, hist_cum);
  
  brightnessEQ(B, tfunc);

  PImage img_RGB = createImage(N, M, RGB);
  convert2Pimage(img_RGB, R, G, B);
  image(img_RGB, N+2, 0, N, M);

  
  double[][] H = new double[N][M];
  convert2double(H, img_in, hue);
 
  double[][] S = new double[N][M];
  convert2double(S, img_in, saturation);

  double[][] V = new double[N][M];
  convert2double(V, img_in, brightness);

  histo(hist_in, V);
  cumHisto(hist_cum, hist_in);
  normalization(tfunc, hist_cum);
  
  brightnessEQ(V, tfunc);

  double[][] R2 = new double[N][M];
  double[][] G2 = new double[N][M];
  double[][] B2 = new double[N][M];

  HSV2RGB(R2, G2, B2, H, S, V);
  PImage img_HSV = createImage(N, M, RGB);
  convert2Pimage(img_HSV, R2, G2, B2);
  image(img_HSV, 2*N+4, 0, N, M);
}
