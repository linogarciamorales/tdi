PImage img_in;
double[][] im_in;

int[] hist_in = new int[256];
int[] hist_out = new int[256];

int N, M;
int[] limits = new int[2];
int histMin = 0, histMax = 0;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+ 2, M);
}

void setup() {
  colorMode(RGB);
  image(img_in, 0, 0, N, M);
  histo(img_in, hist_in);
  
  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in);
  int alpha = 4;
  double[][] im_out = new double[N][M];
  brightnessPosterisation(im_out, im_in, alpha);
  PImage img_out = createImage(N, M, GRAY);
  convert2Pimage(img_out, im_out);
  image(img_out, N+2, 0, N, M);
  histo(img_out, hist_out);

  smooth();
  noLoop();
}

void draw() {
  drawHisto(hist_in, 0);
  drawHisto(hist_out, N+2);
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("histo_detalle_BW_pos.png");    
  }
}
