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
  colorMode(RGB);
  image(img_in, 0, 0, N, M);
  int[] hist_in = new int[256];
  int[] hist_cum = new int[256];
  double[] tfunc = new double[256];
  
  histo(hist_in, img_in, red);
  cumHisto(hist_cum, hist_in);
  normalization(tfunc, hist_cum);
  
  double[][] im_red_out = new double[N][M];
  brightnessEQ(im_red_out, img_in, tfunc, red);
  
  histo(hist_in, img_in, green);
  cumHisto(hist_cum, hist_in);
  normalization(tfunc, hist_cum);
  
  double[][] im_green_out = new double[N][M];
  brightnessEQ(im_green_out, img_in, tfunc, green);

  histo(hist_in, img_in, blue);
  cumHisto(hist_cum, hist_in);
  normalization(tfunc, hist_cum);
  
  double[][] im_blue_out = new double[N][M];
  brightnessEQ(im_blue_out, img_in, tfunc, blue);

  PImage img_out = createImage(N, M, RGB);
  convert2Pimage(img_out, im_red_out, im_green_out, im_blue_out);
  image(img_out, N+2, 0, N, M);

  smooth();
  noLoop();
}

void draw() {
  //drawHisto(hist_in, 0);
  //drawHisto(hist_cum, N+2);
  //drawHisto(hist_out, N+2);
  //save("histo_detalle_BW_nor.png");
}

void keyPressed() {
  if (key == 'q' || key == 'Q') {
    exit();
  }
  if (key == 's' || key == 'S') {
    save("RGB_EQHisto.png");    
  }
}
