PImage img_in;

int[] hist_in = new int[256];
int[] hist_out = new int[256];

int N, M;
int[] limits = new int[2];
int histMin = 0, histMax = 0;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+ 2, M);
}

float[][] colormap = {{255, 0, 240},
                      {90, 255, 255},
                      {255, 255, 90},
                      {0, 255, 255},
                      {255, 240, 0},
                      {255, 0, 255}};

void setup() {
  colorMode(RGB);
  image(img_in, 0, 0, N, M);
  //histo(hist_in, img_in);
  
  double[][] imr = new double[N][M];
  convert2double(imr, img_in, red);
  double[][] img = new double[N][M];
  convert2double(img, img_in, green);
  double[][] imb = new double[N][M];
  convert2double(imb, img_in, blue);

  double[][] imri = complement(imr);
  double[][] imgi = complement(img);
  double[][] imbi = complement(img);
  PImage img_out = createImage(N, M, RGB);
  convert2Pimage(img_out, imri, imgi, imbi);
  image(img_out, N+2, 0, N, M);
  //histo(hist_out, img_out);

  smooth();
  noLoop();
}

void draw() {
  //drawHisto(hist_in, 0);
  //drawHisto(hist_out, N+2);
}

void keyPressed() {
  if (key == 'q' || key == 'Q') {
    exit();
  }
  if (key == 's' || key == 'S') {
    save("RGB_comp.png");    
  }
}
