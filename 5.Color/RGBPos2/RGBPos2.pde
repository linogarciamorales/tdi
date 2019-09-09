PImage img_in;
double[][] im_in;

int N, M;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+ 2, M);
}

float[][] colormap = {{255, 255, 0},
                    {0, 255, 255},
                    {255, 0, 255}};

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
    save("RGB_pos.png");    
  }
}

void code() {
  image(img_in, 0, 0, N, M);
  
  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in);
  double[] levels = {0, 0.25, 0.75, 1.0};
  PImage img_out = createImage(N, M, RGB);
  brightnessPosterisation(img_out, im_in, levels, colormap);
  image(img_out, N+2, 0, N, M);
}
