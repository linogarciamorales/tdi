PImage img_in;
double[][] kernel, kernel1, kernel2;

double[][] gradx = {{ 0, -1.0, 0 }, 
                   { 0, 1.0, 0 }, 
                   { 0, 0, 0 }};

double sigma1, sigma2;
int W;
int N, M;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+ 2, M);
}

void setup() {
  background(255);
  image(img_in, 0, 0, N, M);
  
  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in);
  
  double[][] im_out1 = new double[N][M];

  sigma2 = 2;
  sigma1 = .25;
  W = 5;
  double[][] kernel1 = new double[W][W];
  gaussianKernel(kernel1, sigma1);

  double[][] kernel2 = new double[W][W];
  gaussianKernel(kernel2, sigma2, 1.5);

  double[][] kernel = new double[W][W];
  sub(kernel, kernel1, kernel2);  // first derivative
  conv2D(im_out1, im_in, kernel);  // first derivative

/*
  // Alterative
  double[][] im_out2 = new double[N][M];
  double[][] im_out3 = new double[N][M];

  conv2D(im_out2, im_in, kernel1);  // first derivative
  conv2D(im_out3, im_in, kernel2);
  sub(im_out1, im_out2, im_out3);
*/

  PImage img_out = createImage(N, M, GRAY);
  convert2Pimage(img_out, im_out1);
  image(img_out, N+2, 0, N, M);
  
  smooth();
  noLoop();
} 

void draw() {
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("LP_BW_DoG.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

void sub(double[][] imo, double[][] imi2, double[][] imi1) {
  for(int x=0; x < imo.length; x++) {
    for(int y=0; y < imo[0].length; y++) {
      imo[x][y] = imi2[x][y] - imi1[x][y];
    }
  }
}
