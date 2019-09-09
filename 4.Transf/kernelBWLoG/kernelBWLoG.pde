PImage img_in;
double[][] im_in;
double[][] kernel, kernel1, kernel2;

double[][] laplace = {{ 0,  1, 0 }, 
                     { 1, -4, 1 }, 
                     { 0,  1, 0 }};

double sigma;
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
  
  double[][] im_out = new double[N][M];

  sigma = 1.0;
  W = 5;
  double[][] kernel = new double[W][W];
  LoG(kernel, sigma);

  filter(im_out, im_in, kernel);  // LoG
  //filter(im_out, im_in, laplace);

  PImage img_out = createImage(N, M, GRAY);
  convert2Pimage(img_out, im_out);
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
}

void sub(double[][] imo, double[][] imi2, double[][] imi1) {
  for(int x=0; x < imo.length; x++) {
    for(int y=0; y < imo[0].length; y++) {
      imo[x][y] = imi2[x][y] - imi1[x][y];
    }
  }
}
