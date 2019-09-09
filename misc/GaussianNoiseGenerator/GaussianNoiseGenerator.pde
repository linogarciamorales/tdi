PImage img_in;
double[][] im_in;
int N, M;

double Q = 16;
double[][] kernel = {{ 1/Q, 2/Q, 1/Q }, 
                    { 2/Q, 4/Q, 2/Q }, 
                    { 1/Q, 2/Q, 1/Q }};

void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(3*N+4, M);
}

void setup() {
  background(255);
  smooth();
  noLoop();
} 

void draw() {
  image(img_in, 0, 0, N, M);

  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in);

  double mu_in = mean(im_in);
  double std_in = sqrt((float)var(im_in, mu_in));
  println(std_in);
  float SNR = 30;
  double g = std_in* pow(10,-SNR/10);

  double[][] im_noised = new double[N][M];
  gaussianNoiseGen(im_noised, im_in, g);

  PImage img_noised = createImage(im_noised.length, im_noised[0].length, GRAY);
  convert2Pimage(img_noised, im_noised);
  image(img_noised, N+2, 0, N, M);

  double[][] im_out = new double[N][M];
  conv2D(im_out, im_noised, kernel);
  PImage img_out = createImage(N, M, GRAY);
  convert2Pimage(img_out, im_out);
  image(img_out, 2*N+4, 0, N, M);
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("LP_BW_GaussNoise.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

void gaussianNoiseGen(double[][] im_o, double[][] im_i, double gain) {
  // Iterate over image
  for (int x = 0; x < im_o.length; x++){
    for (int y = 0; y < im_o[0].length; y++) {
      double noise = gain* randomGaussian();
      im_o[x][y] = max(0, min(1.0, (float)(im_i[x][y]+ noise)));
    }
  }
}
