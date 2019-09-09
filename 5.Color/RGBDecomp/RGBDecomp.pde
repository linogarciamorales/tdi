
PImage img_in;
int M, N;
void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle.png");
  img_in.resize(img_in.width/ 4, img_in.height/ 4);
  N = img_in.width; 
  M = img_in.height;
  size(3*N+4, 2*M+2);
}

void setup() {
  background(255);
  smooth();
  noLoop();
} 

void draw() {
  image(img_in, 0, M+2, N, M);

  double[][] imr = new double[N][M];
  convert2double(imr, img_in, red);
  double[][] img = new double[N][M];
  convert2double(img, img_in, green);
  double[][] imb = new double[N][M];
  convert2double(imb, img_in, blue);

  PImage imgr = createImage(N, M, GRAY);
  convert2Pimage(imgr, imr);
  image(imgr, 0, 0, N, M);

  PImage imgg = createImage(N, M, GRAY);
  convert2Pimage(imgg, img);
  image(imgg, N+2, 0, N, M);

  PImage imgb = createImage(N, M, GRAY);
  convert2Pimage(imgb, imb);
  image(imgb, 2*N+4, 0, N, M);
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("RGB_decomposition.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}
