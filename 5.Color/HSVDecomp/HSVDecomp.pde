
PImage img_in;
int M, N;
void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle.png");
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
  //image(img_in, 0, M+2, N, M);

  double[][] H = new double[N][M];
  convert2double(H, img_in, hue);
  double[][] S = new double[N][M];
  convert2double(S, img_in, saturation);
  double[][] V = new double[N][M];
  convert2double(V, img_in, brightness);

  PImage imgh = createImage(N, M, GRAY);
  convert2Pimage(imgh, H);
  image(imgh, 0, 0, N, M);

  PImage imgs = createImage(N, M, GRAY);
  convert2Pimage(imgs, S);
  image(imgs, N+2, 0, N, M);

  PImage imgb = createImage(N, M, GRAY);
  convert2Pimage(imgb, V);
  image(imgb, 2*N+4, 0, N, M);
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("detalle_HSV.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}
