
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

  double[][] imr = new double[N][M];
  convert2double(imr, img_in, red);
  double[][] img = new double[N][M];
  convert2double(img, img_in, green);
  double[][] imb = new double[N][M];
  convert2double(imb, img_in, blue);
  
  double[][] imy = new double[N][M];
  convert2double(imy, img_in, Y);

  PImage imgy = createImage(N, M, GRAY);
  convert2Pimage(imgy, imy);
  image(imgy, 0, 0, N, M);

  double[][] imcb = new double[N][M];
  convert2double(imcb, img_in, Cb);

  PImage imgcb = createImage(N, M, GRAY);
  convert2Pimage(imgcb, imcb);
  image(imgcb, N+2, 0, N, M);

  double[][] imcr = new double[N][M];
  convert2double(imcr, img_in, Cr);

  PImage imgcr = createImage(N, M, GRAY);
  convert2Pimage(imgcr, imcr);
  image(imgcr, 2*N+4, 0, N, M);
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("detalle_YCbCr.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}
