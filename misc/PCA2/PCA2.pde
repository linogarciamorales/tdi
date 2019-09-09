
PImage img_in;
int M, N;
ArrayList<PVector> v = new ArrayList<PVector>();

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
  code();
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_PCA.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

void code() {
  stroke(0, 0, 0);
  strokeWeight(2);
  rect(0, 0, N, M);
  double[][] im = new double[N][M];
  convert2double(im, img_in);

  PImage img = createImage(N, M, GRAY);
  convert2Pimage(img, im);
  image(img, 0, 0, N, M);

  Matrix IM = new Matrix(im);
  SVD s = new SVD(IM);

  Matrix U = s.getU();
  Matrix S = s.getS();
  Matrix V = s.getV();
  
  int k = 64;
  Matrix U2 = U.getMatrix(0, N-1, 0, k-1);
  Matrix S2 = S.getMatrix(0, k-1, 0, k-1);
  Matrix V2 = V.getMatrix(0, M-1, 0, k-1);
  Matrix IM2 = U2.times(S2).times(V2.transpose());
  double[][] im2 = IM2.getArray();
  PImage img2 = createImage(N, M, GRAY);
  convert2Pimage(img2, im2);
  image(img2, N+2, 0, N, M);

  k = 32;
  Matrix U3 = U.getMatrix(0, N-1, 0, k-1);
  Matrix S3 = S.getMatrix(0, k-1, 0, k-1);
  Matrix V3 = V.getMatrix(0, M-1, 0, k-1);
  Matrix IM3 = U3.times(S3).times(V3.transpose());
  double[][] im3 = IM3.getArray();
  PImage img3 = createImage(N, M, GRAY);
  convert2Pimage(img3, im3);
  image(img3, 2*N+4, 0, N, M);
}
