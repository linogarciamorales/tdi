
PImage img_in;
int M, N;
ArrayList<PVector> v = new ArrayList<PVector>();

void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/silueta.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(N, M);
}

void setup() {
  background(255);
  smooth();
  noLoop();
} 

void draw() {
  //test1();
  //test2();
  code();
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_blobs.png");
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
  convert2Pimage(img, invert(im));
  image(img, 0, 0, N, M);
  findRegion(v, im);
  //fillRegion(img, v);
  //image(img, 0, 0, N, M);
  float[] r = getRegionLimits(v);
  stroke(146, 212, 241);
  noFill();
  rect(r[0], r[2], r[1]-r[0], r[3]-r[2]);
  float xc = 0.5* (r[0]+ r[1]);
  float yc = 0.5* (r[2]+ r[3]);
  noStroke();
  fill(146, 212, 241);
  ellipse(xc, yc, 5, 5);
  //float m00 = moments(v, 0, 0);
  float m00 = v.size();
  float m10 = moments(v, 1, 0);
  float m01 = moments(v, 0, 1);
  float xC = m10/ m00;
  float yC = m01/ m00;
  fill(255, 0, 0);
  ellipse(xC, yC, 5, 5);

  double mu20 = centralMoments(v, 2, 0, xC, yC);
  double mu02 = centralMoments(v, 0, 2, xC, yC);
  double mu11 = centralMoments(v, 1, 1, xC, yC);
  double[][] j = {{mu20, mu11}, {mu11, mu02}};
  
  Matrix J = new Matrix(j);
  SVD s = new SVD(J);
  Matrix S = s.getS();
  S.print(3, 3);
  Matrix V = s.getV();
  V.print(3, 3);
  float lambda1 = (float) S.get(0, 0);
  float lambda2 = (float) S.get(1, 1);
  float a = 2* sqrt(lambda1/ m00);
  float b = 2* sqrt(lambda2/ m00);
  float vx = (float) V.get(0, 0);
  float vy = (float) V.get(1, 0);
  float angle = atan(sqrt(vy)/ sqrt(vx));
  stroke(255, 0, 0);
  noFill();
  pushMatrix();
  translate(xC, yC);
  rotate(angle);
  ellipse(0, 0, 2* a, 2* b);
  popMatrix();
}

void findRegion(ArrayList<PVector> a, double[][] im) {
  int rows = im.length;
  int cols = im[0].length;
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      double val = im[row][col];
      if (val > 0) {
        a.add(new PVector(row, col));  
      }
    }
  }
}

void fillRegion(PImage img, ArrayList<PVector> a) {
  int K = a.size();
  for (int k = 0; k < K; k++) {
    PVector v = a.get(k);
    img.set((int) v.x, (int) v.y, color(255, 0, 0));
  }
}

float[] getRegionLimits(ArrayList<PVector> a) {
  int K = a.size();
  float[] l = {-1, -1, -1, -1};    // x.min, x.max, y.min, y.max
  for (int k = 0; k < K; k++) {
    PVector v = a.get(k);
    if (l[0] == -1 || v.x < l[0]) {
      l[0] = v.x;
    }
    if (l[1] == -1 || v.x > l[0]) {
      l[1] = v.x;
    } 
    
    if (l[2] == -1 || v.y < l[2]) {
      l[2] = v.y;
    }
    if (l[3] == -1 || v.y > l[3]) {
      l[3] = v.y;
    }
  }
  return l;
}

float moments(ArrayList<PVector> a, int p, int q) {
  float m = 0;
  int K = a.size();
  for (int k = 0; k < K; k++) {
    PVector v = a.get(k);
    m += pow(v.x, p)* pow(v.y, q);
  }
  return m;
}

float centralMoments(ArrayList<PVector> a, int p, int q, float xc, float yc) {
  float m = 0;
  int K = a.size();
  for (int k = 0; k < K; k++) {
    PVector v = a.get(k);
    m += pow(v.x - xc, p)* pow(v.y - yc, q);
  }
  return m;
}

void test1() {
  int M = 8, N = 5;
  Matrix B = Matrix.random(5, 3);
  Matrix A = Matrix.random(M, N).times(B).times(B.transpose());
  println("A = ");
  A.print(9, 6);
  println("A = U S V^T");
  println();
  SVD s = new SVD(A);
  print("U = ");
  Matrix U = s.getU();
  U.print(9, 6);
  print("Sigma = ");
  Matrix S = s.getS();
  S.print(9, 6);
  print("V = ");
  Matrix V = s.getV();
  V.print(9, 6);
  println("rank = " + s.rank());
  println("condition number = " + s.cond());
  println("2-norm = " + s.norm2());
}

double[][] a = {{1.0453, 1.1066, 0.8294},
                {1.1066, 1.4661, 1.0992},
                {0.8294, 1.0992, 0.8931}};

void test2() {
  Matrix A = new Matrix(a);
  eig e = new eig(A);
  print("U = ");
  Matrix D = e.getD();
  D.print(4, 4);
  print("V = ");
  Matrix V = e.getV();
  V.print(4, 4);
  Matrix M = V.times(V.transpose());
  print("M = ");
  M.print(4, 4);

}
