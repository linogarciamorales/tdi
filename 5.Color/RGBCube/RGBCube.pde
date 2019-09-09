
PImage img_in;
int M, N;
PVector[] pos;
PVector[] colors;
int dim = 250;
boolean flag = true;

void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);  
  N = img_in.width; 
  M = img_in.height;
  pos = new PVector[N* M];
  colors = new PVector[N* M];

  size(512, 512, P3D);
}

void setup() {
  background(255);
  smooth();
  //noLoop();
  initcode();
} 

void draw() {
  if (flag) {
    code();
  }
}

void initcode() {
  stroke(0, 0, 0);
  strokeWeight(2);
  double[][] R = new double[N][M];
  convert2double(R, img_in, red);
  double[][] G = new double[N][M];
  convert2double(G, img_in, green);
  double[][] B = new double[N][M];
  convert2double(B, img_in, green);
  init(R, G, B);
}

void init (double[][] R, double[][] G, double[][] B) {
  for(int row = 0; row < N; row++) {
    for(int col = 0; col < M; col++) {
      int k = row* M + col;
      float r = (float) R[row][col];
      float g = (float) G[row][col];
      float b = (float) B[row][col];
      float x = map(r, 0, 1.0, 0, dim);
      float y = map(b, 0, 1.0, 0, dim);
      float z = map(g, 0, 1.0, 0, dim);
      pos[k] = new PVector(x, y, z);
      colors[k] = new PVector(r* 255, g* 255, b* 255);
    }
  }
}

void code() {
  background(255);
  translate(width/2, height/2);
  scale(1,-1,1); // so Y is up, which makes more sense in plotting
  rotateY(radians(frameCount));
  noFill();
  strokeWeight(1);
  box(dim);

  translate(-dim/2, -dim/2, -dim/2);
  for (int k= 0; k < pos.length; k++) {
    PVector p = pos[k];
    PVector c = colors[k];
    stroke(c.x, c.y, c.z);
    strokeWeight(5);
    point(p.x, p.y, p.z);
  }
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    flag = !flag;
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}
