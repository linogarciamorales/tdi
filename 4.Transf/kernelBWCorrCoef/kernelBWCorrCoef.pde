PImage img_in;
double[][] im_in;
double[][] mask;

int N, M;
int W = 15; // odd

int[] hist_in = new int[256];
int[] hist_out = new int[256];
int[] hist_cum = new int[256];
double[] tfunc = new double[256];

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+2, M);
}

void setup() {
  //test();
  background(255);
  smooth();

  image(img_in, 0, 0, N, M);
  stroke(255, 0, 0);
  noFill();
  int px = 360, py = 70;
  rect(px, py, W, W);
  
  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in);
  
  mask = new double[W][W];
  imcopy(mask, im_in, px, py, W, W);
  double[][] im_corr = new double[im_in.length+mask.length-1][im_in[0].length+mask[0].length-1];
  corrCoeff2D(im_corr, im_in, mask);
  
  PImage img_corr = createImage(im_corr.length, im_corr[0].length, GRAY);
  convert2Pimage(img_corr, im_corr);
  image(img_corr, N+2, 0, N, M);
  double[] pos = new double[2]; 
  pos = findMaxIndex(im_corr);
  rect((float)pos[0]-(W-1)+N+2, (float)pos[1]-(W-1), W, W);
  println(px,",",py);
  println(pos[0]-(W-1),",",pos[1]-(W-1));

  noLoop();
} 

void draw() {
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("LP_BW_CorrCoef.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

double[][] h = {{ -1, -2, -1 }, 
               { 0, 0, 0 }, 
               { 1, 2, 1 }};

double[][] y = {{ 1, 2, 1 }, 
               { 3, 1, 2 }, 
               { 1, 0, 1 }};
               
double[][] b = {{ 1, -1, 3, 2 }, 
               { 2, 1, 2, 4 }, 
               { 1, -1, 2, -2 }, 
               { 3, 1, 2, 2 }};

double[][] c = {{ 2, 2, 3 }, 
               { 3, -1, 2 }, 
               { 1, 2, 1 }};

double[][] x = {{ 1, -1, 3, 2 }, 
               { 2, 1, 2, 4 }, 
               { 1, -1, 2, -2 },
               { 3, 1, 2, 2}};

/*
double[][] x = {{ 1, -1, 3 }, 
               { 2, 1, 2 }, 
               { 1, -1, 2 }};
*/
void test() {
  N = c.length+ c.length-1;
  M = c[0].length+ c[0].length-1;
  
  double[][] im_corr = new double[x.length+y.length-1][x[0].length+y[0].length-1];
  corrCoeff2D(im_corr, x, y);
  printarray(im_corr);
}
