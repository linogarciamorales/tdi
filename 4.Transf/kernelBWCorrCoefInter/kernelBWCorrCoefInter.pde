/*
Donde pincha el raton obtiene un kernel de WxW y 
lo busca en toda la imagen
*/

PImage img_in;
double[][] im_in;
double[][] mask;

int N, M;
int W = 15; // odd

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(N, M);
}

void setup() {
//  test();
  smooth();
  noFill();
  im_in = new double[N][M];
  convert2double(im_in, img_in);
  int px = int(random(0, N-2*W));
  int py = int(random(0, M-2*W));
  searchMask(px, py);
  //noLoop();
} 

void draw() {
  //searchMask(mouseX, mouseY);
}

void keyPressed() {  
  if (key == 's' || key == 'S') {
    save("LP_BW_Sobel.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

void mousePressed() {  
  searchMask(mouseX, mouseY);
}

void searchMask(int px, int py) {
  background(255);
  image(img_in, 0, 0, N, M);
  // red rect
  stroke(255, 0, 0);
  rect(px, py, W, W);
  println("Calculating...");
  
  mask = new double[W][W];
  imcopy(mask, im_in, px, py, W, W);
  double[][] im_corr = new double[im_in.length+mask.length-1][im_in[0].length+mask[0].length-1];
  corrCoeff2D(im_corr, im_in, mask);
  
  PImage img_corr = createImage(im_corr.length, im_corr[0].length, GRAY);
  convert2Pimage(img_corr, im_corr);

  println("Finished.");
  double[] pos = new double[2]; 
  pos = findMaxIndex(im_corr);
  println(px,",",py);
  println(pos[0]-(W-1),",",pos[1]-(W-1));
  
  // yellow rect
  stroke(255, 255, 0);
  rect((float)pos[0]-(W-1), (float)pos[1]-(W-1), W, W);
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
float[][] x = {{ 1, -1, 3 }, 
               { 2, 1, 2 }, 
               { 1, -1, 2 }};
*/
void test() {
  N = c.length+ c.length-1;
  M = c[0].length+ c[0].length-1;
  
  double[][] im_corr = new double[x.length+y.length-1][x[0].length+y[0].length-1];
  corrCoeff2D(im_corr, x, y);
  printarray(im_corr);

/*
  float mu_y = mean(y);
  substract(y, mu_y);
  float var_y = var(y, mean(y));
  printarray(y);
  
  float mu_x = mean(x);  
  substract(x, mu_x);
  float var_x = var(x, mean(x));
  printarray(x);
  
  float num = multiply(x, x, y);
  float denum = sqrt(var_x* var_y);
  println(mu_x, mu_y, var_x, var_y); println();
  println(num/denum);
  */
  
  /*
  float[][] q = new float[N][M];
  corr2DNorm(q, c, c);
  printarray(q);
  println();
*/
}
