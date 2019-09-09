PImage img_in;
double[][] im_in;
double[][] mask;

int N, M;
int W = 25; // odd

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(N, M);
}

void setup() {
  background(255);
  //test();

  image(img_in, 0, 0, N, M);
  stroke(255, 0, 0);
  noFill();
  int px = 100, py = 70;
  rect(px, py, W, W);
  
  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in);
  
  mask = new double[W][W];
  imcopy(mask, im_in, px, py, W, W);
  double[][] im_corr = new double[im_in.length+mask.length-1][im_in[0].length+mask[0].length-1];
  corr2D(im_corr, im_in, mask);
  
  PImage img_corr = createImage(im_corr.length, im_corr[0].length, GRAY);
  convert2Pimage(img_corr, im_corr);

  double[] pos = new double[2]; 
  Limits l = new Limits();
  l.find(im_corr);
  l.getMax(pos);
  println(pos[0], ",", pos[1]);

  stroke(255, 255, 0);
  rect((float)pos[0], (float)pos[1], W, W);

  int p = mask.length;
  int q = floor(p/2.0);
  stroke(0, 255, 0);
  rect((float)pos[0]-q, (float)pos[1]-q, W, W);

  smooth();
  noLoop();
} 

void draw() {
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("LP_BW_Sobel.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

void conv2D(double[][] im_out, double[][] im_in, double[][] kernel) {
  // Loop through every pixel in the image
  int p = kernel.length;
  int q = floor(p/2.0);
  double[][] im_aux = new double[im_in.length+2*(p-1)][im_in[0].length+2*(p-1)];
  double[][] im_aux_out = new double[im_in.length+(p-1)][im_in[0].length+(p-1)];
  //println(im_aux_out.length, im_aux_out[0].length);
  //println();
  
  imcopy(im_aux, p-1, p-1, im_in, 0, 0, im_in.length, im_in[0].length);

  for (int x = q; x < im_aux.length-q; x++) {  // Skip left and right edges
    for (int y = q; y < im_aux[0].length-q; y++) {   // Skip top and bottom edges

      float sum = 0; // Kernel sum for this pixel
      for (int u = -q; u <= q; u++) {
        for (int v = -q; v <= q; v++) {
          // println(x+u, y+v);
          // Calculate the adjacent pixel for this kernel point
          double val = im_aux[x+u][y+v];
          // Multiply adjacent pixels based on the kernel values
          sum += kernel[q-u][q-v] * val;
        }
      }
      // For this pixel in the new image, set the gray value
      // based on the sum from the kernel
      im_aux_out[x-q][y-q] = sum;
    }
  }
  imcopy(im_out, im_aux_out, 0, 0, im_out.length, im_out[0].length);
}

void corr2D(double[][] im_out, double[][] im_in, double[][] kernel) {
  // Loop through every pixel in the image
  int p = kernel.length;
  int q = floor(p/2.0);
  double[][] im_aux = new double[im_in.length+2*(p-1)][im_in[0].length+2*(p-1)];
  imcopy(im_aux, p-1, p-1, im_in, 0, 0, im_in.length, im_in[0].length);
  double[][] im_aux_out = new double[im_in.length+(p-1)][im_in[0].length+(p-1)];
  
  for (int x = q; x < im_aux.length-q; x++) {  // Skip left and right edges
    for (int y = q; y < im_aux[0].length-q; y++) {   // Skip top and bottom edges

      float sum = 0; // Kernel sum for this pixel
      for (int u = -q; u <= q; u++) {
        for (int v = -q; v <= q; v++) {
          // println(x+u, y+v);
          // Calculate the adjacent pixel for this kernel point
          double val = im_aux[x+u][y+v];
          // Multiply adjacent pixels based on the kernel values
          //sum += kernel[q-u][q-v] * val;
          sum += kernel[u+q][v+q] * val;
        }
      }
      // For this pixel in the new image, set the gray value
      // based on the sum from the kernel
      im_aux_out[x-q][y-q] = sum;
    }
  }
  imcopy(im_out, im_aux_out, 0, 0, im_out.length, im_out[0].length);
}

double[][] a = {{ 1, 2, 1 }, 
               { 3, 1, 2 }, 
               { 1, 0, 1 }};
double[][] h = {{ -1, -2, -1 }, 
               { 0, 0, 0 }, 
               { 1, 2, 1 }};

double[][] b = {{ 1, -1, 3, 2 }, 
               { 2, 1, 2, 4 }, 
               { 1, -1, 2, -2 }, 
               { 3, 1, 2, 2 }};

void test() {
  N = a.length+b.length-1;
  M = a[0].length+b[0].length-1;
  
  Limits l = new Limits();
  l.find(a);
  l.print();

  double[][] c = new double[N][M];
  corr2D(c, b, a);
  printarray(c);
  println();

  double[][] d = new double[N][M];
  conv2D(d, b, h);
  printarray(d);
}
