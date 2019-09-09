PImage img_in;

int N, M;
int[] limits = new int[2];
int histMin = 0;
int histMax = 0;

double Q = 49;            
double[][] kernel = {{ 1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q}, 
                     { 1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q}, 
                     { 1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q},
                     { 1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q},
                     { 1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q},
                     { 1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q},
                     { 1.0/Q, 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q , 1.0/Q, 1.0/Q}};

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(3*N+ 4, M);
}

void setup() {
  smooth();
  noLoop();
}

void draw() {
  code();
}

void keyPressed() {
  if (key == 'q' || key == 'Q') {
    exit();
  }
  if (key == 's' || key == 'S') {
    save("KernelRGBHSV_LP.png");    
  }
}

void code() {
  image(img_in, 0, 0, N, M);
  
  double[][] R = new double[N][M];
  convert2double(R, img_in, red);
  double[][] Rf = new double[N][M];
  conv2D(Rf, R, kernel);
  
  double[][] G = new double[N][M];
  convert2double(G, img_in, green);
  double[][] Gf = new double[N][M];
  conv2D(Gf, G, kernel);

  double[][] B = new double[N][M];
  convert2double(B, img_in, blue);
  double[][] Bf = new double[N][M];
  conv2D(Bf, B, kernel);
  
  PImage img_RGB = createImage(N, M, RGB);
  convert2Pimage(img_RGB, Rf, Gf, Bf);
  image(img_RGB, N+2, 0, N, M);

  double[][] H = new double[N][M];
  convert2double(H, img_in, hue);
 
  double[][] S = new double[N][M];
  convert2double(S, img_in, saturation);

  double[][] V = new double[N][M];
  convert2double(V, img_in, brightness);
  double[][] Vf = new double[N][M];
  conv2D(Vf, V, kernel);
  //filter(Vf, V, kernel);


  double[][] R2 = new double[N][M];
  double[][] G2 = new double[N][M];
  double[][] B2 = new double[N][M];

  HSV2RGB(R2, G2, B2, H, S, Vf);
  PImage img_HSV = createImage(N, M, RGB);
  convert2Pimage(img_HSV, R2, G2, B2);
  image(img_HSV, 2*N+4, 0, N, M);
}

void filter(double[][] im_out, double[][] im_in, double[][] kernel) {
  // Loop through every pixel in the image
  int q = floor(kernel.length/2.0);
  int rows = im_out.length;
  int cols = im_out[0].length;
  for (int row = q; row < rows-q; row++) {  // Skip left and right edges
    for (int col = q; col < cols-q; col++) {   // Skip top and bottom edges

      double sum = 0; // Kernel sum for this pixel
      for (int v = -q; v <= q; v++) {
        for (int u = -q; u <= q; u++) {
          // Calculate the adjacent pixel for this kernel point
          double val = im_in[row+u][col+v];
          // Multiply adjacent pixels based on the kernel values
          sum += kernel[q-u][q-v] * val;
        }
      }
      // For this pixel in the new image, set the gray value
      // based on the sum from the kernel
      im_out[row][col] = sum;
    }
  }
}
