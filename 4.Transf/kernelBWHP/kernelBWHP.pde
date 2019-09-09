PImage img_in;

int N, M;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(3*N+ 4, M);
}

void setup() {
  background(255);
  //image(img_in, 0, M+2, N, M);
  
  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in);
  double[][] im_out = new double[N][M];

  double Q = 9.0;
  double[][] kernel1 = {{ -1.0, -1.0, -1.0 }, 
                      {  -1.0, Q, -1.0 }, 
                      {  -1.0, -1.0, -1.0 }};
  Q = 25.0;                      
  double[][] kernel2 = {{ -1.0, -1.0, -1.0 , -1.0, -1.0 }, 
                      {  -1.0, -1.0, -1.0 , -1.0, -1.0 }, 
                      {  -1.0, -1.0, Q , -1.0, -1.0 },
                      {  -1.0, -1.0, -1.0 , -1.0, -1.0 },
                      {  -1.0, -1.0, -1.0 , -1.0, -1.0 }};
                      
  Q = 49.0;            
  double[][] kernel3 = {{ -1.0, -1.0, -1.0 , -1.0, -1.0 , -1.0, -1.0}, 
                      {  -1.0, -1.0, -1.0 , -1.0, -1.0 , -1.0, -1.0}, 
                      {  -1.0, -1.0, -1.0 , -1.0, -1.0 , -1.0, -1.0},
                      {  -1.0, -1.0, -1.0 , Q, -1.0 , -1.0, -1.0},
                      {  -1.0, -1.0, -1.0 , -1.0, -1.0 , -1.0, -1.0},
                      {  -1.0, -1.0, -1.0 , -1.0, -1.0 , -1.0, -1.0},
                      {  -1.0, -1.0, -1.0 , -1.0, -1.0 , -1.0, -1.0}};
                    
  filter(im_out, im_in, kernel1);
  PImage img_out1 = createImage(N, M, GRAY);
  convert2Pimage(img_out1, im_out);
  image(img_out1, 0, 0, N, M);
  
  filter(im_out, im_in, kernel2);
  PImage img_out2 = createImage(N, M, GRAY);
  convert2Pimage(img_out2, im_out);
  image(img_out2, N+2, 0, N, M);

  filter(im_out, im_in, kernel3);
  PImage img_out3 = createImage(N, M, GRAY);
  convert2Pimage(img_out3, im_out);
  image(img_out3, 2*N+4, 0, N, M);

  smooth();
  noLoop();
} 

void draw() {
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("LP_BW_HP.png");
  }
}

void filter(double[][] im_out, double[][] im_in, double[][] kernel) {
  // Loop through every pixel in the image
  int q = floor(kernel.length/2.0);
  for (int x = q; x < N-q; x++) {  // Skip left and right edges
    for (int y = q; y < M-q; y++) {   // Skip top and bottom edges

      double sum = 0; // Kernel sum for this pixel
      for (int v = -q; v <= q; v++) {
        for (int u = -q; u <= q; u++) {
          // Calculate the adjacent pixel for this kernel point
          double val = im_in[x+u][y+v];
          // Multiply adjacent pixels based on the kernel values
          sum += kernel[q-u][q-v] * val;
        }
      }
      // For this pixel in the new image, set the gray value
      // based on the sum from the kernel
      im_out[x][y] = sum;
    }
  }
}

void convert2double(double[][] im, PImage img) {
  int rows = im.length;
  for(int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col= int(k/ rows)) {
    im[row][col] = map(red(img.pixels[k]), 0, 255, 0, 1.0);
  }
}

void convert2Pimage(PImage img, double[][] im) {
  //img.loadPixels();
  int rows = im.length;
  for (int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col = int(k / rows)) {
    double Y = im[row][col]* 255;
    //img.pixels[k] = color(Y);
    img.set(row, col, color((float)Y));
  }
  //img.updatePixels();
}
