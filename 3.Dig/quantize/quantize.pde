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
  background(0);
  smooth();
  noLoop();
} 

void draw() {
  code();
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_quantize.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

void code() {
  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in);

  double[][] im4 = pixelate(im_in, 4);
  PImage img1 = createImage(N, M, GRAY);
  convert2Pimage(img1, im4);
  image(img1, 0, 0, N, M);
  
  double[][] im8 = pixelate(im_in, 8);
  PImage img2 = createImage(N, M, GRAY);
  convert2Pimage(img2, im8);
  image(img2, N+2, 0, N, M);

  double[][] im16 = pixelate(im_in, 16);
  PImage img3 = createImage(N, M, GRAY);
  convert2Pimage(img3, im16);
  image(img3, 2*N+4, 0, N, M);
}

double[][] pixelate(double[][] in, int df) {
  int rows = in.length;
  int cols = in[0].length;
  double[][] out = new double[rows][cols];
  double s = 0;
  for (int row = 0; row < rows-df; row+=df) {    
    for (int col = 0; col < cols-df; col+=df) {
      s = in[row][col];      
      for (int u = 0; u < df; u++) {
        for (int v = 0; v < df; v++) {
          out[row+u][col+v] = s;
        }
      }
    }
  }
  return out;
}
