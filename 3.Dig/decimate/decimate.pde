PImage img_in;

int N, M;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/plate.png");
  N = img_in.width; 
  M = img_in.height;
  size(N+ N+ 4, M);
}

void setup() {
  smooth();
  noLoop();
  background(255);
}

void draw() {
  code();
}

void keyPressed() {
  if (key == 'q' || key == 'Q') {
    exit();
  }
  if (key == 's' || key == 'S') {
    save("RGB_Dec.png"); 
  }
}

double Q = 5* 5;
double[][] kernel = {{ 1.0/Q, 1.0/Q, 1.0/Q, 1.0/Q, 1.0/Q }, 
                     { 1.0/Q, 1.0/Q, 1.0/Q, 1.0/Q, 1.0/Q }, 
                     { 1.0/Q, 1.0/Q, 1.0/Q, 1.0/Q, 1.0/Q }, 
                     { 1.0/Q, 1.0/Q, 1.0/Q, 1.0/Q, 1.0/Q }, 
                     { 1.0/Q, 1.0/Q, 1.0/Q, 1.0/Q, 1.0/Q }};

void code() {
  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in, red);
  PImage img_in = createImage(N, M, GRAY);
  convert2Pimage(img_in, im_in);
  image(img_in, 0, 0, N, M);

  double[][] D = diezm(im_in, 2);
  int N2 = D.length;
  int M2 = D[0].length;

  PImage img_out = createImage(N2, M2, GRAY);
  convert2Pimage(img_out, D);
  image(img_out, N+2, 0, N2, M2);

  double[][] im_out = new double[N][M];
  conv2D(im_out, im_in, kernel);
  double[][] E = diezm(im_out, 2);
  PImage img_out2 = createImage(N2, M2, GRAY);
  convert2Pimage(img_out2, E);
  image(img_out2, N+ N2+ 4, 0, N2, M2);
}

double[][] diezm(double[][] in, int df) {
  int rows = in.length;
  int cols = in[0].length;
  int new_rows = (int)(float(rows)/ float(df));
  int new_cols = (int)(float(cols)/ float(df));
  println(rows, cols, new_rows, new_cols);
  double[][] out = new double[new_rows][new_cols];
  int new_row = -1, new_col;
  double s = 0;
  for (int row = 0; row < rows-df; row+=df) {    
    new_row++;
    new_col = 0;
    for (int col = 0; col < cols-df; col+=df) {
      new_col++;
      s = 0;
      for (int u = 0; u < df; u++) {
        for (int v = 0; v < df; v++) {
          s = in[row+u][col+v];
        }
      }
      if( (new_row < new_rows) && (new_col < new_cols)) {
        out[new_row][new_col] = s;// df2;
      }
    }
  }
  return out;
}
