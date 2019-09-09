PImage img_in;

int N, M;
int K = 3;
boolean flag = false;
double threshold = 0.07;

void settings() {
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(2*N+ 2, M);
}

int[][] colormap1 = {{150, 110, 80},
                     {190, 190, 190},
                     {75, 65, 45}};
  
int[][] colormap2 = {{255, 255, 0},
                     {0, 255, 255},
                     {255, 0, 255}};

void setup() {
  colorMode(RGB);
  noStroke();
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
    save("RGB_posKMeans.png");    
  }
}

void code() {
  image(img_in, 0, 0, N, M);
  double[][] R = new double[N][M];
  convert2double(R, img_in, red);

  double[][] G = new double[N][M];
  convert2double(G, img_in, green);

  double[][] B = new double[N][M];
  convert2double(B, img_in, blue);
  
  ArrayList<double[][]> rgb = new ArrayList<double[][]>();
  rgb.add(R);
  rgb.add(G);
  rgb.add(B);

  ArrayList<double[]> C_old = new ArrayList<double[]>();
  initArray(C_old, K);
  ArrayList<double[]> C_new = new ArrayList<double[]>();
  
  double[] C1 = {colormap1[0][0]/255.0, colormap1[0][1]/255.0, colormap1[0][2]/255.0};
  double[] C2 = {colormap1[1][0]/255.0, colormap1[1][1]/255.0, colormap1[1][2]/255.0};
  double[] C3 = {colormap1[2][0]/255.0, colormap1[2][1]/255.0, colormap1[2][2]/255.0};
  C_new.add(C1);
  C_new.add(C2);
  C_new.add(C3);
  
  leyend(N-75, 10, 20, colormap1);
  
  ArrayList<double[][]> D = new ArrayList<double[][]>();
  initMatrix(D, K);

  ArrayList<ArrayList<PVector>> C = new ArrayList<ArrayList<PVector>>();
  initArrayToArray(C, K);

  while(!flag) {
    updateCentroids(C_old, C_new);
    dist(D, C_old, rgb);
  
    clean(C);
    pixel2cluster(C, D);
    getCentroids(C_new, C, rgb);
    flag = dist(C_new, C_old, threshold);
    println(flag);
  }  
  for (int k=0; k < 3; k++) {
    fillRegion(img_in, C.get(k), colormap1[k]);
  }

  image(img_in, N+2, 0, N, M);
  leyend((2*N+2)-75, 10, 20, colormap1);
}

void leyend(int row, int col, int b, int[][] colormap) {
  int L = colormap.length;
  for(int k = 0; k < L; k++) {
    fill(colormap[k][0], colormap[k][1], colormap[k][2]);
    rect(row+ k*(b+ 2), col, b, b);
  }
}
