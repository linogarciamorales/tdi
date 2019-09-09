PImage img_in;
double[][] im_in;
int N, M;

void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 4, img_in.height/ 4);
  N = img_in.width; 
  M = img_in.height;
  size(N, M);
}

void setup() {
  smooth();
  noLoop();
} 

void draw() {
  image(img_in, 0, 0, N, M);

  double[][] im_in = new double[N][M];
  convert2double(im_in, img_in);

  int p = N*M;
  double[] v_in = new double[p];
  vectorize(v_in, im_in);
  sort(v_in);

  double mu_in = mean(im_in);
  double std_in = sqrt((float)(var(im_in, mu_in)));
  double m_in = median(v_in);
  println("mean: "+mu_in);
  println("std: "+std_in);
  double min_in = findMin(im_in);
  double max_in = findMax(im_in);
  println("min: "+min_in);
  println("median: "+m_in);
  println("max: "+max_in);
  double q1 = findQuartile(v_in, 25);
  println("Q1: "+q1);
  double q3 = findQuartile(v_in, 75);
  println("Q3: "+q3);
  double IQR = q3- q1;
  println("IQR: "+IQR);
  stroke(255, 0, 0);
  fill(255, 0, 0);
  int ymax = int(map((float) q3, 0, 1, M, 0));
  int x = floor(N/2);
  line(x, 0, x, ymax);
  text(String.format("%.2f", max_in), x+10, 10);
  int ymin = int(map((float) q1, 0, 1, M, 0));
  line(x, ymin, x, M);
  text(String.format("%.2f", min_in), x+10, M-10);
  noFill();
  text(String.format("%.2f", q1), x+15, ymin+5);
  text(String.format("%.2f", q3), x+15, ymax+5);
  rect(x-10, ymax, 20, ymin-ymax);
  int y = int(map((float) m_in, 0, 1, M, 0));
  line(x-10, y, x+10, y);
  text(String.format("%.2f", m_in), x+15, y+5);
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_Stat.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

void sort(double[] x) {
  for (int index = 0; index < x.length - 1; index++) {
    for (int j = index + 1; j < x.length; j++) {
      if (x[j] < x[index]) {
        double temp = x[j];
        x[j] = x[index];
        x[index] = temp;
      }
    }
  }
}

double median(double[] v) {
  int p = v.length;
  //sort(v);
  return v[(int)ceil(p/2.0)];
}

double median(double[][] win) {
  int p = win.length* win[0].length;
  double[] win_sort = new double[p];
  vectorize(win_sort, win);
  sort(win_sort);
  return win_sort[(int)ceil(p/2.0)];
}

void vectorize(double[] v, double[][] m) {
  int rows = m.length;
  int cols = m[0].length;
  for (int x = 0; x < rows; x++) {
    for (int y = 0; y < cols; y++) {
      v[y*rows + x] = m[x][y];   
    }
  }
}

double findMin(double[][] im) {
  double min = Double.MAX_VALUE; //lowest possible value of an double.
  for(int x=0; x < im.length; x++) {
    for(int y=0; y < im[0].length; y++) {
      if(im[x][y] < min) {
        min = im[x][y];
      }
    }
  }
  return min;
}

double findMax(double[][] im) {
  double max = Double.MIN_VALUE; //lowest possible value of an double.
  for(int x=0; x < im.length; x++) {
    for(int y=0; y < im[0].length; y++) {
      if(im[x][y] > max) {
        max = im[x][y];
      }
    }
  }
  return max;
}

double findQuartile(double[] v, double q) {
  //sort(v);
  int n = (int) Math.round(v.length * q / 100);      
  return v[n];
}
