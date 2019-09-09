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

void imcopy(double[][] mask, double[][] im_in, int u, int v, int w, int h) {
  for(int x = 0; x < w; x++) {
    for(int y = 0; y < h; y++) {
      mask[x][y] = im_in[x+u][y+v];
    }
  }
}

void imcopy(double[][] mask, int u1, int v1, double[][] im_in, int u2, int v2, int w, int h) {
  for(int x = 0; x < w; x++) {
    for(int y = 0; y < h; y++) {
      mask[x+u1][y+v1] = im_in[x+u2][y+v2];
    }
  }
}


double[] maxIndex(double[][] im) {
  double[] pos = new double[2];
  pos[0] = pos[1] = -1;
  double max = Double.MIN_VALUE; //lowest possible value of an float.
  for(int x=0; x < im.length; x++) {
    for(int y=0; y < im[0].length; y++) {
      if(im[x][y] > max) {
        pos[0] = x;
        pos[1] = y;
        max = im[x][y];
      }
    }
  }
  return pos;
}

void printarray(double[][] a) {
  for(int x = 0; x < a.length; x++) {
    for(int y = 0; y < a[0].length; y++) {
      print(a[x][y]+ " ");
    }
    println();
  }
  println();
}

void histo(int[] hist, PImage img) {
  // Calculate the histogram
  for(int k = 0; k < img.pixels.length; k++) {
      int bright = int(brightness(img.pixels[k]));
      hist[bright]++;
  }
}

void drawHisto(int[] hist, int xdisp) {
  float histMax = max(hist);
  stroke(255);
  // Draw half of the histogram (skip every second value)
  for (int i = 0; i < N; i += 2) {
    // Map i (from 0..img.width) to a location in the histogram (0..255)
    int which = int(map(i, 0, N, 0, 255));
    // Convert the histogram value to a location between 
    // the bottom and the top of the picture
    int y = int(map(hist[which], 0, histMax, M, 0));
    line(i+xdisp, M, i+xdisp, y);
  }
}

void cumHisto(int[] ohist, int[] ihist) {
  ohist[0] = ihist[0];
  for(int i = 1; i < ihist.length; i++) {
    ohist[i] = ihist[i] + ohist[i-1];
  }
}

void normalization(double[] tf, int[] hist) {
  // normalization
  int P = max(hist);
  int Q = min(hist);
  for(int i = 0; i < hist.length; i++) {
    tf[i] = (hist[i] - Q)/ (P - Q);
  }  
}

void brightnessEQ(double[][] im, PImage img, double[] tf) {
  // equalization
  for(int k=0,x=0,y=0; k < img.pixels.length; k++,x=k%im.length,y=int(k/im.length)) {
    int bright = int(brightness(img.pixels[k]));
    im[x][y] = tf[bright];
  }
}

static class Limits {
  static double max = Double.MIN_VALUE;
  static int xmax, ymax;
  static double min = Double.MAX_VALUE;
  static int xmin, ymin;

  void find(double[][] im) {
    for (int x=0; x < im.length; x++) {
      for (int y=0; y < im[0].length; y++) {
        if (im[x][y] > max) {
          max = im[x][y];
          xmax = x;
          ymax = y;
        } 
        if (im[x][y] < min) {
          min = im[x][y];
          xmin = x;
          ymin = y;
        }
      }
    }
  }
  
  void print() {
    println("min:", min,"(",xmin,",", ymin,")");
    println("max:", max,"(",xmax,",", ymax,")");
  }
  
  void getMax(double[] pos) {
    pos[0] = xmax;
    pos[1] = ymax;
  }
}
