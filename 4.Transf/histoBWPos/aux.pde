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

void histo(PImage img, int[] hist) {
  // Calculate the histogram
  for(int k = 0; k < img.pixels.length; k++) {
      int bright = int(brightness(img.pixels[k]));
      hist[bright]++;
  }
}

void drawHisto(int[] hist, int xdisp) {
  int histMax = max(hist);
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

int maxIndex(int[] hist) {
  int pos = -1;
  int max = Integer.MIN_VALUE; //lowest possible value of an int.
  for(int i=0; i < hist.length; i++) {
    if(hist[i] > max) {
      pos = i;
      max = hist[i];
    }
  }
  return pos;
}

int minIndex(int[] hist) {
  int pos = -1;
  int min = Integer.MAX_VALUE; //lowest possible value of an int.
  for(int i = 0; i < hist.length; i++) {
    if(hist[i] < min) {
      pos = i;
      min = hist[i];
    }
  }
  return pos;
}

void brightnessPosterisation(double[][] out, double[][] in, int n) {
  double a = 1.0/ n;
  double[] levels = new double[n+1];
  int l;
  
  // init
  for(l = 0; l < n; l++) {
    levels[l] = l * a;
  }
  levels[n] = 1.0;
  //println(levels);
  
  // posterisation
  for(int k = 0, x = 0, y = 0; k < N* M; k++, x = k % N, y = int(k / N)) {
    for(l = 0; l < n; l++) {
      if(in[x][y] > levels[l] && in[x][y] < levels[l + 1]) {
        out[x][y] = levels[l];
        break;
      }
    }
  }
}

void constrainLimits(double[][] a) {
  for(int k = 0, x = 0, y = 0; k < N* M; k++, x = k % N, y = int(k / N)) {
      if(a[x][y] > 1.0) {
        a[x][y] = 1.0;
      }
      if(a[x][y] < 0.0) {
        a[x][y] = 0.0;
      }
  }
}
