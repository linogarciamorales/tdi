void convert2double(double[][] im, PImage img) {
  img.loadPixels();
  int rows = im.length;
  for(int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col= int(k/ rows)) {
    im[row][col] = map(red(img.pixels[k]), 0, 255, 0, 1.0);
  }
}

void convert2Pimage(PImage img, double[][] im) {
  int rows = im.length;
  for (int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col = int(k / rows)) {
    double Y = im[row][col]* 255;
    //img.pixels[k] = color(Y);
    img.set(row, col, color((float)Y));
  }
  img.updatePixels();
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
  double max = Double.MIN_VALUE; //lowest possible value of an double.
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
  double histMax = max(hist);
  stroke(255);
  // Draw half of the histogram (skip every second value)
  for (int i = 0; i < N; i += 2) {
    // Map i (from 0..img.width) to a location in the histogram (0..255)
    int which = int(map(i, 0, N, 0, 255));
    // Convert the histogram value to a location between 
    // the bottom and the top of the picture
    int y = int(map((float)hist[which], 0, (float)histMax, M, 0));
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
    tf[i] = (double)(hist[i] - Q)/ (P - Q);
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

      double sum = 0; // Kernel sum for this pixel
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

      double sum = 0; // Kernel sum for this pixel
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

void gaussianKernel(double[][] kernel, double sigma) {
  int W = kernel.length;
  int q = floor(kernel.length/2.0);
  double sum = 0;
  for (int v=-q, y=0; v <= q; v++, y=v+q) {
    for (int u=-q, x=0; u <= q; u++, x=u+q) {
      //double A = 1.0/(2.0* PI* sigma*sigma);
      double A = 1.0;
      kernel[x][y] = A* exp((float)(-(u*u+v*v)/(2*sigma*sigma)));
      sum += kernel[x][y];
    }
  }
  
  for(int x = 0; x < W; x++) {
    for(int y = 0; y < W; y++) {
      kernel[x][y] /= sum;
    }
  }
}

void gaussianKernel(double[][] kernel, double sigma, double K) {
  int W = kernel.length;
  int q = floor(kernel.length/2.0);
  double sigma2 = sigma*sigma;
  double K2 = K*K;
  double sum = 0;
  for (int v=-q, y=0; v <= q; v++, y=v+q) {
    for (int u=-q, x=0; u <= q; u++, x=u+q) {
      double A = 1.0/(2.0* PI* sigma2* K2);
      //double A = 1.0/(K2);
      kernel[x][y] = A* exp((float)(-(u*u+v*v)/(2*sigma2* K2)));
      sum += kernel[x][y];
    }
  }
  
  for(int x = 0; x < W; x++) {
    for(int y = 0; y < W; y++) {
      kernel[x][y] /= sum;
    }
  }
}

void LoG(double[][] kernel, double sigma) {
  int q = floor(kernel.length/2.0);
  double sigma2 = sigma*sigma;
  double sigma4 = sigma2*sigma2;
  double sum = 0;
  for (int v=-q, y=0; v <= q; v++, y=v+q) {
    for (int u=-q, x=0; u <= q; u++, x=u+q) {
      double r2 = u*u + v*v;
      double r = sqrt((float)r2);
      //kernel[x][y] = -(1.0-r/sigma2)* exp(-(r2)/(2*sigma2));
      kernel[x][y] = (r2-2*sigma2/sigma4)* exp((float)(-(r2)/(2*sigma2)));
      sum += kernel[x][y];
    }
  }
/*
  int W = kernel.length;
  for(int x = 0; x < W; x++) {
    for(int y = 0; y < W; y++) {
      kernel[x][y] /= sum;
    }
  }
*/
}
