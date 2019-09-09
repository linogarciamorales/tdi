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

double[] findMaxIndex(double[][] im) {
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

double mean(double[][] a) {
  double mu = 0;
  for(int x=0; x < a.length; x++) {
    for(int y=0; y < a[0].length; y++) {
      mu += a[x][y];
    }
  }
  mu /= a.length*a[0].length;
  return mu;
}

double var(double[][] a, double mu) {
  double d2 = 0;
  for(int x=0; x < a.length; x++) {
    for(int y=0; y < a[0].length; y++) {
      double xy = a[x][y]- mu;
      d2 += xy* xy;
    }
  }
  d2 /= a.length*a[0].length;
  return d2;
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

void substract(double[][] c, double[][] a, double[][] b) {
  for(int x=0; x < c.length; x++) {
    for(int y=0; y < c[0].length; y++) {
      c[x][y] = a[x][y]- b[x][y];
    }
  }
}

void substract(double[][] c, double cte) {
  for(int x=0; x < c.length; x++) {
    for(int y=0; y < c[0].length; y++) {
      c[x][y] = c[x][y]- cte;
    }
  }
}

float multiply(double[][] c, double[][] a, double[][] b) {
  float sum = 0;
  for(int x=0; x < c.length; x++) {
    for(int y=0; y < c[0].length; y++) {
      c[x][y] = a[x][y]* b[x][y];
      sum += c[x][y];
    }
  }
  return sum;
}

void corrCoeff2D(double[][] im_out, double[][] im_in, double[][] kernel) {
  // Loop through every pixel in the image
  int p = kernel.length;
  int q = floor(p/2.0);
  double[][] im_aux = new double[im_in.length+2*(p-1)][im_in[0].length+2*(p-1)];
  imcopy(im_aux, p-1, p-1, im_in, 0, 0, im_in.length, im_in[0].length);
  double[][] im_aux_out = new double[im_in.length+(p-1)][im_in[0].length+(p-1)];
  
  double mu_kernel = mean(kernel);
  substract(kernel, mu_kernel);
  double var_kernel = var(kernel, mean(kernel));

  for (int x = q; x < im_aux.length-q; x++) {  // Skip left and right edges
    for (int y = q; y < im_aux[0].length-q; y++) {   // Skip top and bottom edges
      // get win
      double[][] win = new double[kernel.length][kernel[0].length];
      for (int u = -q; u <= q; u++) {
        for (int v = -q; v <= q; v++) {
          // println(x+u, y+v);
          // Calculate the adjacent pixel for this kernel point
          win[u+q][v+q] = im_aux[x+u][y+v];
        }
      }
      //printarray(win);
      double mu_win = mean(win);
      substract(win, mu_win);
      double var_win = var(win, mean(win));
  
      double num = multiply(win, win, kernel);
      double denum = sqrt((float)(var_win* var_kernel));
      //println(mu_win, mu_kernel, var_win, var_kernel); println();

      // For this pixel in the new image, set the gray value
      // based on the sum from the kernel
      im_aux_out[x-q][y-q] = num/denum;
    }
  }
  imcopy(im_out, im_aux_out, 0, 0, im_out.length, im_out[0].length);
}
