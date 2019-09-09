void convert2double(double[][] im, PImage img) {
  int rows = im.length;
  img.loadPixels();
  for(int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col= int(k/ rows)) {
    im[row][col] = map(brightness(img.pixels[k]), 0, 255, 0, 1.0);
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

double[][] zeros(int rows, int cols) {
  double[][] out = new double[rows][cols];
  for(int y = 0; y < rows; y++) {
    for(int x = 0; x < cols; x++) {
      out[x][y] = 0;
    }
  }
  return out;
}

void arraycopy(double[] out, double[] in, int start, int end) {
  for(int x = 0; x+start < end; x++) {
      out[x] = in[x+start];
  }
}

void arraycopy(double[] out, int start_out, double[] in, int start_in, int end) {
  for(int x = 0; x < end; x++) {
      out[x+start_out] = in[x+start_in];
  }
}

void imcopy(double[][] out, double[][] in, int start_row, int start_col, int end_row, int end_col) {
  for(int x = 0; x < end_row; x++) {
    for(int y = 0; y < end_col; y++) {
      out[x][y] = in[x+start_row][y+start_col];
    }
  }
}

void imcopy(double[][] out, int start_row_out, int start_col_out, 
            double[][] in, int start_row_in, int start_col_in, int end_row, int end_col) {
  for(int x = 0; x < end_row; x++) {
    for(int y = 0; y < end_col; y++) {
      out[x+start_row_out][y+start_col_out] = in[x+start_row_in][y+start_col_in];
    }
  }
}

void array2im(double[][] out, double[] in, int rows, int cols) {
  for (int k=0, x=0, y=0; k < rows* cols; k++, y= k% rows, x= int(k/ rows)) {
    out[x][y] = in[k];
  }
}

void im2array(double[] out, double[][] in) {
  int n = in.length;
  int m = in[0].length;
  for (int k = 0, x = 0, y = 0; k < n*m; k++, y=k% n, x=int(k/n)) {
    out[k] = in[x][y];
  }
}

void imrow2array(double[] out, double[][] in, int row) {
  int n = in[0].length;
  for (int y = 0; y < n; y++) {
    out[y] = in[row][y];
  }
}

void imcol2array(double[] out, double[][] in, int col) {
  int n = in.length;
  for (int x = 0; x < n; x++) {
    out[x] = in[x][col];
  }
}

double[] findMaxIndex(double[][] im) {
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

double findMax(double[][] im) {
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
  return max;
}

void normalize(double[][] a) {
  double max = findMax(a);
  
  int rows = a.length;
  int cols = a[0].length;
  for(int x=0; x < rows; x++) {
    for(int y=0; y < cols; y++) {
      a[x][y] = a[x][y]/max;
    }
  }
}

void printarray(double[] a) {
  for(int x = 0; x < a.length; x++) {
    print(a[x]+ " ");
  }
  println();
}

void printarray(double[][] a) {
  for(int x = 0; x < a.length; x++) {
    double[] b = new double[a[0].length];
    imrow2array(b, a, x);
    printarray(b);
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
  int n = im.length;
  for(int k=0,x=0,y=0; k < img.pixels.length; k++,x=k%n,y=int(k/n)) {
    int bright = int(brightness(img.pixels[k]));
    im[x][y] = tf[bright];
  }
}

double mean(double[][] a) {
  int rows = a.length;
  int cols = a[0].length;
  double mu = 0;
  for(int x=0; x < rows; x++) {
    for(int y=0; y < cols; y++) {
      mu += a[x][y];
    }
  }
  mu /= (rows* cols);
  return mu;
}

double var(double[][] a, double mu) {
  int rows = a.length;
  int cols = a[0].length;
  double d2 = 0;
  for(int x=0; x < rows; x++) {
    for(int y=0; y < cols; y++) {
      double xy = a[x][y]- mu;
      d2 += xy* xy;
    }
  }
  d2 /= (rows* cols);
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

double multiply(double[][] c, double[][] a, double[][] b) {
  double sum = 0;
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

double[][] expanded(double[][] in) {
  int rows = in.length;
  int cols = in[0].length;
  int N = max((int) pow(2, ceil((float)log2(rows))), (int) pow(2, ceil((float)log2(cols))));
  double[][] out = zeros(N, N);
  imcopy(out, in, 0, 0, rows, cols);
  return out;
}

double[][] reflect(double[][] in) {
  int rows = in.length;
  int cols = in[0].length;
  int N = max((int) pow(2, ceil((float)log2(rows))), (int) pow(2, ceil((float)log2(cols))));
  double[][] out = new double[N][N];
  
  imcopy(out, in, 0, 0, rows, cols);
  for(int x = rows, p = -2; x < N; x++, p--) {
    for(int y = 0; y < cols; y++) {
      out[x][y] = in[rows+p][y];
    }
  }
  for(int x = 0; x < rows; x++) {
    for(int y = cols, q = -2; y < N; y++, q--) {
      out[x][y] = in[x][cols+q];
    }
  }
  for(int x = rows, p = -2; x < N; x++, p--) {
    for(int y = cols, q = -2; y < N; y++, q--) {
      out[x][y] = in[cols+q][rows+p];
    }
  }
  return out;
}

double log10(int n) {
  return log(n)/ log(10.0);
}

double log2(int n) {
  return log(n)/ log(2.0);
}

void butterflies(double[][] in, int n, int nu, boolean forward) {
  int n2 = n/ 2;
  int nu1 = nu - 1;
  double t_re, t_im, p, arg, c, s;

  // Here I check if I'm going to do the direct transform or the inverse transform.
  double constant;
  if (forward)
    constant = -2 * PI/ n;
  else
    constant = 2 * PI/ n;

  // First phase - calculation
  int k = 0;
  for (int l = 1; l <= nu; l++) {
    while (k < n) {
      for (int i = 1; i <= n2; i++) {
        p = bitreverseReference(k >> nu1, nu);
        // direct FFT or inverse FFT
        arg = constant * p;
        c = cos((float)arg);
        s = sin((float)arg);
        t_re = in[0][k + n2] * c - in[1][k + n2] * s;
        t_im = in[1][k + n2] * c + in[0][k + n2] * s;
        in[0][k + n2] = in[0][k] - t_re;
        in[1][k + n2] = in[1][k] - t_im;
        in[0][k] += t_re;
        in[1][k] += t_im;
        k++;
      }
      k += n2;
    }
    k = 0;
    nu1--;
    n2 /= 2;
  }
}

void butterflies(double[] re, double[] im, int n, int nu, boolean forward) {
  int n2 = n/ 2;
  int nu1 = nu - 1;
  double t_re, t_im, p, arg, c, s;

  // Here I check if I'm going to do the direct transform or the inverse transform.
  double constant;
  if (forward)
    constant = -2 * PI/ n;
  else
    constant = 2 * PI/ n;

  // First phase - calculation
  int k = 0;
  for (int l = 1; l <= nu; l++) {
    while (k < n) {
      for (int i = 1; i <= n2; i++) {
        p = bitreverseReference(k >> nu1, nu);
        // direct FFT or inverse FFT
        arg = constant * p;
        c = cos((float)arg);
        s = sin((float)arg);
        t_re = re[k + n2] * c - im[k + n2] * s;
        t_im = im[k + n2] * c + re[k + n2] * s;
        re[k + n2] = re[k] - t_re;
        im[k + n2] = im[k] - t_im;
        re[k] += t_re;
        im[k] += t_im;
        k++;
      }
      k += n2;
    }
    k = 0;
    nu1--;
    n2 /= 2;
  }
}

int bitreverseReference(int j, int nu) {
  // The reference bitreverse function.
  int j2;
  int j1 = j;
  int k = 0;
  for (int i = 1; i <= nu; i++) {
    j2 = j1 / 2;
    k = 2 * k + j1 - 2 * j2;
    j1 = j2;
  }
  return k;
}

void scramble(double[][] in, int n, int nu) {
  // Second phase - recombination
  int k = 0;
  int r;
  double t_re, t_im;
  while (k < n) {
    r = bitreverseReference(k, nu);
    if (r > k) {
      t_re = in[0][k];
      t_im = in[1][k];
      in[0][k] = in[0][r];
      in[1][k] = in[1][r];
      in[0][r] = t_re;
      in[1][r] = t_im;
    }
    k++;
  }
}

double[][] fft(double[][] in, boolean forward) {
    int n = in[0].length;

    // If n is a power of 2, then ld is an integer (_without_ decimals)
    double ld = log2(n);

    // Here I check if n is a power of 2. If exist decimals in ld, I quit
    // from the function returning null.
    if (((int) ld) - ld != 0) {
        System.out.println("The number of elements is not a power of 2.");
        return null;
    }

    double[][] t = new double[2][n];
    imcopy(t, in, 0, 0, in.length, in[0].length);

    int nu = (int) ld;
    butterflies(t, n, nu, forward);
    scramble(t, n, nu);

    double[][] out = new double[2][t[0].length];
    double radice = 1.0 / n;
    for (int i = 0; i < out[0].length; i++) {
        if (forward) {
          out[0][i] = t[0][i];
          out[1][i] = t[1][i];
        } else {
          out[0][i] = t[0][i] * radice;
          out[1][i] = t[1][i] * radice;
        }
    }
    return out;
}

void scramble(double[] re, double[] im, int n, int nu) {
  // Second phase - recombination
  int k = 0;
  int r;
  double t_re, t_im;
  while (k < n) {
    r = bitreverseReference(k, nu);
    if (r > k) {
      t_re = re[k];
      t_im = im[k];
      re[k] = re[r];
      im[k] = im[r];
      re[r] = t_re;
      im[r] = t_im;
    }
    k++;
  }
}

void fft(double[] re, double[] im, boolean forward) {
    int n = re.length;

    // If n is a power of 2, then ld is an integer (_without_ decimals)
    double ld = log2(n);

    // Here I check if n is a power of 2. If exist decimals in ld, I quit
    // from the function returning null.
    if (((int) ld) - ld != 0) {
        System.out.println("The number of elements is not a power of 2.");
        return;
    }

    double[] t_re = new double[n];
    double[] t_im = new double[n];
    arraycopy(t_re, re, 0, n);
    arraycopy(t_im, im, 0, n);

    int nu = (int) ld;
    butterflies(t_re, t_im, n, nu, forward);
    scramble(t_re, t_im, n, nu);

    double radice = 1.0 / n;
    for (int i = 0; i < n; i++) {
        if (forward) {
          re[i] = t_re[i];
          im[i] = t_im[i];
        } else {
          re[i] = t_re[i] * radice;
          im[i] = t_im[i] * radice;
        }
    }
}

void fft2(double[] re, double[] im, boolean forward) {
  // im cuadrada
  // parte real e imaginaria vectorizada
  int W = (int) sqrt(re.length); 
  double[] t_re = new double[W];
  double[] t_im = new double[W];
  
  // FFT filas
  for(int x = 0; x < W; x++) {
    int offset = x*W;
    arraycopy(t_re, 0, re, offset, W);
    arraycopy(t_im, 0, im, offset, W);
    fft(t_re, t_im, forward);
    arraycopy(re, offset, t_re, 0, W);
    arraycopy(im, offset, t_im, 0, W);
    //printarray(t_re);
  }

  //println();
  // FFT columnas
  for(int y = 0; y < W; y++) {
    int index = y;
    for(int x = 0; x < W; x++ ) {
      t_re[x] = re[index];
      t_im[x] = im[index];
      index += W;
    }
    //printarray(t_re);
    fft(t_re, t_im, forward);
    index = y;
    for(int x = 0; x < W; x++ ) {
      re[index] = t_re[x];
      im[index] = t_im[x];
      index += W;
    }    
  }
  //println();
}

void fft2(double[][] re, double[][] im) {
  int rows = re.length;
  int cols = re[0].length;
  
  double[] r = new double[rows* cols];
  im2array(r, re);  // vectoriza matriz
  double[] i = new double[rows* cols];
  im2array(i, im);                    // vectoriza matriz
  
  fft2(r, i, true);
  array2im(re, r, rows, cols);
  array2im(im, i, rows, cols);
}

void ifft2(double[][] re, double[][] im) {
  int rows = re.length;
  int cols = re[0].length;
  
  double[] r = new double[rows* cols];
  im2array(r, re);  // vectoriza matriz
  double[] i = new double[rows* cols];
  im2array(i, im);                    // vectoriza matriz
  
  fft2(r, i, false);
  array2im(re, r, rows, cols);
  array2im(im, i, rows, cols);
}

void abs(double[][] m, double[][] re, double[][] im) {
  int rows = re.length;
  int cols = re[0].length;
  for(int x=0; x < rows; x++) {
    for(int y=0; y < cols; y++) {
      m[x][y] = sqrt((float)(re[x][y]* re[x][y] + im[x][y]* im[x][y])); 
    }
  }  
}

void angle(double[][] a, double[][] re, double[][] im) {
  int rows = re.length;
  int cols = re[0].length;
  for(int x=0; x < rows; x++) {
    for(int y=0; y < cols; y++) {
      a[x][y] = atan2((float)im[x][y], (float)re[x][y]); 
    }
  }  
}
