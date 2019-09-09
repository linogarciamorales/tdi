String bilinear = new String("bilinear");
String biquadratic = new String("biquadratic");
String bicubic = new String("bicubic");

double range(double x) {
  if(x < 0) {
    x = 0.0;
  }
  if(x > 1.0) {
    x = 1.0;
  }
  return x;
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

double[][] zeros(int rows, int cols) {
  double[][] out = new double[rows][cols];
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      out[row][col] = 0;
    }
  }
  return out;
}

double[][] ones(int rows, int cols) {
  double[][] out = new double[rows][cols];
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      out[row][col] = 1;
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

void imcopy(double[][] out, double[][] in) {
  int rows = out.length;
  int cols = out[0].length;
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      out[row][col] = in[row][col];
    }
  }
}

void imcopy(double[][] out, double[][] in, int start_row, int start_col, int end_row, int end_col) {
  for(int row = 0; row < end_row; row++) {
    for(int col = 0; col < end_col; col++) {
      out[row][col] = in[row+start_row][col+start_col];
    }
  }
}

void imcopy(double[][] out, int start_row_out, int start_col_out, 
            double[][] in, int start_row_in, int start_col_in, int end_row, int end_col) {
  for(int row = 0; row < end_row; row++) {
    for(int col = 0; col < end_col; col++) {
      out[row+start_row_out][col+start_col_out] = in[row+start_row_in][col+start_col_in];
    }
  }
}

void array2im(double[][] out, double[] in, int rows, int cols) {
  for (int k=0, row=0, col=0; k < rows* cols; k++, row = int(k / rows), col= k% rows) {
    out[row][col] = in[k];
  }
}

void im2array(double[] out, double[][] in) {
  int rows = in.length;
  int cols = in[0].length;
  for (int k = 0, row = 0, col = 0; k < rows*cols; k++, row = int(k / rows), col= k% rows) {
    out[k] = in[row][col];
  }
}

void imrow2array(double[] out, double[][] in, int row) {
  int cols = in[0].length;
  for (int col = 0; col < cols; col++) {
    out[col] = in[row][col];
  }
}

void imcol2array(double[] out, double[][] in, int col) {
  int rows = in.length;
  for (int row = 0; row < rows; row++) {
    out[row] = in[row][col];
  }
}

double[] findMaxIndex(double[][] im) {
  double[] pos = new double[2];
  pos[0] = pos[1] = -1;
  double max = Double.MIN_VALUE; //lowest possible value of an double.
  int rows = im.length;
  int cols = im[0].length;
  for(int row=0; row < rows; row++) {
    for(int col=0; col < cols; col++) {
      if(im[row][col] > max) {
        pos[0] = row;
        pos[1] = col;
        max = im[row][col];
      }
    }
  }
  return pos;
}

double findMax(double[][] im) {
  double[] pos = new double[2];
  pos[0] = pos[1] = -1;
  double max = Double.MIN_VALUE; //lowest possible value of an double.
  int rows = im.length;
  int cols = im[0].length;
  for(int row=0; row < rows; row++) {
    for(int col=0; col < cols; col++) {
      if(im[row][col] > max) {
        pos[0] = row;
        pos[1] = col;
        max = im[row][col];
      }
    }
  }
  return max;
}

void normalize(double[][] a) {
  double max = findMax(a);
  
  int rows = a.length;
  int cols = a[0].length;
  for(int row=0; row < rows; row++) {
    for(int col=0; col < cols; col++) {
      a[row][col] = a[row][col]/max;
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

void printarray(double[][] a, int row_start, int col_start, int rows, int cols) {
  for(int row = row_start; row < rows+row_start; row++) {
    double[] r = new double[cols];
    for(int col = col_start, k = 0; col < cols+col_start; col++, k++) {
        r[k] = a[row][col];
    }
    printarray(r);
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
  for (int i = 0; i < hist.length; i += 2) {
    // Map i (from 0..img.width) to a location in the histogram (0..255)
    int which = int(map(i, 0, hist.length, 0, 255));
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
  int rows = im.length;
  for(int k=0,row=0,col=0; k < img.pixels.length; k++, row = int(k / rows), col= k% rows) {
    int bright = int(brightness(img.pixels[k]));
    im[row][col] = tf[bright];
  }
}

void brightnessThr(double[][] out, double[][] in, double th) {
  // thresholding
  for(int k = 0, x = 0, y = 0; k < N* M; k++, x = k % N, y = int(k / N)) {
    if(in[x][y] <= th) {
      out[x][y] = 0.0;
    } else { 
      out[x][y] = 1.0;
    }
  }
}

void brightnessThrInv(double[][] out, double[][] in, double th) {
  // thresholding
  for(int k = 0, x = 0, y = 0; k < N* M; k++, x = k % N, y = int(k / N)) {
    if(in[x][y] >= th) {
      out[x][y] = 0.0;
    } else { 
      out[x][y] = 1.0;
    }
  }
}

double mean(double[][] a) {
  int rows = a.length;
  int cols = a[0].length;
  double mu = 0;
  for(int row=0; row < rows; row++) {
    for(int col=0; col < cols; col++) {
      mu += a[row][col];
    }
  }
  mu /= (rows* cols);
  return mu;
}

double var(double[][] a, double mu) {
  int rows = a.length;
  int cols = a[0].length;
  double d2 = 0;
  for(int row=0; row < rows; row++) {
    for(int col=0; col < cols; col++) {
      double xy = a[row][col]- mu;
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

  for (int row = q; row < im_aux.length-q; row++) {  // Skip left and right edges
    for (int col = q; col < im_aux[0].length-q; col++) {   // Skip top and bottom edges

      double sum = 0; // Kernel sum for this pixel
      for (int u = -q; u <= q; u++) {
        for (int v = -q; v <= q; v++) {
          // println(x+u, y+v);
          // Calculate the adjacent pixel for this kernel point
          double val = im_aux[row+u][col+v];
          // Multiply adjacent pixels based on the kernel values
          sum += kernel[q-u][q-v] * val;
        }
      }
      // For this pixel in the new image, set the gray value
      // based on the sum from the kernel
      im_aux_out[row-q][col-q] = sum;
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
  
  for (int row = q; row < im_aux.length-q; row++) {  // Skip left and right edges
    for (int col = q; col < im_aux[0].length-q; col++) {   // Skip top and bottom edges

      double sum = 0; // Kernel sum for this pixel
      for (int u = -q; u <= q; u++) {
        for (int v = -q; v <= q; v++) {
          // println(x+u, y+v);
          // Calculate the adjacent pixel for this kernel point
          double val = im_aux[row+u][col+v];
          // Multiply adjacent pixels based on the kernel values
          //sum += kernel[q-u][q-v] * val;
          sum += kernel[u+q][v+q] * val;
        }
      }
      // For this pixel in the new image, set the gray value
      // based on the sum from the kernel
      im_aux_out[row-q][col-q] = sum;
    }
  }
  imcopy(im_out, im_aux_out, 0, 0, im_out.length, im_out[0].length);
}

void substract(double[][] c, double[][] a, double[][] b) {
  int rows = c.length;
  int cols = c[0].length;
  for(int row=0; row < rows; row++) {
    for(int col=0; col < cols; col++) {
      c[row][col] = a[row][col]- b[row][col];
    }
  }
}

void substract(double[][] c, double cte) {
  int rows = c.length;
  int cols = c[0].length;
  for(int row=0; row < rows; row++) {
    for(int col=0; col < cols; col++) {
      c[row][col] = c[row][col]- cte;
    }
  }
}

double multiply(double[][] c, double[][] a, double[][] b) {
  int rows = c.length;
  int cols = c[0].length;
  double sum = 0;
  for(int row=0; row < rows; row++) {
    for(int col=0; col < cols; col++) {
      c[row][col] = a[row][col]* b[row][col];
      sum += c[row][col];
    }
  }
  return sum;
}

void multiply(double[][] a_re, double[][] a_im, double[][] b_re, double[][] b_im) {
  int rows = a_re.length;
  int cols = a_re[0].length;
  for(int row=0; row < rows; row++) {
    for(int col=0; col < cols; col++) {
      double xx = a_re[row][col] * b_re[row][col];
      double yy = a_im[row][col] * b_im[row][col];
      double xy = a_re[row][col] * b_im[row][col];
      double yx = a_im[row][col] * b_re[row][col];
      a_re[row][col] = xx - yy;
      a_im[row][col] = xy + yx;
    }
  }
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

  for (int row = q; row < im_aux.length-q; row++) {  // Skip left and right edges
    for (int col = q; col < im_aux[0].length-q; col++) {   // Skip top and bottom edges
      // get win
      double[][] win = new double[kernel.length][kernel[0].length];
      for (int u = -q; u <= q; u++) {
        for (int v = -q; v <= q; v++) {
          // println(x+u, y+v);
          // Calculate the adjacent pixel for this kernel point
          win[u+q][v+q] = im_aux[row+u][col+v];
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
      im_aux_out[row-q][col-q] = num/denum;
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

double[][] imreflect(double[][] in) {
  int rows = in.length;
  int cols = in[0].length;
  int N = max((int) pow(2, ceil((float)log2(rows))), (int) pow(2, ceil((float)log2(cols))));
  double[][] out = new double[N][N];
  
  imcopy(out, in, 0, 0, rows, cols);
  for(int row = rows, p = -2; row < N; row++, p--) {
    for(int col = 0; col < cols; col++) {
      out[row][col] = in[rows+p][col];
    }
  }
  for(int row = 0; row < rows; row++) {
    for(int col = cols, q = -2; col < N; col++, q--) {
      out[row][col] = in[row][cols+q];
    }
  }
  for(int row = rows, p = -2; row < N; row++, p--) {
    for(int col = cols, q = -2; col < N; col++, q--) {
      out[row][col] = in[cols+q][rows+p];
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

void gaussianNoiseGen(double[][] im_o, double[][] im_i, double gain) {
  // Iterate over image
  for (int x = 0; x < im_o.length; x++){
    for (int y = 0; y < im_o[0].length; y++) {
      double noise = gain* randomGaussian();
      im_o[x][y] = max(0, min(1.0, (float)(im_i[x][y]+ noise)));
    }
  }
}

double[][] shift(double[][] in) {
  int rows = in.length;
  int cols = in[0].length;
  int half_rows = rows/ 2;
  int half_cols = cols/ 2;
  double[][] out = zeros(rows, cols);
  imcopy(out, 0, 0, in, half_rows, half_cols, half_rows, half_cols);    // 0, 0
  imcopy(out, 0, half_cols, in, half_rows, 0, half_rows, half_cols);    // 0, 1
  imcopy(out, half_rows, 0, in, 0, half_cols, half_rows, half_cols);    // 1, 0
  imcopy(out, half_rows, half_cols, in, 0, 0, half_rows, half_cols);    // 1, 1
  return out;
}

double[] matrix_x_vector(double[][] M, double[] v) {
  int rows = M.length;
  int cols = M[0].length;
  double[] out = new double[rows];
  for(int row=0; row < rows; row++) {
    for(int col=0; col < cols; col++) {
      out[row] += M[row][col] * v[col];
    }
  }
  return out;
}

double linterp(double a, double b, double c) {
  return a + (b - a) * c;
}
 
double linterp(double a, double b, double c, double d) {
  return a + (b - c) * d;
}

double bqinterp(double a, double b, double c, double d, double q) {
  double d0 = a- b;
  double d2 = c- b;
  double d3 = d- b;
  double a0 = b;
  double a1 = -(1.0/3.0)* d0 + d2 - (1.0/6.0)* d3;
  double a2 = (1.0/2.0) * d0 + (1.0/2.0)* d2;
  double a3 = -(1.0/6.0)* d0 -(1.0/2.0)* d2 + (1.0/6.0)* d3;
  double q2 = q* q;
  double q3 = q2* q;
  return a0 + a1* q + a2* q2+ a3* q3;
}

double[] affine(double[][] m, double[] xy) {
  double kl[] = matrix_x_vector(m, xy);
  return kl;  
}

void transform(double[][] out, double[][] in, double[][] m, String method) {
  int rows_in = in.length;
  int cols_in = in[0].length;
  int rows = out.length;
  int cols = out[0].length;
  double Y;

  for (int x = 0; x < rows; x++) {
    for (int y = 0; y < cols; y++) {
      Y = 1;
      double[] xy = {x, y, 1};
      double[] kl = affine(m, xy);
      int xi = (int) kl[0];
      int yi = (int) kl[1];
      if ((xi >= 0 && xi <= rows_in-2) && (yi >= 0 && yi <= cols_in-2)) {
        double c00 = in[xi][yi];
        double c10 = in[xi+ 1][yi];
        double c01 = in[xi][yi+ 1];
        double c11 = in[xi+ 1][yi+ 1];
        double dx = kl[0] - xi;
        double dy = kl[1] - yi;
        if(method.equals(bilinear)) {
          Y = linterp(linterp(c00, c10, dx), linterp(c01, c11, dx), dy);
        }
        if (xi <= rows_in-3 && yi <= cols_in-3) {
          double c02 = in[xi][yi+ 2];
          double c12 = in[xi+ 1][yi+ 2];
          double c20 = in[xi+ 2][yi];
          double c21 = in[xi+ 2][yi+ 1];
          double c22 = in[xi+ 2][yi+ 2];
          if(method.equals(biquadratic)) {
            double qx = 0.5* dx* dx;  
            double qy = 0.5* dy* dy;  
            double ca = linterp(c10, c20, c00, dx) + linterp(0, c00, 2*c10-c20, qx);
            double cb = linterp(c11, c21, c01, dx) + linterp(0, c01, 2*c11-c21, qx);
            double cc = linterp(c12, c22, c02, dx) + linterp(0, c02, 2*c12-c22, qx);
            Y = linterp(cb, cc, ca, dy) + linterp(0, ca, 2*cb-cc, qy);
          }
          if (xi <= rows_in-4 && yi <= cols_in-4) {
            if(method.equals(bicubic)) {
              double c03 = in[xi][yi+ 3];
              double c13 = in[xi+ 1][yi+ 3];
              double c23 = in[xi+ 2][yi+ 3];
              double c33 = in[xi+ 3][yi+ 3];
              double c32 = in[xi+ 3][yi+ 2];
              double c31 = in[xi+ 3][yi+ 1];
              double c30 = in[xi+ 3][yi];
              double ca = bqinterp(c00, c10, c20, c30, dx);
              double cb = bqinterp(c01, c11, c21, c31, dx);
              double cc = bqinterp(c02, c12, c22, c32, dx);
              double cd = bqinterp(c03, c13, c23, c33, dx);
              Y = bqinterp(ca, cb, cc, cd, dy);
            }
          }
        }
      }
      out[x][y] = Y;
    }
  }
}

void dilation(double[][] im_out, double[][] im_in, double[][] kernel) {
  // Loop through every pixel in the image
  int p = kernel.length;
  int q = floor(p/2.0);
  double[][] im_aux = new double[im_in.length+2*(p-1)][im_in[0].length+2*(p-1)];
  double[][] im_aux_out = zeros(im_in.length+(p-1), im_in[0].length+(p-1));
  
  imcopy(im_aux, p-1, p-1, im_in, 0, 0, im_in.length, im_in[0].length);

  for (int row = q; row < im_aux.length-q; row++) {  // Skip left and right edges
    for (int col = q; col < im_aux[0].length-q; col++) {   // Skip top and bottom edges

      boolean sum = false; // Kernel sum for this pixel
      for (int u = -q; u <= q; u++) {
        for (int v = -q; v <= q; v++) {
          // Calculate the adjacent pixel for this kernel point
          boolean val = im_aux[row+u][col+v] != 0;
          boolean kval = kernel[q-u][q-v] !=0;
          // Multiply adjacent pixels based on the kernel values
          sum |= (kval & val);
        }
      }
      // For this pixel in the new image, set the gray value
      // based on the sum from the kernel
      im_aux_out[row-q][col-q] = sum?1.0:0.0;
    }
  }
  imcopy(im_out, im_aux_out, 1, 1, im_out.length, im_out[0].length);
}

void erosion(double[][] im_out, double[][] im_in, double[][] kernel) {
  dilation(im_out, invert(im_in), reflect(kernel));
  cmp(im_out);
}

void cmp(double[][] im) {
  int rows = im.length;
  int cols = im[0].length;
  
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      double val = im[row][col];
      if(val == 1.0) {
        im[row][col] = 0;
      } else {
        im[row][col] = 1;
      }
    }
  }
}

double[][] reflect(double[][] in) {
  int rows = in.length;
  int cols = in[0].length;
  double[][] out = new double[rows][cols];
  
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      out[row][col] = in[rows-row-1][cols-col-1];
    }
  }
  return out;
}

double[][] invert(double[][] in) {
  int rows = in.length;
  int cols = in[0].length;
  double[][] out = new double[rows][cols];
  
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      double val = in[row][col];
      if(val == 1.0) {
        out[row][col] = 0;
      } else {
        out[row][col] = 1;
      }
    }
  }
  return out;
}

double[][] xor(double[][] in1, double[][] in2) {
  int rows = in1.length;
  int cols = in1[0].length;
  double[][] out = new double[rows][cols];
  
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      boolean val1 = in1[row][col] != 0;
      boolean val2 = in2[row][col] != 0;
      boolean val = val1^val2;
      out[row][col] = val?1.0:0.0;
    }
  }
  return out;
}

void closing(double[][] im_out, double[][] im_in, double[][] m) {
  int rows = im_out.length;
  int cols = im_out[0].length;
  double[][] im_aux = new double[rows][cols];
  
  dilation(im_aux, im_in, m);
  erosion(im_out, im_aux, m);
}

void opening(double[][] im_out, double[][] im_in, double[][] m) {
  int rows = im_out.length;
  int cols = im_out[0].length;
  double[][] im_aux = new double[rows][cols];
  
  erosion(im_aux, im_in, m);
  dilation(im_out, im_aux, m);
}

PImage sum(PImage chn1, PImage chn2) {
  int rows = chn1.width;
  int cols = chn1.height;
  double[][] ch1 = new double[rows][cols];
  convert2double(ch1, chn1);
  double[][] ch2 = new double[rows][cols];
  convert2double(ch2, chn2);
  
  double[][] imo = new double[rows][cols];
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      imo[row][col] = range(ch1[row][col]+ ch2[row][col]);
    }
  }
  
  PImage img = createImage(rows, cols, RGB);
  convert2Pimage(img, imo);
  return img;
}

double[][] sub(double[][] ch1, double[][] ch2) {
  int rows = ch1.length;
  int cols = ch1[0].length;
  
  double[][] imo = new double[rows][cols];
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      imo[row][col] = range(ch1[row][col]- ch2[row][col]);
    }
  }
  return imo;
}

PImage sub(PImage chn1, PImage chn2) {
  int rows = chn1.width;
  int cols = chn1.height;

  double[][] ch1 = new double[rows][cols];
  convert2double(ch1, chn1);
  double[][] ch2 = new double[rows][cols];
  convert2double(ch2, chn2);
/*  
  double[][] imo = new double[rows][cols];
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      imo[row][col] = range(ch1[row][col]- ch2[row][col]);
    }
  }
*/
  double[][] imo = sub(ch1, ch2);
  PImage img = createImage(rows, cols, RGB);
  convert2Pimage(img, imo);
  return img;
}

PImage div(PImage chn1, PImage chn2) {
  int rows = chn1.width;
  int cols = chn1.height;
  double[][] ch1 = new double[rows][cols];
  convert2double(ch1, chn1);
  double[][] ch2 = new double[rows][cols];
  convert2double(ch2, chn2);
  
  double[][] imo = new double[rows][cols];
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      imo[row][col] = range(ch1[row][col]/ ch2[row][col]);
    }
  }
  
  PImage img = createImage(rows, cols, RGB);
  convert2Pimage(img, imo);
  return img;
}

PImage brightnessThr(PImage imgi,double th) {
  // thresholding
  int rows = imgi.width;
  int cols = imgi.height;
  double[][] imi = new double[rows][cols];
  convert2double(imi, imgi);

  double[][] imo = new double[rows][cols];
  for(int row = 0; row < rows; row++) {
    for(int col = 0; col < cols; col++) {
      if(imi[row][col] <= th) {
        imo[row][col] = 0.0;
      } else { 
        imo[row][col] = 1.0;
      }
    }
  }
  
  PImage img = createImage(rows, cols, RGB);
  convert2Pimage(img, imo);
  return img;
}
