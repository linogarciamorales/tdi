String red = new String("red");
String green = new String("green");
String blue = new String("blue");
String hue = new String("hue");
String saturation = new String("saturation");
String brightness = new String("brightness");

String bilinear = new String("bilinear");
String biquadratic = new String("biquadratic");
String bicubic = new String("bicubic");

void convert2double(double[][] im, PImage img, String chn) {
  img.loadPixels();
  int rows = im.length;
  for(int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col= int(k/ rows)) {
    if(chn.equals(red)) {
      im[row][col] = map(red(img.pixels[k]), 0, 255, 0, 1.0);
    }
    if(chn.equals(green)) {
      im[row][col] = map(green(img.pixels[k]), 0, 255, 0, 1.0);
    }
    if(chn.equals(blue)) {
      im[row][col] = map(blue(img.pixels[k]), 0, 255, 0, 1.0);
    }
    // HSV
    if(chn.equals(hue)) {
      im[row][col] = map(hue(img.pixels[k]), 0, 255, 0, 1.0);
    }
    if(chn.equals(saturation)) {
      im[row][col] = map(saturation(img.pixels[k]), 0, 255, 0, 1.0);
    }
    if(chn.equals(brightness)) {
      im[row][col] = map(brightness(img.pixels[k]), 0, 255, 0, 1.0);
    }
  }
}

void convert2double(double[][] im, PImage img) {
  img.loadPixels();
  int rows = im.length;
  for(int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col= int(k/ rows)) {
    im[row][col] = map(brightness(img.pixels[k]), 0, 255, 0, 1.0);
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

double range(double x) {
  if(x < 0) {
    x = 0.0;
  }
  if(x > 1.0) {
    x = 1.0;
  }
  return x;
}

void convert2Pimage(PImage img, double[][] im) {
  //img.loadPixels();
  int rows = im.length;
  for (int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col = int(k / rows)) {
    double Y = range(im[row][col])* 255;
    //img.pixels[k] = color(Y);
    img.set(row, col, color((float)Y));
  }
  //img.updatePixels();
}

void convert2Pimage(PImage img, double[][] im_red, double[][] im_green, double[][] im_blue) {
  int rows = im_red.length;
  for (int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col = int(k / rows)) {
    double R = range(im_red[row][col])* 255;
    double G = range(im_green[row][col])* 255;
    double B = range(im_blue[row][col])* 255;
    img.set(row, col, color((float) R, (float) G, (float) B));
  }
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

double linterp(double s, double e, double t) {
  return s + (e - s) * t;
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

double[][] scale(double[][] in, double scaleX, double scaleY) {
  int rows_in = in.length;
  int cols_in = in[0].length;
  int rows = (int) (rows_in * scaleX);
  int cols = (int) (cols_in * scaleY);
  double[][] m = zeros(3, 3);
  m[0][0] = ((double)rows_in - 1)/rows;
  m[1][1] = ((double)cols_in - 1)/cols;
  m[2][2] = 1;
  return m;
}

double[] HSV2RGB(double[] hsv) {
  double[] rgb = {0, 0, 0};
  //double H = hsv[0]/ 360f;
  //double S = hsv[1]/ 100f;
  //double V = hsv[2]/ 100f;
  double H = hsv[0];
  double S = hsv[1];
  double V = hsv[2];
  
  double R = 0, G = 0, B = 0;
  
  if (S == 0) {
    R = V;
    G = V;
    B = V;
  } else {
    double hue = H * 6;
    if (hue == 6)
      hue = 0; // H must be < 1
    int i = (int) Math.floor(hue); // Or ... var_i  
                                                // floor( var_h )
    double var_1 = V * (1 - S);
    double var_2 = V * (1 - S * (hue - i));
    double var_3 = V * (1 - S * (1 - (hue - i)));

    if (i == 0) {
      R = V;
      G = var_3;
      B = var_1;
    } else if (i == 1) {
      R = var_2;
      G = V;
      B = var_1;
    } else if (i == 2) {
      R = var_1;
      G = V;
      B = var_3;
    } else if (i == 3) {
      R = var_1;
      G = var_2;
      B = V;
    } else if (i == 4) {
      R = var_3;
      G = var_1;
      B = V;
    } else {
      R = V;
      G = var_1;
      B = var_2;
    }
  }
          
  rgb[0] = R;
  rgb[1] = G;
  rgb[2] = B;
  return rgb;
}

void HSV2RGB(double[][] R, double[][] G, double[][] B, double[][] H, double[][] S, double[][] V) {
  int rows = R.length;
  int cols = R[0].length;
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      double hue = H[row][col];
      double sat = S[row][col];
      double val = V[row][col];
      double[] hsv = {hue, sat, val};
      double[] rgb = HSV2RGB(hsv);
      R[row][col] = rgb[0];
      G[row][col] = rgb[1];
      B[row][col] = rgb[2];
    }
  }
}

void brightnessThr(double[][] out, double[][] in, double th) {
  // thresholding
  int rows = out.length;
  int cols = out[0].length;
  for(int k = 0, row = 0, col = 0; k < rows* cols; k++, row= k% rows, col= int(k/ rows)) {
    if(in[row][col] <= th) {
      out[row][col] = 0.0;
    } else { 
      out[row][col] = 1.0;
    }
  }
}

double[][] dist(double[] avr, double[][] R, double[][] G, double[][] B) {
  int rows = R.length;
  int cols = R[0].length;
  double[][] d = new double[rows][cols];
  float ar = (float) avr[0];
  float ag = (float) avr[1];
  float ab = (float) avr[2];
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      float xr = (float) R[row][col];
      float xg = (float) G[row][col];
      float xb = (float) B[row][col];
      d[row][col] = sqrt((xr- ar)* (xr- ar) + (xg- ag)* (xg- ag) + (xb- ab)* (xb- ab));
    }
  }
  return d;
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
