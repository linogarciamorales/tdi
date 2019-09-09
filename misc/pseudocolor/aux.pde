void convert2double(double[][] im, PImage img) {
  int rows = im.length;
  for(int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col= int(k/ rows)) {
    im[row][col] = map(red(img.pixels[k]), 0, 255, 0, 1.0);
  }
}

void convert2Pimage(PImage img, double[][] im) {
  img.loadPixels();
  int rows = im.length;
  for (int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col = int(k / rows)) {
    double Y = im[row][col]* 255;
    //img.pixels[k] = color(Y);
    img.set(row, col, color((float)Y));
  }
  img.updatePixels();
}

void convert2Pimage(PImage img, float[][] im_red, float[][] im_green, float[][] im_blue) {
  int rows = im_red.length;
  for (int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col = int(k / rows)) {
    float R = im_red[row][col]* 255;
    float G = im_green[row][col]* 255;
    float B = im_blue[row][col]* 255;
    img.set(row, col, color(R, G, B));
  }
  img.updatePixels();
}

void brightnessEQ(double[][] im, PImage img, double[] tf) {
  // equalization
  for(int k = 0, x = 0, y = 0; k < img.pixels.length; k++, x = k % im.length, y = int(k / im.length)) {
    int bright = int(brightness(img.pixels[k]));
    im[x][y] = tf[bright];
  }
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

void normalization(double[] tf, int[] ihist) {
  // normalization
  int P = max(ihist);
  int Q = min(ihist);
  for(int i = 0; i < ihist.length; i++) {
    tf[i] = (double)(ihist[i] - Q)/ (P - Q);
  }  
}

float[] HSV2RGB(float[] hsv) {
  float[] rgb = {0, 0, 0};
  //double H = hsv[0]/ 360f;
  //double S = hsv[1]/ 100f;
  //double V = hsv[2]/ 100f;
  float H = hsv[0];
  float S = hsv[1];
  float V = hsv[2];
  
  float R = 0, G = 0, B = 0;
  
  if (S == 0) {
    R = V;
    G = V;
    B = V;
  } else {
    float hue = H * 6;
    if (hue == 6)
      hue = 0; // H must be < 1
    int i = (int) Math.floor(hue); // Or ... var_i  
                                                // floor( var_h )
    float var_1 = V * (1 - S);
    float var_2 = V * (1 - S * (hue - i));
    float var_3 = V * (1 - S * (1 - (hue - i)));

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

void HSV2RGB(float[][] R, float[][] G, float[][] B, float[][] H, float[][] S, float[][] V) {
  int rows = R.length;
  int cols = R[0].length;
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      float hue = H[row][col];
      float sat = S[row][col];
      float val = V[row][col];
      float[] hsv = {hue, sat, val};
      float[] rgb = HSV2RGB(hsv);
      R[row][col] = rgb[0];
      G[row][col] = rgb[1];
      B[row][col] = rgb[2];
    }
  }
}
