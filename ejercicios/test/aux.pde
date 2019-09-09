String red = new String("red");
String green = new String("green");
String blue = new String("blue");
String hue = new String("hue");
String saturation = new String("saturation");
String brightness = new String("brightness");

void convert2double(double[][] im, PImage img, String chn) {
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
  int rows = im.length;
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

void convert2Pimage(PImage img, double[][] im_red, double[][] im_green, double[][] im_blue) {
  int rows = im_red.length;
  for (int k = 0, row = 0, col = 0; k < img.pixels.length; k++, row= k% rows, col = int(k / rows)) {
    double R = im_red[row][col]* 255;
    double G = im_green[row][col]* 255;
    double B = im_blue[row][col]* 255;
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
