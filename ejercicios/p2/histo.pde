int[] histo(PImage img) {
  img.loadPixels();
  int[] hist = new int[256];
  // Calculate the histogram
  for(int k = 0; k < img.pixels.length; k++) {
      int bright = int(brightness(img.pixels[k]));
      hist[bright]++; 
  }
  return hist;
}

double[] histoNormalized(PImage img) {
  int rows = img.height;
  int cols = img.width;
  int[] hist = histo(img);
  int L = hist.length;
  double[] histn = new double[L];
  for(int k = 0; k < L; k++) {
      histn[k] = ((double) hist[k])/(rows*cols);
  }
  return histn;
}

void cumHisto(int[] ohist, int[] ihist) {
  ohist[0] = ihist[0];
  for(int i = 1; i < ihist.length; i++) {
    ohist[i] = ihist[i] + ohist[i-1];
  }
}

double[] cumHisto(double[] ihist) {
  double[] ohist = new double[ihist.length];
  ohist[0] = ihist[0];
  for(int i = 1; i < ihist.length; i++) {
    ohist[i] = ihist[i] + ohist[i-1];
  }
  return ohist;
}

void brightnessEQ(double[][] im, PImage img, double[] tf) {
  int rows = im.length;
  int cols = im[0].length;
  int MN = rows* cols;
  img.loadPixels();
  // equalization
  for(int k = 0, row = 0, col = 0; k < MN; k++, row = k % rows, col = int(k / rows)) {
    int bright = int(brightness(img.pixels[k]));
    im[row][col] = tf[bright];
  }
}

PImage histoEQ(PImage imgi) {
  PImage imgo = createImage(imgi.width, imgi.height, GRAY);
 
  double[] hist_in = histoNormalized(imgi);
  double[] hist_cum = cumHisto(hist_in);
  
  double[][] im_in = new double[imgi.width][imgi.height];
  convert2double(im_in, imgi);
  double[][] im_out = new double[imgi.width][imgi.height];
  brightnessEQ(im_out, imgi, hist_cum);
  convert2Pimage(imgo, im_out);
  return imgo;
}

int[] matchHistograms(int[] ht, int[] hs) {
  int K = ht.length;
  int[] fsh = new int[K];
  int[] cdft = new int[K];
  int[] cdfs = new int[K];
  cumHisto(cdft, ht);
  cumHisto(cdfs, hs);
  for(int k = 0; k < K; k++) {
    int j = K-1;
    do {
      fsh[k] = j--;
    } while((j>=0) && (cdft[k] <= cdfs[j]));
  }
  return fsh;
}

void matchHisto(PImage imgt, PImage imgs) {
  int[] ht = histo(imgt);
  int[] hs = histo(imgs);
  int[] fsh = matchHistograms(ht, hs);
  imgt.loadPixels();
  int rows = imgt.width;
  for (int k = 0, row = 0, col = 0; k < imgt.pixels.length; k++, row= k% rows, col = int(k / rows)) {
      int bright = int(brightness(imgt.pixels[k]));
      imgt.set(row, col, fsh[bright]);
  }
  imgt.updatePixels();
}
