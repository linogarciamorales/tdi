import gab.opencv.*;

OpenCV cv;
Histogram grayHist, rHist, gHist, bHist;

PImage img;

int w = 200;
int h = 200;

void setup() {
  size(640, 420);
  img = loadImage("../../img/lena.jpg");
  cv = new OpenCV(this, img);

  grayHist = cv.findHistogram(cv.getGray(), 256);
  rHist = cv.findHistogram(cv.getR(), 256);
  gHist = cv.findHistogram(cv.getG(), 256);
  bHist = cv.findHistogram(cv.getB(), 256);
  smooth();
  noLoop();
}

void draw() {
  background(0);
  image(img, 10, 0,  img.width/2.5, img.height/2.5);
  
  stroke(125); noFill();  
  rect(430, 5, w, h);
  
  fill(125); noStroke();
  grayHist.draw(430, 5, w, h);

  stroke(255, 0, 0); noFill();  
  rect(10, 210, w, h);
  
  fill(255, 0, 0); noStroke();
  rHist.draw(10, 210, w, h);

  stroke(0, 255, 0); noFill();  
  rect(220, 210, w, h);
  
  fill(0, 255, 0); noStroke();
  gHist.draw(220, 210, w, h);

  stroke(0, 0, 255); noFill();  
  rect(430, 210, w, h);
  
  fill(0, 0, 255); noStroke();
  bHist.draw(430, 210, w, h);
}
