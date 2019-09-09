// Declaring a variable of type PImage
// PImage is a class available from the Processing core library.
PImage img; 
int N, M;
Boolean showFlag = false;

void settings() {
  //fullScreen();
  // Make a new instance of a PImage by loading an image file
  img = loadImage("../../img/mapa_antiguo.jpg");
  N = img.width/4; 
  M = img.height/4;
  size(N, M);
}

void setup() {
  rectMode(CORNER);  // Default rectMode is CORNER
  noStroke();
  smooth();
  textAlign(CENTER, CENTER);
}

void draw() {
  //background(0);
  // The image() function displays the image at a location
  // in this case the point (0,0).
  image(img, 0, 0, width, height);
  show(mouseX, mouseY);
}

void show(int x, int y) {
  if(!showFlag)
    return;
  int Q = 5;
  int a = 39;
  int ofs = 2;
  loadPixels();

  for(int u = 0; u < Q; u++) { 
    for(int v = 0; v < Q; v++) { 
      int loc = (x + u- ofs) + (y + v- ofs) * width;
      if(loc > 0) {
        float X = x+ u* a;
        float Y = y+ v* a;
        float red = red(pixels[loc]);
        fill(255, 0, 0);
        rect(X, Y, a, a/3);
        fill(0);
        text(str(int(red)), X+ 0.5* a, Y+ a/6);
        float green = green(pixels[loc]);
        fill(0, 255, 0);
        rect(X, Y+ a/3, a, a/3);
        fill(0);
        text(str(int(green)), X+ 0.5* a, Y+ 3*a/6);
        float blue = blue(pixels[loc]);
        fill(0, 0, 255);
        rect(X, Y+ 2*a/3, a, a/3);
        fill(255);
        text(str(int(blue)), X+ 0.5* a, Y+ 5*a/6);
      }
    }
  }
  //updatePixels();
}

void mouseClicked() {
  showFlag = !showFlag;
}
