PImage img; 
int N, M;
Boolean showFlag = false;

void settings() {
  //fullScreen();
  // Make a new instance of a PImage by loading an image file
  img = loadImage("../../img/mapa_antiguo_BW.png");
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
  loadPixels();
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
  int a = 40;
  int ofs = 2;
  loadPixels();

  for(int u = 0; u < Q; u++) { 
    for(int v = 0; v < Q; v++) { 
      int loc = (x + u- ofs) + (y + v- ofs) * width;
      if(loc > 0) {
        float value = brightness(pixels[loc]);
        fill(value);  // Set fill to white
        float X = x+ u* a;
        float Y = y+ v* a;
        rect(X, Y, a, a);
        fill(255);
        text(str(int(value)), X+ 0.5* a, Y+ 0.5* a);
      }
    }
  }
  //updatePixels();
}

void mouseClicked() {
  showFlag = !showFlag;
}
