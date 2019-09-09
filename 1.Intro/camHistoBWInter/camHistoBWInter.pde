import processing.video.*;

Capture cam;
int histMax = 0;
int[] hist = new int[256];

void setup() {
  size(640, 480);

  String[] cameras = Capture.list();
  
  if (cameras.length == 0) {
    println("There are no cameras available for capture.");
    exit();
  } else {
    println("Available cameras:");
    for (int i = 0; i < cameras.length; i++) {
      println(cameras[i]);
    }
    
    // The camera can be initialized directly using an 
    // element from the array returned by list():
    cam = new Capture(this, cameras[0]);
    cam.start();     
  }      
}

void draw() {
  loadPixels();
  if (cam.available() == true) {
    cam.read();
  }
  cam.filter(GRAY);
  image(cam, 0, 0);
  // The following does the same, and is faster when just drawing the image
  // without any additional resizing, transformations, or tint.
  //set(0, 0, cam);
  histo();
  histMax = max(hist);
  drawHisto();
}

void histo() {
  // Calculate the histogram
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      int bright = int(brightness(get(i, j)));
      hist[bright]++; 
    }
  }
}

void drawHisto() {
  stroke(255);
  // Draw half of the histogram (skip every second value)
  for (int i = 0; i < width; i += 2) {
    int which = int(map(i, 0, width, 0, 255));
    // Convert the histogram value to a location between 
    // the bottom and the top of the picture
    int y = int(map(hist[which], 0, histMax, height, 0));
    line(i, height, i, y);
  }
}
