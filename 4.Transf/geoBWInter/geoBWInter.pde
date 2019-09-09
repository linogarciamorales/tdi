PImage img_in;
float[][] im_in;
int N, M;

// left mouse button + mouse drag = rotate
// right mouse button + mouse drag = translate
// mouse scroll wheel = scale
PVector movement = new PVector();
PVector rotation = new PVector();
PVector velocity = new PVector();
PVector position = new PVector(0, 0);

float rotationSpeed = 0.035;
float movementSpeed = 0.05;
float scaleSpeed = 0.05;
float fScale = 2;

void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/detalle_BW.png");
  img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(N, M, P2D);
  smooth();
}

void setup() {
  background(255);
} 

void draw() {
   if (mousePressed) {
    if (mouseButton==LEFT) {
      velocity.add( (pmouseY-mouseY) * 0.01, (mouseX-pmouseX) * 0.01, 0);
    }
    if (mouseButton==RIGHT) {
      movement.add( (mouseX-pmouseX) * movementSpeed, (mouseY-pmouseY) * movementSpeed, 0);
    }
  }
  velocity.mult(0.95);
  rotation.add(velocity);
  movement.mult(0.95);
  position.add(movement);
  
  background(0);
  //lights();

  //pushMatrix();
  //translate(position.x, position.y, position.z);
  translate(position.x, position.y);
  //rotateX(rotation.x*rotationSpeed);
  //rotateY(rotation.y*rotationSpeed);
  rotate(rotation.x*rotationSpeed);
  scale(fScale);
  
  image(img_in, 0, 0, N, M);
  //popMatrix();
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("geoBW_Inter.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

void mouseWheel(MouseEvent event) {
  float delta = event.getCount();
  fScale -= delta * scaleSpeed;
  fScale = max(0.5, fScale);
}
