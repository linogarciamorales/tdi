int N=128;

void settings() {
  size(N, N);
}

void setup() { 
  background(0);
  smooth();
  fill(255);
  noStroke();
  
  ellipse(64, 64, 4, 4);
  
  save("LP.png");
 
  noLoop();
}

void draw() {
}
