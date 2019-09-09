class Blackhole {
  PVector location;
  PVector center;
 
  float radius = 0;
  int numPoints = 4;
  float angle = TWO_PI/(float) numPoints;
  float phi = PI/10f;
  float [][] xyArray;
 
  Blackhole(float w, float h, float r, int n) {
    location = new PVector(w, h);
    numPoints = n;
    radius = r;
  }
 
  void divide() {
    xyArray = new float [numPoints][3]; 
    for (int i=0; i<numPoints; i++) { 
      float x = radius* sin(angle* i + phi)+ location.x;
      float y = radius* cos(angle* i + phi)+ location.y;
      xyArray[i][0] = x; 
      xyArray[i][1] = y;
    }
  }
 
  void display() {
    fill(0);
    for ( int i = 0; i < numPoints; i++) {
      float x = xyArray[i][0];
      float y = xyArray[i][1];
      ellipse(x, y, 30, 30);
    }
  }
}

void noiseRect(float x, float y, float w, float h, float W, float c) {
  fill(c);
  float px = random(x+W-1, x+w-W-1);
  float py = random(y+W-1, y+h-W-1);
  rect(px, py, W, W);
}
