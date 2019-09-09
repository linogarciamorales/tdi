

void setup() {
  //size(256, 120);
  size(256, 256);
  noStroke();
  smooth();
  noLoop();
}

void draw() {
  background(0);
  //code1();
  //code2();
  //code3();
  //code4();
  code5();
}

void keyPressed() {
  if (key == 'q' || key == 'Q') {
    exit();
  }
  if (key == 's' || key == 'S') {
    save("../../img/patt5.png");    
  }
}

void code1() {
  float ofs = 3.4;
  float W = 100;
  ellipse (width/ofs, height/2, W, W);
  ellipse (width*(1-1/ofs), height/2, W, W);
  
  float w = 30;
  float r = w/1.1;
  Blackhole b1 = new Blackhole(width/ofs, height/2, r, 4);
  b1.divide();
  b1.display();
  
  Blackhole b2 = new Blackhole(width*(1-1/ofs), height/2, r, 4);
  b2.divide();
  b2.display();
}

void code2() {
  fill(255);
  ellipse(width/2, height/2, height-75, height-75);
  fill(0);
  rect(0, 0, width, height/2);
  fill(255);
  ellipse(width/2.6, height/2, height/2, height/2);
}

void code3() {
  fill(255);
  float w1 = 30;
  rect(w1, w1, width-2*w1, height-2*w1);
  fill(0);
  float w2 = 20;
  rect(w1+w2, w1+w2, width-2*w1-2*w2, height-2*w1-2*w2);
}

void code4() {
  fill(255);
  float w1 = 30;
  rect(w1, w1, width-2*w1, height-2*w1);
  fill(0);
  float w2 = 20;
  rect(w1+w2, w1+w2, width-2*w1-2*w2, height-2*w1-2*w2);
  int N = 3;
  float w3 = 5;
  for(int k= 0; k < N; k++) {
    noiseRect(w1, w1, w2, height-2*w1, w3, 0);
    noiseRect(width-w1-w2, w1, w2, height-2*w1, w3, 0);
    noiseRect(w1, w1, width-2*w1, w2, w3, 0);
    noiseRect(w1, height-w1-w2, width-2*w1, w2, w3, 0);
  }
}

void code5() {
  fill(255);
  float w1 = 30;
  rect(w1, w1, width-2*w1, height-2*w1);
  fill(0);
  float w2 = 20;
  rect(w1+w2, w1+w2, width-2*w1-2*w2, height-2*w1-2*w2);
  int N = 3;
  float w3 = 5;
  for(int k= 0; k < N; k++) {
    noiseRect(w1+w2, w1+w2, width-2*(w1+w2), height-2*(w1+w2), w3, 255);
    noiseRect(0, 0, w1, height, w3, 255);
    noiseRect(width-w1, 0, w1, height, w3, 255);
    noiseRect(0, 0, width, w1, w3, 255);
    noiseRect(0, height-w1, width, w1, w3, 255);
  }
}
