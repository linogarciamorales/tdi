float minimo = 0;
float maximo = 255;
int Q = 10;
float escala[] = new float[Q];

int w = 100, h = 100;

void settings() {
  size(w*Q, h);
}

void setup() {
  noStroke();
  noLoop();

  setLinealScale();
  //setWeberFechnerScale(5);
}

void draw() {
  background(0);
  for(int k = 0; k < Q; k++) {
    fill(escala[k]);
    rect(k*w, 0, w, h);
  }
  save("escala.png");
}

void setLinealScale() {
  float intervalo = (maximo - minimo)/ (Q-1);
  println(intervalo);

  for(int k = 0; k < Q; k++) {
    escala[k] = intervalo*k;
    print(escala[k]+" ");
  }
}

void setWeberFechnerScale(float eps) {
  float intervalo = maximo - minimo;
  float K = 1.0/ (Q-1);
  float r = exp(K * log(intervalo/ eps));

  escala[0] = minimo;
  escala[1] = eps;
  println();
  for(int k = 1; k < Q; k++) {
    escala[k] = eps* pow(r, k-1);
    print(escala[k]+" ");
  }
  println(eps* pow(r, Q-1));
}
