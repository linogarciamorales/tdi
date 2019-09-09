// author: thomas diewald
// date: 21.06.2013
// 
// PCA - Principal Component Analysis
//
// in short:
// the eigenvectors (sorted by their eigenvalue) of the covarianz-matrix, which
// is built of the data-vectors(R,G,B), make the tranformation matrix.
// so the first column (eigenvector with the biggest eigenvalue) is the
// most significant component.
//
// - 2D -
// '1' ... source image
// '2' ... destination image
// - 3D -
// '3' ... original RGB values 
// '4' ... transformed RGB values (eigenvector matrix)
// '5' ... back-transformed RGB values (eigenvector matrix)
//
// 'q' | 'w' ... change scale of 1st component
// 'a' | 's' ... change scale of 2nd component
// 'y' | 'x' ... change scale of 3rd component
// 
// to see the affect of scaling choose display '2', '4' or '5'
//

import peasy.PeasyCam;

  
  PCA pca;

  PeasyCam cam;
  
  PImage img_src;
  PImage img_cmp;
  
  DwVector[] data_src; // source data
  DwVector[] data_trn; // transformed data
  DwVector[] data_cmp; // compressed data, after backtransformation and dimension reduction
  
  int DISPLAY = 2;
  
  // scaling of the data-vectors
  float SX = 1;
  float SY = 1;
  float SZ = 0;

  boolean online = true;
  
  public void setup() {
    size(512, 512, P3D);

    cam = new PeasyCam(this, 300);
    
    
    // loading src image / creating dst image
    img_src = loadImage("../../img/lena.jpg");
    img_cmp = createImage(img_src.width, img_src.height, ARGB);
    
    
    long timer;
    
 
    // creating data vectors from rgb pixels
    timer = System.currentTimeMillis();
    System.out.println(">> creating rgb vectors");

    data_src = new DwVector[img_src.pixels.length];
    img_src.loadPixels();
    for(int i = 0; i < img_src.pixels.length; i++){
      int argb = img_src.pixels[i];
      int r = (argb>>16)&0xFF;
      int g = (argb>> 8)&0xFF;
      int b = (argb    )&0xFF;
      data_src[i] = new DwVector(r,g,b);
    }
    timer = System.currentTimeMillis()-timer;
    System.out.println("   time: "+timer+" ms");
    System.out.println("   num vectors: "+data_src.length);
    System.out.println("   len vec: "+data_src[0].v.length);
    
    
    // PCA
    timer = System.currentTimeMillis();
    System.out.println(">> creating PCA (Principal Axis Analysis)");
    pca = new PCA( data_src ).compute();
    timer = System.currentTimeMillis()-timer;
    System.out.println("   time: "+timer+" ms");
    
    // reducing dimensions and transform data
    reduceData(1,1,0);
    
    // create image of (reduced) transformed data
    updateDstImage();
    
    
    textFont(createFont("Calibri", 12));
    textMode(SCREEN);
   
  }
  
  

  
  
  void reduceData(float sx, float sy, float sz){
    long timer;
    
    // transform data
    timer = System.currentTimeMillis();
    System.out.println(">> transforming data");
    data_trn = pca.transformData(pca.data, false);
    timer = System.currentTimeMillis()-timer;
    System.out.println("   time: "+timer+" ms");
    
    // reducing dimensions, or apply anything else  here
    for(int i = 0; i < data_trn.length; i++){
      data_trn[i].v[0] *= sx;
      data_trn[i].v[1] *= sy;
      data_trn[i].v[2] *= sz;
    }
    
    
    // transform data back
    data_cmp = pca.transformData(data_trn, !false);
  }
  
  
  
  void updateDstImage(){
    for(int i = 0; i < img_cmp.pixels.length; i++){
      int r = (int)(data_cmp[i].v[0] + pca.mean.v[0]);
      int g = (int)(data_cmp[i].v[1] + pca.mean.v[1]);
      int b = (int)(data_cmp[i].v[2] + pca.mean.v[2]);
      
      if(r < 0 ) r = 0; else if(r > 255 ) r = 255;
      if(g < 0 ) g = 0; else if(g > 255 ) g = 255;
      if(b < 0 ) b = 0; else if(b > 255 ) b = 255;

      img_cmp.pixels[i] = 0xFF000000 | r<<16 | g << 8 | b;
    }
    img_cmp.updatePixels();
  }
  


  public void draw(){
    background(255);
    
    
    if( DISPLAY >= 3 ){
     
      float hs = 255*0.5f;
      pushMatrix();
      {
        translate(hs-pca.mean.v[0],
                  hs-pca.mean.v[1],
                  hs-pca.mean.v[2] );
        stroke(0);
        strokeWeight(1);
        noFill();
        box(255);
      }
      popMatrix();
      strokeWeight(1);
      gizmo(hs);
      
      if( DISPLAY == 3 ) drawDataVectors(pca, data_src, data_src);
      if( DISPLAY == 4 ) drawDataVectors(pca, data_trn, data_src);
      if( DISPLAY == 5 ) drawDataVectors(pca, data_cmp, data_src);
      
      strokeWeight(3);
      drawEigenVector(pca.evec[0], 360, 0xFFFF0000);
      drawEigenVector(pca.evec[1], 240, 0xFF00FF00);
      drawEigenVector(pca.evec[2], 120, 0xFF0000FF);
      
    } else {
      cam.beginHUD();
      if( DISPLAY == 1 ) image(img_src, 0, 0);
      if( DISPLAY == 2 ) image(img_cmp, 0, 0);
      cam.endHUD();
    }
    
    cam.beginHUD();
    
    
    String txt_display = "";
    if( DISPLAY == 1 ) txt_display = "[1] source image";
    if( DISPLAY == 2 ) txt_display = "[2] destination image";
    if( DISPLAY == 3 ) txt_display = "[3] original RGB";
    if( DISPLAY == 4 ) txt_display = "[4] transformed RGB";
    if( DISPLAY == 5 ) txt_display = "[5] backtransformed RGB";
    
    String txt = txt_display;
    fill(0);
    text(txt, 20, 20);
    cam.endHUD();
    
  }
  
  
  void drawDataVectors(PCA pca, DwVector[] pos, DwVector[] col){
    strokeWeight(1);
    for(int i = 0; i < pos.length; i++){
      DwVector data_pos = pos[i];
      DwVector data_col = col[i]; 
      float x = data_pos.v[0];
      float y = data_pos.v[1];
      float z = data_pos.v[2];

      float r = data_col.v[0]+ pca.mean.v[0];
      float g = data_col.v[1]+ pca.mean.v[1];
      float b = data_col.v[2]+ pca.mean.v[2];

      stroke(r,g,b);
      point(x,y,z);
    }
  }
  
  void drawEigenVector(DwEigenVector evec, float len, int col){
    float x = (float) evec.evec[0] * len;
    float y = (float) evec.evec[1] * len;
    float z = (float) evec.evec[2] * len;
    stroke(col);
    line(0,0,0, x,y,z);
  }
  
  
  void gizmo(float s){
    stroke(155,0,0); line(0,0,0, s,0,0);
    stroke(0,155,0); line(0,0,0, 0,s,0);
    stroke(0,0,155); line(0,0,0, 0,0,s);
  }
  
  
  public void keyPressed(){
    if( online && key == ESC ) key = 0;
    
    if( key == 'q' || key == 'w' ||
        key == 'a' || key == 's' ||
        key == 'y' || key == 'x')
    {
      float s = .1f;
      if( key == 'q') SX -= s;
      if( key == 'w') SX += s;
      if( key == 'a') SY -= s;
      if( key == 's') SY += s;
      if( key == 'y') SZ -= s;
      if( key == 'x') SZ += s;
      
      reduceData(SX, SY, SZ);
      updateDstImage();
    }
  }
  
  public void keyReleased(){
    
    if( online && key == ESC ) key = 0;
    
    if( key == '1' ) DISPLAY = 1;
    if( key == '2' ) DISPLAY = 2;
    if( key == '3' ) DISPLAY = 3;
    if( key == '4' ) DISPLAY = 4;
    if( key == '5' ) DISPLAY = 5;
    
    if( !online && key =='s'){
      save("data/lena_PCA_transformation.png");
    }
    
  }
