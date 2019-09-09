import blobDetection.*;
BlobDetection theBlobDetection;

int N, M;
PImage img_in;

void settings() { 
  // Make a new instance of a PImage by loading an image file
  img_in = loadImage("../../img/full_silueta.png");
  //img_in.resize(img_in.width/ 2, img_in.height/ 2);
  N = img_in.width; 
  M = img_in.height;
  size(N, M);
}

void setup() {
  noLoop();
}

void draw() {
  //test();
  code();
}

void keyPressed() {
  if (key == 's' || key == 'S') {
    save("BW_blobs2.png");
  }
  if (key == 'q' || key == 'Q') {
    exit();
  }
}

void code() {
  double[][] im = new double[N][M];
  convert2double(im, img_in, red);
  
  theBlobDetection = new BlobDetection(N, M);
  theBlobDetection.setPosDiscrimination(false);
  theBlobDetection.setThreshold(0.38f);
  //theBlobDetection.setThreshold(0.5f);
  theBlobDetection.computeBlobs(img_in.pixels);
  theBlobDetection.computeTriangles();

  image(img_in, 0, 0, N, M);
  drawBlobsAndEdges(true, true, true, N, M);
}

void drawBlobsAndEdges(boolean drawBlobs, boolean drawEdges, boolean drawTriangles, int rows, int cols) {
  noFill();
  Blob b;
  EdgeVertex eA, eB, eC;
  BlobTriangle tr;
  for (int n=0 ; n < theBlobDetection.getBlobNb() ; n++) {
    b = theBlobDetection.getBlob(n);
    if (b != null) {
      // Edges
      if (drawEdges) {
        strokeWeight(2);
        stroke(146, 212, 241);
        for (int m=0; m < b.getEdgeNb(); m++) {
          eA = b.getEdgeVertexA(m);
          eB = b.getEdgeVertexB(m);
          if (eA !=null && eB !=null)
            line(eA.x*rows, eA.y*cols, eB.x*rows, eB.y*cols);
        }
      }

      // Blobs
      if (drawBlobs) {
        strokeWeight(1);
        stroke(255, 0, 0);
        rect(b.xMin* rows, b.yMin* cols, b.w* rows, b.h* cols);
      }
      
      if (drawTriangles) {
        println(b.getTriangleNb());
        for (int t=0 ; t < b.getTriangleNb() ; t++) {
          tr = b.getTriangle(t);
          if(tr != null) {
            eA = b.getTriangleVertexA(tr);
            eB = b.getTriangleVertexB(tr);
            eC = b.getTriangleVertexC(tr);
            if (eA !=null && eB !=null && eC != null) {
               beginShape();
               vertex(eA.x* rows, eA.y* cols);
               vertex(eB.x* rows, eB.y* cols);
               vertex(eC.x* rows, eC.y* cols);
               endShape(CLOSE);
            }
          }
        }
      }
    }
  }
}

void test() {
  int M = 8, N = 5;
  Matrix B = Matrix.random(5, 3);
  Matrix A = Matrix.random(M, N).times(B).times(B.transpose());
  println("A = ");
  A.print(9, 6);
  println("A = U S V^T");
  println();
  SVD s = new SVD(A);
  print("U = ");
  Matrix U = s.getU();
  U.print(9, 6);
  print("Sigma = ");
  Matrix S = s.getS();
  S.print(9, 6);
  print("V = ");
  Matrix V = s.getV();
  V.print(9, 6);
  println("rank = " + s.rank());
  println("condition number = " + s.cond());
  println("2-norm = " + s.norm2());
}
