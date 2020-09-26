#include "homography.h"
#include "matrix.h"

using namespace std;


void applyHomography(const Image &source, const Matrix &H, Image &out, bool bilinear) {
    // // --------- HANDOUT  PS06 ------------------------------
    // Transform image source using the homography H, and composite in onto out.
    // if bilinear == true, using bilinear interpolation. Use nearest neighbor
    // otherwise.
    Matrix H1 = H.inverse();
    for(int x = 0; x < out.width(); x++) {
        for(int y = 0; y < out.height(); y++) {
            // transform
            // x,y = output
            // xx,yy = input
            Vec3f dest(x, y, 1);
            Vec3f src = H1 * dest;
            float xx = src(0) / src(2);
            float yy = src(1) / src(2);
            int xx_int = round(xx);
            int yy_int = round(yy);

            // check bounds
            if(xx_int < 0 || xx_int >= source.width()
            || yy_int < 0 || yy_int >= source.height())
                continue;

            // modify out
            if(bilinear) {
                for(int c = 0; c < out.channels(); c++) {
                    out(x, y, c) = interpolateLin(source, xx, yy, c, true);
                }
            } else {
                for(int c = 0; c < out.channels(); c++) {
                    out(x, y, c) = source(xx_int, yy_int, c);
                }
            }
        }
    }
}




Matrix computeHomography(const CorrespondencePair correspondences[4]) {
    // --------- HANDOUT  PS06 ------------------------------
    // Compute a homography from 4 point correspondences.
    // x, y
    float x0 = correspondences[0].point1(0) / correspondences[0].point1(2);
    float y0 = correspondences[0].point1(1) / correspondences[0].point1(2);
    float x1 = correspondences[1].point1(0) / correspondences[1].point1(2);
    float y1 = correspondences[1].point1(1) / correspondences[1].point1(2);
    float x2 = correspondences[2].point1(0) / correspondences[2].point1(2);
    float y2 = correspondences[2].point1(1) / correspondences[2].point1(2);
    float x3 = correspondences[3].point1(0) / correspondences[3].point1(2);
    float y3 = correspondences[3].point1(1) / correspondences[3].point1(2);

    // x', y'
    float xp0 = correspondences[0].point2(0) / correspondences[0].point2(2);
    float yp0 = correspondences[0].point2(1) / correspondences[0].point2(2);
    float xp1 = correspondences[1].point2(0) / correspondences[1].point2(2);
    float yp1 = correspondences[1].point2(1) / correspondences[1].point2(2);
    float xp2 = correspondences[2].point2(0) / correspondences[2].point2(2);
    float yp2 = correspondences[2].point2(1) / correspondences[2].point2(2);
    float xp3 = correspondences[3].point2(0) / correspondences[3].point2(2);
    float yp3 = correspondences[3].point2(1) / correspondences[3].point2(2);

    // Ax should equal (w0'x0', w0'y0', ...)T
    Matrix A(8, 8);
    A <<
        x0, y0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, x0, y0, 1, 0, 0,
        x1, y1, 1, 0, 0, 0, 0, 0,
        0, 0, 0, x1, y1, 1, 0, 0,
        x2, y2, 1, 0, 0, 0, 0, 0,
        0, 0, 0, x2, y2, 1, 0, 0,
        x3, y3, 1, 0, 0, 0, 0, 0,
        0, 0, 0, x3, y3, 1, 0, 0;

    // B = (x0', y0', x1', y1', ...)T
    Matrix B(8, 1);
    B << xp0, yp0, xp1, yp1, xp2, yp2, xp3, yp3;

    // AAx = (w0'-1, w0'-1, w1'-1, w1'-1, ...)T
    Matrix AA(8, 8);
    AA <<
        0, 0, 0, 0, 0, 0, x0, y0,
        0, 0, 0, 0, 0, 0, x0, y0,
        0, 0, 0, 0, 0, 0, x1, y1,
        0, 0, 0, 0, 0, 0, x1, y1,
        0, 0, 0, 0, 0, 0, x2, y2,
        0, 0, 0, 0, 0, 0, x2, y2,
        0, 0, 0, 0, 0, 0, x3, y3,
        0, 0, 0, 0, 0, 0, x3, y3;

    // AAx + (11111111)T = (w0', w0', w1', w1', ...)T
    // diagonal(B) * (AAx + (11111111)T) = (w0'x0', w0'y0', ...)T = Ax
    // therefore (A - (diagonal(B) * AA))x = diagonal(B) * (11111111)T = B
    // compute correct A for solving
    A = A - (B.asDiagonal() * AA);

    // Ax = B, solve for x
    Matrix x = A.fullPivLu().solve(B);
    Matrix out(3, 3);
    out <<
        x(0), x(1), x(2),
        x(3), x(4), x(5),
        x(6), x(7), 1;
    return out;
}

Matrix computeHomographySVD(const CorrespondencePair correspondences[4]) {
    // x, y
    float x0 = correspondences[0].point1(0) / correspondences[0].point1(2);
    float y0 = correspondences[0].point1(1) / correspondences[0].point1(2);
    float x1 = correspondences[1].point1(0) / correspondences[1].point1(2);
    float y1 = correspondences[1].point1(1) / correspondences[1].point1(2);
    float x2 = correspondences[2].point1(0) / correspondences[2].point1(2);
    float y2 = correspondences[2].point1(1) / correspondences[2].point1(2);
    float x3 = correspondences[3].point1(0) / correspondences[3].point1(2);
    float y3 = correspondences[3].point1(1) / correspondences[3].point1(2);

    // x', y'
    float xp0 = correspondences[0].point2(0) / correspondences[0].point2(2);
    float yp0 = correspondences[0].point2(1) / correspondences[0].point2(2);
    float xp1 = correspondences[1].point2(0) / correspondences[1].point2(2);
    float yp1 = correspondences[1].point2(1) / correspondences[1].point2(2);
    float xp2 = correspondences[2].point2(0) / correspondences[2].point2(2);
    float yp2 = correspondences[2].point2(1) / correspondences[2].point2(2);
    float xp3 = correspondences[3].point2(0) / correspondences[3].point2(2);
    float yp3 = correspondences[3].point2(1) / correspondences[3].point2(2);

    Matrix A(8, 9);
    A <<
        -x0, -y0, -1, 0, 0, 0, xp0*x0, xp0*y0, xp0,
        0, 0, 0, -x0, -y0, -1, yp0*x0, yp0*y0, yp0,
        -x1, -y1, -1, 0, 0, 0, xp1*x1, xp1*y1, xp1,
        0, 0, 0, -x1, -y1, -1, yp1*x1, yp1*y1, yp1,
        -x2, -y2, -1, 0, 0, 0, xp2*x2, xp2*y2, xp2,
        0, 0, 0, -x2, -y2, -1, yp2*x2, yp2*y2, yp2,
        -x3, -y3, -1, 0, 0, 0, xp3*x3, xp3*y3, xp3,
        0, 0, 0, -x3, -y3, -1, yp3*x3, yp3*y3, yp3;

    // Ax = 0, solve for x
    // solution is the last column of V matrix in the SVD
    Eigen::JacobiSVD<Matrix> svd(A, Eigen::ComputeFullV);
    //cout << svd.singularValues() << endl;
    //cout << svd.matrixV() << endl;
    Vector x = svd.matrixV().col(svd.matrixV().cols()-1);

    Matrix out(3, 3);
    out <<
        x(0), x(1), x(2),
        x(3), x(4), x(5),
        x(6), x(7), x(8);
    return out;
}


BoundingBox computeTransformedBBox(int imwidth, int imheight, Matrix H) {
    // --------- HANDOUT  PS06 ------------------------------
    // Predict the bounding boxes that encompasses all the transformed
    // coordinates for pixels frow and Image with size (imwidth, imheight)

    // corners of image
    Vec3f p0(0, 0, 1);
    Vec3f p1(imwidth-1, 0, 1);
    Vec3f p2(0, imheight-1, 1);
    Vec3f p3(imwidth-1, imheight-1, 1);

    // transform corners
    Vec3f q0 = H * p0;
    Vec3f q1 = H * p1;
    Vec3f q2 = H * p2;
    Vec3f q3 = H * p3;

    // convert homogenous to regular coordinates
    q0 = q0 / q0(2);
    q1 = q1 / q1(2);
    q2 = q2 / q2(2);
    q3 = q3 / q3(2);

    float min_x = min(min(q0(0), q1(0)), min(q2(0), q3(0)));
    float min_y = min(min(q0(1), q1(1)), min(q2(1), q3(1)));
    float max_x = max(max(q0(0), q1(0)), max(q2(0), q3(0)));
    float max_y = max(max(q0(1), q1(1)), max(q2(1), q3(1)));

    return BoundingBox(round(min_x), round(max_x), round(min_y), round(max_y));
}


BoundingBox bboxUnion(BoundingBox B1, BoundingBox B2) {
    // --------- HANDOUT  PS06 ------------------------------
    // Compute the bounding box that tightly bounds the union of B1 an B2.
    return BoundingBox(min(B1.x1, B2.x1), max(B1.x2, B2.x2), min(B1.y1, B2.y1), max(B1.y2, B2.y2));

}


Matrix makeTranslation(BoundingBox B) {
    // --------- HANDOUT  PS06 ------------------------------
    // Compute a translation matrix (as a homography matrix) that translates the
    // top-left corner of B to (0,0).
    Matrix tr(3, 3);
    tr <<
        1, 0, -B.x1,
        0, 1, -B.y1,
        0, 0, 1;
    return tr;
}


Image stitch(const Image &im1, const Image &im2, const CorrespondencePair correspondences[4]) {
    // --------- HANDOUT  PS06 ------------------------------
    // Transform im1 to align with im2 according to the set of correspondences.
    // make sure the union of the bounding boxes for im2 and transformed_im1 is
    // translated properly (use makeTranslation)

    Matrix H = computeHomographySVD(correspondences);
    BoundingBox im2_bbox    = BoundingBox(0, im2.width()-1, 0, im2.height()-1);
    BoundingBox im1_tr_bbox = computeTransformedBBox(im1.width(), im1.height(), H);
    BoundingBox bbox_union  = bboxUnion(im1_tr_bbox, im2_bbox);
    Matrix tr = makeTranslation(bbox_union);

    Image out(bbox_union.x2 - bbox_union.x1 + 1, bbox_union.y2 - bbox_union.y1 + 1, im1.channels());
    applyHomography(im1, tr*H, out, true);
    applyHomography(im2, tr,   out, true);
    return out;
}

// debug-useful
Image drawBoundingBox(const Image &im, BoundingBox bbox) {
    // // --------- HANDOUT  PS06 ------------------------------
    /*
      ________________________________________
     / Draw me a bounding box!                \
     |                                        |
     | "I jumped to my                        |
     | feet, completely thunderstruck. I      |
     | blinked my eyes hard. I looked         |
     | carefully all around me. And I saw a   |
     | most extraordinary small person, who   |
     | stood there examining me with great    |
     | seriousness."                          |
     \              Antoine de Saint-Exupery  /
      ----------------------------------------
             \   ^__^
              \  (oo)\_______
                 (__)\       )\/\
                     ||----w |
                     ||     ||
    */

    Image out = im;

    for(int x = bbox.x1; x <= bbox.x2; x++) {
        for(int c = 0; c < im.channels(); c++) {
            out(x, bbox.y1, c) = 1;
            out(x, bbox.y2, c) = 1;
        }
    }
    for(int y = bbox.y1; y <= bbox.y2; y++) {
        for(int c = 0; c < im.channels(); c++) {
            out(bbox.x1, y, c) = 1;
            out(bbox.x2, y, c) = 1;
        }
    }

    return out;
}

void applyHomographyFast(const Image &source, const Matrix &H, Image &out, bool bilinear) {
    // // --------- HANDOUT  PS06 ------------------------------
    // Same as apply but change only the pixels of out that are within the
    // predicted bounding box (when H maps source to its new position).
    BoundingBox bbox = computeTransformedBBox(source.width(), source.height(), H);

    Matrix H1 = H.inverse();
    for(int x = bbox.x1; x <= bbox.x2; x++) {
        for(int y = bbox.y1; y <= bbox.y2; y++) {
            // transform
            // x,y = output
            // xx,yy = input
            Vec3f dest(x, y, 1);
            Vec3f src = H1 * dest;
            float xx = src(0) / src(2);
            float yy = src(1) / src(2);
            int xx_int = round(xx);
            int yy_int = round(yy);

            // check bounds
            if(xx_int < 0 || xx_int >= source.width()
            || yy_int < 0 || yy_int >= source.height())
                continue;

            // modify out
            if(bilinear) {
                for(int c = 0; c < out.channels(); c++) {
                    out(x, y, c) = interpolateLin(source, xx, yy, c, true);
                }
            } else {
                for(int c = 0; c < out.channels(); c++) {
                    out(x, y, c) = source(xx_int, yy_int, c);
                }
            }
        }
    }
}
