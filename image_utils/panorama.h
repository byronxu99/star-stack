/* --------------------------------------------------------------------------
 * File:    panorama.h
 * Created: 2015-10-29
 * --------------------------------------------------------------------------
 *
 *
 *
 * ------------------------------------------------------------------------*/


#ifndef __panorama__h
#define __panorama__h

#include "Image.h"
#include "matrix.h"
#include "filtering.h"
#include "homography.h"
#include "basicImageManipulation.h"
#include <iostream>
#include <algorithm>    // std::max, random_shuffle
#include <cmath>

using namespace std;


// panorama.h
// Assignment 7

// ---------------------------------------------------------
// Point and Feature classes, already implemented for Pset07
class Point {
public:
    int x, y; // point info
    Point(int xp, int yp); // Constructors
    Point();              // Constructors
    void print() const; // helpful printing
    Vec3f toHomogenousCoords() const;
};


class Feature {
public:
    Feature(const Point &ptp, const Image &descp); // Constructor
    const Point& point() const; // get point (this calls point cc)
    const Image& desc() const; // get Image (this maybe calls image cc)
    void print() const; // pretty printing
private:
    Point pt;
    Image dsc; // 9x9 descriptor
};


class FeatureCorrespondence {
public:
    FeatureCorrespondence(const Feature &f0p, const Feature &f1p); // Constructors
    vector<Feature> features() const; //vector of features [f0, f1]
    const Feature& feature(int i) const;    // return f0 if i == 0, f1 if i ==1
    void print() const; // helpful printing
    CorrespondencePair toCorrespondencePair() const;
private:
    Feature f0, f1; // corresp info
};
// ---------------------------------------------------------



// Pset07: Harris Corners
Image computeTensor(const Image &im, float sigmaG = 1, float factorSigma = 4);
Image cornerResponse(const Image &im,
                     float k = 0.15f,
                     float sigmaG = 1,
                     float factorSigma = 4);
vector<Point> HarrisCorners(const Image &im,
                            float k = 0.15f,
                            float sigmaG = 1,
                            float factorSigma = 4,
                            float maxiDiam = 7,
                            float boundarySize = 5);

// Pset07: FeatureCorrespondences
Image descriptor(const Image &blurredIm, const Point &p, float radiusDescriptor=4);
vector<Feature> computeFeatures(const Image &im,
                                const vector<Point> &cornersL,
                                float sigmaBlurDescriptor = 0.5,
                                float radiusDescriptor = 4);
float l2Features(const Feature &f1, const Feature &f2);
vector<FeatureCorrespondence> findCorrespondences(
    const vector<Feature> &listFeatures1,
    const vector<Feature> &listFeatures2,
    float threshold = 1.7f);

// Pset07: RANSAC
vector<CorrespondencePair> getListOfPairs(
    const vector<FeatureCorrespondence> &listOfCorrespondences);
vector<bool> inliers(const Matrix &H,
                     const vector<FeatureCorrespondence> &listOfCorrespondences,
                     float epsilon = 4);
Matrix RANSAC(const vector<FeatureCorrespondence> &listOfCorrespondences,
              int Niter = 200,
              float epsilon=4);
vector<FeatureCorrespondence> sampleFeatureCorrespondences(
    vector<FeatureCorrespondence> listOfCorrespondences);

// PSet07: Final stitching
Image autostitch(const Image &im1,
                 const Image &im2,
                 float blurDescriptor = 0.5,
                 float radiusDescriptor = 4);



// ---------------------------------------------------------

// potentially useful function, optional to implement
Image getBlurredLumi(const Image &im, float sigmaG);
int countBoolVec(const vector<bool> &ins);

// ---------------------------------------------------------

/***********************************************************************
 * Helper Functions, already implemented in Pset07 *
 ***********************************************************************/

void drawLine(const Point &p1,
              const Point &p2,
              Image &im,
              const vector<float> &color = vector<float>{1.0, 1.0, 1.0});

Image visualizeCorners(const Image &im,
                       const vector<Point> &pts,
                       int rad = 2,
                       const vector<float> &color = vector<float>{1.0, 1.0, 1.0});

Image visualizeFeatures(const Image &im,
                        const vector<Feature> &LF,
                        float radiusDescriptor=4);

Image visualizePairs(const Image &im1,
                     const Image &im2,
                     const vector<FeatureCorrespondence> &corr);

Image visualizePairsWithInliers(const Image &im1,
                                const Image &im2,
                                const vector<FeatureCorrespondence> &corr,
                                const vector<bool> &ins);

vector<Image> visualizeReprojection(const Image &im1,
                                    const Image &im2,
                                    const Matrix &H,
                                    const vector<FeatureCorrespondence> &corr,
                                    const vector<bool> &ins);

#endif
