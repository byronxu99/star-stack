#ifndef __ASTRO__H
#define __ASTRO__H

#include "matrix.h"
#include "Image.h"
#include "basicImageManipulation.h"
#include "filtering.h"
#include "homography.h"
#include "panorama.h"

// detect stars in a grayscale image
Image find_stars_matte(const Image &im, float threshold = 0.3f, float sigma = 2.0f);

// detect stars in a color image
Image find_stars_color(const Image &im_color, float mean_factor = 2.0f, float sigma = 1.0f);

// make an image that shows detected stars
Image visualize_stars(const Image &im, float r = 1.0f, float g = 0.0f, float b = 0.0f);

// remove stars from an alpha matte
Image remove_stars_from_matte(const Image &matte, const Image &stars, float sigma = 2.0f);

// compute harris corners with a star filter
vector<Point> HarrisCorners_stars(const Image &im,
                            float mean_factor = 2.0f,
                            float sigmaS = 1.0f,
                            float k = 0.15f,
                            float sigmaG = 2.0f,
                            float factorSigma = 4,
                            float maxiDiam = 7,
                            float boundarySize = 5,
                            int   star_window_radius = 2);

// compute feature correspondences with a star filter
vector<FeatureCorrespondence> findCorrespondences_stars(
        const vector<Feature> &listFeatures1,
        const vector<Feature> &listFeatures2,
        float distance_sq, float threshold_sq);

// compute a homography that maps one image of stars to another
Matrix compute_star_homography(const Image &im_from, const Image &im_to);

// apply a homography to an image "in-place"
// where the output image has the same size as the input
Image align_with_homography(const Image &im, const Matrix &H);

// stack a series of images
Image stack_images(const std::vector<Image> &images, Vec3f (*stack_fn)(const std::vector<Vec3f>&), bool remove_unknown = true);

// stacking functions that take a set of pixels and return a single pixel
Vec3f stack_max(const std::vector<Vec3f> &pixels);
Vec3f stack_mean(const std::vector<Vec3f> &pixels);
Vec3f stack_mean_remove_outliers(const std::vector<Vec3f> &pixels);

// compute pairwise homographies between images
std::vector<Matrix> compute_pairwise_homographies(const std::vector<Image> &images, int reference_index = 0);

// convert pairwise homographies to absolute homographies from a reference value
std::vector<Matrix> compute_absolute_homographies(const std::vector<Matrix> &pairwise, int reference_index = 0);

// stack images after applying a corresponding homography to each one
Image stack_images_with_homographies(const std::vector<Image> &images, const std::vector<Matrix> &homographies, Vec3f (*stack_fn)(const std::vector<Vec3f>&), bool remove_unknown = true);

// stack images after applying a corresponding homography to each one and with a matte as a guide
Image stack_images_with_homographies_and_matte(const std::vector<Image> &images, const std::vector<Matrix> &homographies, const Image &matte, Vec3f (*stack_fn)(const std::vector<Vec3f>&), bool remove_unknown = true);

#endif
