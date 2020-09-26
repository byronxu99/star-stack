#ifndef __CFM__H
#define __CFM__H

#include "matrix.h"
#include "Image.h"
#include "basicImageManipulation.h"
#include "filtering.h"

// convert a WxHxC image to a (WH)xC 2D matrix
Matrix img2eigen(const Image &im);

// return a (window area)x(channels) matrix around pixel i,j
Matrix get_window(const Matrix &I, int i, int j, int radius, int image_width);

// compute average for each channel, returning a 1xC matrix
Matrix mu_window(const Matrix &W);

// subtract mu from each row of W
Matrix window_minus_mu(const Matrix &W, const Matrix &mu);

// compute CxC covariance matrix
Matrix sigma_window(const Matrix &W, const Matrix &mu);

// compute the matting laplacian matrix given image and mask of known pixels
SparseMatrix matting_laplacian(const Image &im, const Image &known, float epsilon = 1e-6, int radius = 1);

// compute matte given image and trimap
Image cfm_with_trimap(const Image &im, const Image &trimap, float lambda = 100.f, float epsilon = 1e-6, int radius = 1);

// iterative closed-form matting by scaling image
Image iterative_cfm(const Image &im, const Image &trimap, float threshold = 0.01f, float decimation = 0.5f, int max_size = 2000, float lambda = 100.f, float epsilon = 1e-6, int radius = 1);

// convert a set of scribbles to a trimap
Image scribble2trimap(const Image &im, const Image &scribble);

// generate a trimap guess for landscape photos
Image sky_ground_trimap(const Image &im, float frac_sky = 0.25f, float frac_ground = 0.1f);

// threshold a matte
Image threshold_matte(const Image &matte, float lo = 0.1f, float hi = 0.9f);

// compute foreground and background given matte
void compute_fg_bg(const Image &im, const Image &matte, Image &fg_out, Image &bg_out, float fg_threshold = 1.0f, float bg_threshold = 0.0f);

// blend foreground and background given a matte
Image blend_fg_bg(const Image &matte, const Image &fg, const Image &bg);

#endif
