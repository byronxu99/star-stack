#include "astro.h"

using std::cout;
using std::endl;

// detect stars in a grayscale image
// we can model stars as points with a gaussian distribution in intensity
Image find_stars_matte(const Image &im, float threshold, float sigma)
{
    // convolve the image with a gaussian, computing how strongly the region
    // around each pixel correlates with the gaussian distribution.
    // this is equivalent to blurring the image!
    Image conv_gaussian = gaussianBlur_separable(im, sigma);
    //(conv_gaussian/conv_gaussian.max()).write("./Output/conv_gaussian.png");

    // get the maximum value in each region the size of the gaussian kernel
    Image max_conv_gaussian = maximum_filter(conv_gaussian, 3.0f*sigma);
    //(max_conv_gaussian/max_conv_gaussian.max()).write("./Output/max_conv_gaussian.png");

    int window_radius = 3.0f * sigma;
    Image out(im.width(), im.height(), 1);
    for(int y=0; y<im.height(); y++) {
        for(int x=0; x<im.width(); x++) {
            out(x, y) = 0.0f;

            // for each local maxima
            if(conv_gaussian(x, y) > 0.0f && conv_gaussian(x, y) == max_conv_gaussian(x, y)) {
                // correlate with gaussian
                float sum = 0;
                float norm = 0;
                for(int dy=-window_radius; dy<=window_radius; dy++) {
                    for(int dx=-window_radius; dx<=window_radius; dx++) {
                        // check bounds
                        if(x+dx < 0 || x+dx >= im.width() || y+dy < 0 || y+dy >= im.height())
                            continue;

                        // compute gaussian
                        float weight = exp(-(dx*dx) / (2.0f *sigma*sigma)) * exp(-(dy*dy) / (2.0f *sigma*sigma));
                        sum += im(x+dx, y+dy) * weight;
                        norm += weight;
                    }
                }
                sum /= norm;

                // check to see that it's under threshold
                // a small point like a star will have a low sum
                // a large white area will have a high sum and won't pass threshold
                if(sum <= threshold) {
                    out(x, y) = 1.0f;
                }
            }
        }
    }

    return out;
}

// detect stars in a color image
Image find_stars_color(const Image &im_color, float mean_factor, float sigma)
{
    // get grayscale image
    Image im = lumiChromi(im_color)[0];
    float mean = im.mean();

    // convolve with gaussian
    Image conv_gaussian = gaussianBlur_separable(im, sigma);

    // get the maximum value in each region the size of the gaussian kernel
    Image max_conv_gaussian = maximum_filter(conv_gaussian, 3.0f*sigma);

    int window_radius = 3.0f * sigma;
    Image out(im.width(), im.height(), 1);
    for(int y=0; y<im.height(); y++) {
        for(int x=0; x<im.width(); x++) {
            out(x, y) = 0.0f;

            // for each local maxima
            if(conv_gaussian(x, y) > 0.0f && conv_gaussian(x, y) == max_conv_gaussian(x, y)) {
                // correlate with gaussian
                float sum = 0;
                float norm = 0;
                for(int dy=-window_radius; dy<=window_radius; dy++) {
                    for(int dx=-window_radius; dx<=window_radius; dx++) {
                        // check bounds
                        if(x+dx < 0 || x+dx >= im.width() || y+dy < 0 || y+dy >= im.height())
                            continue;

                        // compute gaussian
                        float weight = exp(-(dx*dx) / (2.0f *sigma*sigma)) * exp(-(dy*dy) / (2.0f *sigma*sigma));
                        sum += im(x+dx, y+dy) * weight;
                        norm += weight;
                    }
                }
                sum /= norm;

                // check to see that the correlation is strong enough
                // by doing so we remove background noise
                if(sum > mean_factor*mean) {
                    out(x, y) = 1.0f;
                }
            }
        }
    }

    return out;
}

// make an image that shows detected stars
Image visualize_stars(const Image &im, float r, float g, float b)
{
    static const int window_radius = 1;
    Image out(im.width(), im.height(), 3);
    for(int y=0; y<im.height(); y++) {
        for(int x=0; x<im.width(); x++) {
            if(im(x, y) > 0.0f) {
                for(int dy=-window_radius; dy<=window_radius; dy++) {
                    for(int dx=-window_radius; dx<=window_radius; dx++) {
                        if(x+dx < 0 || x+dx >= im.width() || y+dy < 0 || y+dy >= im.height())
                            continue;
                        out(x+dx, y+dy, 0) = r;
                        out(x+dx, y+dy, 1) = g;
                        out(x+dx, y+dy, 2) = b;
                    }
                }
            }
        }
    }
    return out;
}

// remove stars from an alpha matte
Image remove_stars_from_matte(const Image &matte, const Image &stars, float sigma)
{
    int window_radius = 3.0f * sigma;
    Image out(matte.width(), matte.height(), 1);
    for(int y=0; y<matte.height(); y++) {
        for(int x=0; x<matte.width(); x++) {
            out(x, y) = matte(x, y);
        }
    }

    for(int y=0; y<matte.height(); y++) {
        for(int x=0; x<matte.width(); x++) {
            // set the region around each star to zero
            if(stars(x, y)) {
                for(int dy=-window_radius; dy<=window_radius; dy++) {
                    for(int dx=-window_radius; dx<=window_radius; dx++) {
                        // check bounds
                        if(x+dx < 0 || x+dx >= matte.width() || y+dy < 0 || y+dy >= matte.height())
                            continue;

                        // weigh by inverted unnormalized gaussian
                        // equal to 0 at center and approaches 1 at edges
                        float weight = exp(-(dx*dx) / (2.0f *sigma*sigma)) * exp(-(dy*dy) / (2.0f *sigma*sigma));
                        out(x+dx, y+dy) = out(x+dx, y+dy) * (1.0f-weight);
                    }
                }
            }
        }
    }
    return out;
}


// compute harris corners with a star filter
vector<Point> HarrisCorners_stars(const Image &im, float mean_factor, float sigmaS, float k, float sigmaG, float factorSigma, float maxiDiam, float boundarySize, int star_window_radius)
{
    Image R = cornerResponse(im, k, sigmaG, factorSigma);
    Image S = find_stars_color(im, mean_factor, sigmaS);

    // max filter
    Image max_R = maximum_filter(R, maxiDiam);

    vector<Point> out;
    int nearby_star_radius = std::min((int)boundarySize, star_window_radius);
    for(int y = boundarySize; y < R.height()-boundarySize; y++) {
        for(int x = boundarySize; x < R.width()-boundarySize; x++) {
            // find local maxima of corner response
            if(R(x, y, 0) == max_R(x, y, 0)) {
                // only include local maxima near stars
                bool has_nearby_star = false;
                for(int dy=-nearby_star_radius; dy<=nearby_star_radius; dy++) {
                    for(int dx=-nearby_star_radius; dx<=nearby_star_radius; dx++) {
                        if(S(x+dx, y+dy)) {
                            has_nearby_star = true;
                        }
                    }
                }
                if(has_nearby_star) {
                    out.emplace_back(x, y);
                }
            }
        }
    }

    return out;
}

// compute feature correspondences with a star filter
// our assumption here is that stars only move a little between images
// so we check for pairs of features within a squared distance of distance_sq
vector<FeatureCorrespondence> findCorrespondences_stars(const vector<Feature> &listFeatures1, const vector<Feature> &listFeatures2, float distance_sq, float threshold_sq)
{
    vector<FeatureCorrespondence> correspondences;

    for(auto feature1 : listFeatures1) {
        for(auto feature2 : listFeatures2) {
            Point p1 = feature1.point();
            Point p2 = feature2.point();
            float d_sq = pow(p1.x-p2.x, 2) + pow(p1.y-p2.y, 2);
            if(d_sq < distance_sq && l2Features(feature1, feature2) < threshold_sq) {
                correspondences.emplace_back(feature1, feature2);
            }
        }
    }

    return correspondences;
}

// compute a homography that maps one image of stars to another
Matrix compute_star_homography(const Image &im_from, const Image &im_to)
{
    // compute harris corners
    vector<Point> h1 = HarrisCorners_stars(im_from);
    vector<Point> h2 = HarrisCorners_stars(im_to);

    // compute features
    float sigmaBlurDescriptor = 2.0f;
    float radiusDescriptor = 5;
    vector<Feature> f1 = computeFeatures(im_from, h1);//, sigmaBlurDescriptor, radiusDescriptor);
    vector<Feature> f2 = computeFeatures(im_to, h2);  //, sigmaBlurDescriptor, radiusDescriptor);

    //vector<FeatureCorrespondence> corr = findCorrespondences(f1, f2, 1.2f);
    vector<FeatureCorrespondence> corr = findCorrespondences_stars(f1, f2, 1000.0f, 1000.0f);

    // compute homography
    Matrix H = RANSAC(corr, 50000, 1);

    return H;
}

// apply a homography to an image "in-place"
// where the output image has the same size as the input
Image align_with_homography(const Image &im, const Matrix &H)
{
    Image out(im.width(), im.height(), im.channels());
    out = out * 0.0f + -1.0f; // use -1 for "missing" values

    applyHomography(im, H, out, true);
    return out;
}


// stack a series of images
// all images are assumed to have the same size and have 3 channels
Image stack_images(const std::vector<Image> &images, Vec3f (*stack_fn)(const std::vector<Vec3f>&), bool remove_unknown)
{
    int n_images = images.size();
    int width = images[0].width();
    int height = images[0].height();
    int channels = images[0].channels();

    if(channels != 3)
        throw InvalidArgument();

    Image out(width, height, channels);

    // stack each pixel
    for(int y=0; y<height; y++) {
        for(int x=0; x<width; x++) {
            // create a list of pixels
            std::vector<Vec3f> pixels;
            pixels.reserve(n_images);

            for(const Image &im : images) {
                // check if pixel is a known value
                if(remove_unknown && im(x, y, 0) <= 0 && im(x, y, 1) <= 0 && im(x, y, 2) <= 0)
                    continue;

                // add pixel to list of pixels
                pixels.emplace_back(im(x, y, 0), im(x, y, 1), im(x, y, 2));
            }

            // stack the pixels
            if(pixels.size() > 0) {
                Vec3f stack = stack_fn(pixels);
                out(x, y, 0) = stack(0);
                out(x, y, 1) = stack(1);
                out(x, y, 2) = stack(2);
            } else {
                out(x, y, 0) = 0.0f;
                out(x, y, 1) = 0.0f;
                out(x, y, 2) = 0.0f;
            }
        }
    }

    return out;
}

// stacking functions that take a set of pixels and return a single pixel
Vec3f stack_max(const std::vector<Vec3f> &pixels)
{
    Vec3f max = pixels[0];
    for(const Vec3f &p : pixels) {
        max = max.array().max(p.array());
    }
    return max;
}

Vec3f stack_mean(const std::vector<Vec3f> &pixels)
{
    Vec3f sum = Vec3f::Zero();
    for(const Vec3f &p : pixels) {
        sum = sum + p;
    }
    return sum / (float)pixels.size();
}

Vec3f stack_mean_remove_outliers(const std::vector<Vec3f> &pixels)
{
    // we must have at least three samples to compute outliers
    if(pixels.size() < 3)
        return stack_mean(pixels);

    float N = pixels.size();
    Vec3f sum = Vec3f::Zero();
    Vec3f sum_sq = Vec3f::Zero();
    for(const Vec3f &p : pixels) {
        sum += p;
        sum_sq += p.array().square().matrix();
    }
    Vec3f mean = sum / N;
    Vec3f mean_sq = sum_sq / N;
    Vec3f variance = mean_sq - mean.array().square().matrix();
    variance = variance * N / (N-1.0f); // correct for division by N-1 instead of N

    Vec3f sum_without_outliers = Vec3f::Zero();
    for(const Vec3f &p : pixels) {
        Vec3f sq_dist_to_mean = (p-mean).array().square();
        // if all components are less than variance
        if(sq_dist_to_mean.cwiseMin(4.0f*variance) == sq_dist_to_mean) {
            sum_without_outliers += p;
        }
    }
    return sum_without_outliers / (float)pixels.size();
}


// compute pairwise homographies between images
// in the output, out[i] maps the images[i] to the image closer to reference
// out[reference] = identity matrix
std::vector<Matrix> compute_pairwise_homographies(const std::vector<Image> &images, int reference_index)
{
    int n_images = images.size();
    std::vector<Matrix> homographies(n_images, Matrix::Identity(3, 3));

    // compute homographies from reference to start
    for(int i=reference_index-1; i>=0; i--) {
        homographies[i] = compute_star_homography(images[i], images[i+1]);
    }

    // compute homographies from reference to end
    for(int i=reference_index+1; i<n_images; i++) {
        homographies[i] = compute_star_homography(images[i], images[i-1]);
    }

    return homographies;
}

// convert pairwise homographies to absolute homographies from a reference value
// out[i] maps image[i] to image[reference]
// out[reference] = identity matrix
std::vector<Matrix> compute_absolute_homographies(const std::vector<Matrix> &pairwise, int reference_index)
{
    int n_images = pairwise.size();
    std::vector<Matrix> absolute(n_images, Matrix::Identity(3, 3));

    // compute homographies from reference to start
    for(int i=reference_index-1; i>=0; i--) {
        absolute[i] = absolute[i+1] * pairwise[i];
    }

    // compute homographies from reference to end
    for(int i=reference_index+1; i<n_images; i++) {
        absolute[i] = absolute[i-1] * pairwise[i];
    }

    return absolute;
}

// stack images after applying a corresponding homography to each one
Image stack_images_with_homographies(const std::vector<Image> &images, const std::vector<Matrix> &homographies, Vec3f (*stack_fn)(const std::vector<Vec3f>&), bool remove_unknown)
{
    if(images.size() != homographies.size())
        throw InvalidArgument();

    int n_images = images.size();
    std::vector<Image> aligned;

    for(int i=0; i<n_images; i++) {
        aligned.push_back(align_with_homography(images[i], homographies[i]));
    }

    return stack_images(aligned, stack_fn, remove_unknown);
}

// stack images after applying a corresponding homography to each one and with a matte as a guide
// we only stack pixels that are fully background according to the matte
Image stack_images_with_homographies_and_matte(const std::vector<Image> &images, const std::vector<Matrix> &homographies, const Image &matte, Vec3f (*stack_fn)(const std::vector<Vec3f>&), bool remove_unknown)
{
    if(images.size() != homographies.size())
        throw InvalidArgument();

    int width = images[0].width();
    int height = images[0].height();
    int channels = images[0].channels();
    int n_images = images.size();

    // create a mask of only pure background regions
    Image mask(width, height, channels);
    for(int c=0; c<channels; c++)
        for(int y=0; y<height; y++)
            for(int x=0; x<width; x++)
                mask(x, y, c) = (matte(x, y) == 0.0f)? 1.0f : 0.0f;

    // align images
    std::vector<Image> aligned;
    for(int i=0; i<n_images; i++) {
        aligned.push_back(align_with_homography(images[i] * mask, homographies[i]));
    }

    return stack_images(aligned, stack_fn, remove_unknown);
}

