#include "matrix.h"
#include "Image.h"
#include "homography.h"
#include "panorama.h"
#include "timing.h"
#include "cfm.h"
#include "astro.h"
#include "tests.h"

using std::cout;
using std::endl;

Image make_simple_test_img()
{
    Image out(4, 4, 3);
    for(int x=0; x<out.width(); x++) {
        for(int y=0; y<out.height(); y++) {
            out(x, y, 0) = x + y*out.width();
            out(x, y, 1) = 0;
            out(x, y, 2) = 0;
        }
    }
    return out;
}

void test_img2eigen()
{
    Image im = make_simple_test_img();
    Matrix I = img2eigen(im);
    cout << "rows: " << I.rows() << endl;
    cout << "cols: " << I.cols() << endl;
    cout << I << endl;
}

void test_window()
{
    Image im = make_simple_test_img();
    Matrix I = img2eigen(im);
    Matrix W = get_window(I, 1, 2, 1, 4);
    cout << W << endl;
    cout << endl;

    cout << "Channel averages" << endl;
    Matrix mu = mu_window(W);
    cout << mu << endl;
    cout << endl;

    cout << "Channel covariance matrix" << endl;
    Matrix sigma = sigma_window(W, mu);
    cout << sigma << endl;
    cout << endl;
}

void test_laplacian()
{
    Image im = make_simple_test_img();
    Image known = im*0.0f;
    SparseMatrix L = matting_laplacian(im, known);
    cout << L << endl;
}

void test_scribble2trimap()
{
    Image im("./TestPhotos/dandelion_clipped.png");
    Image scribble("./TestPhotos/dandelion_clipped_m.png");

    Image trimap = scribble2trimap(im, scribble);
    trimap.write("./TestPhotos/dandelion_trimap.png");
}

void test_sky_ground_trimap()
{
    Image im("./TestPhotos/ec_test.png");

    Image trimap = sky_ground_trimap(im, 0.5, 0.05);
    trimap.write("./TestPhotos/ec_test_trimap.png");
}

void print_image(const Image &im)
{
    for(int y=0; y<im.height(); y++) {
        for(int x=0; x<im.width(); x++) {
            if(im.channels() < 3) {
                printf("%.2f ", im(x, y));
            } else {
                printf("(%.2f,%.2f,%.2f) ", im(x, y, 0), im(x, y, 1), im(x, y, 2));
            }
        }
        printf("\n");
    }
}

void test_cfm_simple()
{

    Image im = make_simple_test_img();
    cout << "Test image" << endl;
    print_image(im);

    Image trimap(im.width(), im.height(), 1);
    trimap = trimap * 0.0f + 0.5f;
    trimap(0, 0) = 1.0f;
    trimap(3, 3) = 0.0f;
    cout << "Test trimap" << endl;
    print_image(trimap);

    Image matte = cfm_with_trimap(im, trimap);
    cout << "Test matte" << endl;
    print_image(matte);
}

void test_cfm()
{

    cout << "testing cfm..." << endl;
    Image im("./TestPhotos/dandelion_clipped.png");
    Image scribble("./TestPhotos/dandelion_clipped_m.png");

    Image trimap = scribble2trimap(im, scribble);
    unsigned long s = millisecond_timer();
    Image matte  = cfm_with_trimap(im, trimap);
    unsigned long time = millisecond_timer() - s;
    cout << "time: " << time << endl;
    matte.write("./TestPhotos/dandelion_matte.png");

    int x=0, y=0;
    cout << "Values at (0, 0):" << endl;
    printf("(%.2f,%.2f,%.2f)\n", im(x, y, 0), im(x, y, 1), im(x, y, 2));
    printf("%.2f\n", trimap(x, y));
    printf("%.2f\n", matte(x, y));
}

void test_cfm_iterative()
{
    cout << "testing iterative cfm..." << endl;

    Image im("./TestPhotos/dandelion_clipped.png");
    Image scribble("./TestPhotos/dandelion_clipped_m.png");

    Image trimap = scribble2trimap(im, scribble);
    unsigned long s = millisecond_timer();
    Image matte  = iterative_cfm(im, trimap);
    unsigned long time = millisecond_timer() - s;
    cout << "time: " << time << endl;
    matte.write("./TestPhotos/dandelion_matte_iterative.png");

    int x=0, y=0;
    cout << "Values at (0, 0):" << endl;
    printf("(%.2f,%.2f,%.2f)\n", im(x, y, 0), im(x, y, 1), im(x, y, 2));
    printf("%.2f\n", trimap(x, y));
    printf("%.2f\n", matte(x, y));
}

void test_fg_bg()
{
    Image im("./TestPhotos/dandelion_clipped.png");
    Image scribble("./TestPhotos/dandelion_clipped_m.png");

    Image trimap = scribble2trimap(im, scribble);
    Image matte  = cfm_with_trimap(im, trimap);

    Image fg(im.width(), im.height(), im.channels());
    Image bg(im.width(), im.height(), im.channels());
    compute_fg_bg(im, matte, fg, bg);
    fg.write("./TestPhotos/dandelion_fg.png");
    bg.write("./TestPhotos/dandelion_bg.png");

    Image blend = blend_fg_bg(matte, fg, bg);
    blend.write("./TestPhotos/dandelion_blend.png");
}

void test_ec()
{
    Image im("./TestPhotos/ec_test.png");

    Image trimap = sky_ground_trimap(im, 0.6, 0.05);
    //Image matte  = cfm_with_trimap(im, trimap);
    //matte.write("./TestPhotos/ec_matte.png");
    Image matte2 = iterative_cfm(im, trimap, 0.05f);
    matte2.write("./TestPhotos/ec_matte_iterative.png");
}

void test_wallace()
{
    Image im("./TestPhotos/wallace_test.png");
    Image trimap = sky_ground_trimap(im, 0.25, 0.1);

    unsigned long s = millisecond_timer();
    Image matte  = cfm_with_trimap(im, trimap);
    unsigned long time = millisecond_timer() - s;
    cout << "time for direct solve: " << time << endl;
    matte.write("./TestPhotos/wallace_matte.png");

    s = millisecond_timer();
    Image matte2  = iterative_cfm(im, trimap, 0.0001f);
    time = millisecond_timer() - s;
    cout << "time for iterative solve: " << time << endl;
    matte2.write("./TestPhotos/wallace_matte_iterative.png");
}

void test_stata()
{
    Image im("./TestPhotos/stata_test.png");
    Image scribble("./TestPhotos/stata_test_scribble.png");

    Image trimap = scribble2trimap(im, scribble);
    Image matte  = iterative_cfm(im, trimap, 0.001);
    matte.write("./TestPhotos/stata_matte.png");

    //Image matte_th = threshold_matte(matte, 0.3, 0.9);
    //matte_th = gaussianBlur_separable(matte_th, 0.5);
    //matte_th.write("./TestPhotos/state_matte_th.png");

    Image fg(im.width(), im.height(), im.channels());
    Image bg(im.width(), im.height(), im.channels());
    compute_fg_bg(im, matte, fg, bg);
    fg.write("./TestPhotos/stata_fg.png");
    bg.write("./TestPhotos/stata_bg.png");

    Image blend = blend_fg_bg(matte, fg, bg);
    blend.write("./TestPhotos/stata_blend.png");

    Image dome("./TestPhotos/dome.png");
    Image blend2 = blend_fg_bg(matte, fg, dome);
    blend2.write("./TestPhotos/stata_dome.png");
}

void test_asb()
{
    Image im("./TestPhotos/asb_test.png");
    Image scribble("./TestPhotos/asb_test_scribble.png");

    Image trimap = scribble2trimap(im, scribble);
    Image matte = iterative_cfm(im, trimap);
    matte.write("./TestPhotos/asb_matte.png");

    Image fg(im.width(), im.height(), im.channels());
    Image bg(im.width(), im.height(), im.channels());
    compute_fg_bg(im, matte, fg, bg);
    fg.write("./TestPhotos/asb_fg.png");
    bg.write("./TestPhotos/asb_bg.png");

    Image blend = blend_fg_bg(matte, fg, bg);
    blend.write("./TestPhotos/asb_blend.png");
}

void test_find_remove_stars()
{
    // this requires first running test_wallace()
    //Image matte("./TestPhotos/ec_matte_iterative.png");
    Image matte("./TestPhotos/wallace_matte.png");
    Image stars = find_stars_matte(matte);
    visualize_stars(stars).write("./TestPhotos/wallace_matte_stars.png");
    Image no_stars = remove_stars_from_matte(matte, stars);
    no_stars.write("./TestPhotos/wallace_matte_stars_removed.png");

    /*
    auto corners = HarrisCorners(matte);
    visualizeCorners(matte, corners).write("./TestPhotos/wallace_corners.png");
    auto features = computeFeatures(matte, corners);
    visualizeFeatures(matte, features).write("./TestPhotos/wallace_features.png");
    */
}

void test_homography_svd()
{
    int w = 143;
    int h = 66;

    // from pset 6
    CorrespondencePair corresp[4] = {
        CorrespondencePair( 0,0,1, 96,171,1),
        CorrespondencePair( w-1,0,1, 235,174,1),
        CorrespondencePair( w-1,h-1,1, 235,232,1),
        CorrespondencePair( 0,h-1,1, 95,238,1)
    };

    // Compute homography.
    Matrix Hcomputed = computeHomography(corresp);
    cout << "Computed homography for bus example from given point pairs" << endl;
    cout << Hcomputed << endl;

    // Compute homography with SVD
    Matrix Hcomputed_svd = computeHomographySVD(corresp);
    cout << "Computed homography for bus example from given point pairs" << endl;
    cout << Hcomputed_svd << endl;
    cout << "normalized to bottom right = 1" << endl;
    cout << Hcomputed_svd / Hcomputed_svd(2, 2) << endl;
}

void test_correspondence()
{
    Image stars1("./TestPhotos/DSC_3659_bg.png");
    Image stars2("./TestPhotos/DSC_3660_bg.png");

    // try finding stars
    //cout << "test finding stars" << endl;
    //Image stars = find_stars_color(stars1);
    //visualize_stars(stars).write("./TestPhotos/stars1-stars.png");

    // compute harris corners
    cout << "computing harris corners" << endl;
    vector<Point> h1 = HarrisCorners_stars(stars1);
    vector<Point> h2 = HarrisCorners_stars(stars2);
    visualizeCorners(stars1, h1).write("./TestPhotos/stars1-corners.png");
    visualizeCorners(stars2, h2).write("./TestPhotos/stars2-corners.png");

    // compute features
    cout << "computing features" << endl;
    float sigmaBlurDescriptor = 2.0f;
    float radiusDescriptor = 5;
    vector<Feature> f1 = computeFeatures(stars1, h1, sigmaBlurDescriptor, radiusDescriptor);
    vector<Feature> f2 = computeFeatures(stars2, h2, sigmaBlurDescriptor, radiusDescriptor);
    visualizeFeatures(stars1, f1).write("./TestPhotos/stars1-features.png");
    visualizeFeatures(stars2, f2).write("./TestPhotos/stars2-features.png");

    cout << "computing correspondences" << endl;
    //vector<FeatureCorrespondence> corr = findCorrespondences(f1, f2, 1.2f);
    vector<FeatureCorrespondence> corr = findCorrespondences_stars(f1, f2, 1000.0f, 1000.0f);
    visualizePairs(stars1, stars2, corr).write("./TestPhotos/stars-featcorr.png");

    // compute homography
    cout << "computing homography" << endl;
    float epsilon = 4;
    Matrix H = RANSAC(corr, 10000, epsilon);
    vector<bool> ins = inliers(H, corr, epsilon);
    visualizePairsWithInliers(stars1, stars2, corr, ins).write("./TestPhotos/stars-RANSAC-featcorr-inliers.png");
}

void test_align_stars()
{
    Image stars1("./TestPhotos/DSC_3659_bg.png");
    Image stars2("./TestPhotos/DSC_3660_bg.png");
    Image stars3("./TestPhotos/DSC_3661_bg.png");

    // align stars1 to stars2
    Matrix H1 = compute_star_homography(stars1, stars2);
    Image stars1_aligned = align_with_homography(stars1, H1);

    // align stars3 to stars2
    Matrix H3 = compute_star_homography(stars3, stars2);
    Image stars3_aligned = align_with_homography(stars3, H3);

    stars1_aligned.write("./TestPhotos/align-stars1.png");
    stars3_aligned.write("./TestPhotos/align-stars3.png");
    stars2.write("./TestPhotos/align-stars2.png");
}

void test_stack_stars()
{
    Image stars1("./TestPhotos/DSC_3659_bg.png");
    Image stars2("./TestPhotos/DSC_3660_bg.png");
    Image stars3("./TestPhotos/DSC_3661_bg.png");

    cout << "stacking unaligned images" << endl;
    std::vector<Image> unaligned{stars1, stars2, stars3};
    Image stack_with_max  = stack_images(unaligned, stack_max);
    Image stack_with_mean = stack_images(unaligned, stack_mean);
    stack_with_max.write("./TestPhotos/stack_unaligned_max.png");
    stack_with_mean.write("./TestPhotos/stack_unaligned_mean.png");

    cout << "aligning stars in images" << endl;
    // align stars1 to stars2
    Matrix H1 = compute_star_homography(stars1, stars2);
    Image stars1_aligned = align_with_homography(stars1, H1);

    // align stars3 to stars2
    Matrix H3 = compute_star_homography(stars3, stars2);
    Image stars3_aligned = align_with_homography(stars3, H3);

    cout << "stacking aligned images" << endl;
    std::vector<Image> aligned{stars1_aligned, stars2, stars3_aligned};
    stack_with_max  = stack_images(aligned, stack_max);
    stack_with_mean = stack_images(aligned, stack_mean);
    stack_with_max.write("./TestPhotos/stack_aligned_max.png");
    stack_with_mean.write("./TestPhotos/stack_aligned_mean.png");
}
