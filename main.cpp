#include <iostream>
#include <iomanip>

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


void matte_ec()
{
    Image scribble("./Photos/ec/scribble.png");
    Image im_ref("./Photos/ec/DSC_3659.png");
    Image trimap = scribble2trimap(im_ref, scribble);

    for(int i=3659; i<=3674; i++) {
        cout << "processing image i=" << i << endl;
        std::ostringstream fname;
        fname << "./Photos/ec/DSC_";
        fname << std::setfill('0') << std::setw(4);
        fname << i;

        Image im(fname.str() + ".png");
        Image matte = iterative_cfm(im, trimap);
        matte.write(fname.str() + "_matte.png");

        Image fg(im.width(), im.height(), im.channels());
        Image bg(im.width(), im.height(), im.channels());
        compute_fg_bg(im, matte, fg, bg);
        fg.write(fname.str() + "_fg.png");
        bg.write(fname.str() + "_bg.png");
    }
}

void matte_stata()
{
    Image scribble("./Photos/stata/scribble.png");
    Image im_ref("./Photos/stata/foreground.png");
    Image trimap = scribble2trimap(im_ref, scribble);

    std::string filename = "./Photos/stata/foreground";
    Image matte_ref = iterative_cfm(im_ref, trimap);
    matte_ref.write(filename + "_matte.png");

    Image fg_ref(im_ref.width(), im_ref.height(), im_ref.channels());
    Image bg_ref(im_ref.width(), im_ref.height(), im_ref.channels());
    compute_fg_bg(im_ref, matte_ref, fg_ref, bg_ref);
    fg_ref.write(filename + "_fg.png");
    bg_ref.write(filename + "_bg.png");

    for(int i=3781; i<=3791; i++) {
        cout << "processing image i=" << i << endl;
        std::ostringstream fname;
        fname << "./Photos/stata/DSC_";
        fname << std::setfill('0') << std::setw(4);
        fname << i;

        Image im(fname.str() + ".png");
        Image matte = iterative_cfm(im, trimap);
        matte.write(fname.str() + "_matte.png");

        Image fg(im.width(), im.height(), im.channels());
        Image bg(im.width(), im.height(), im.channels());
        compute_fg_bg(im, matte, fg, bg);
        fg.write(fname.str() + "_fg.png");
        bg.write(fname.str() + "_bg.png");
    }
}

void matte_stata_one_step()
{
    Image scribble("./Photos/stata_2/scribble.png");
    Image im_ref("./Photos/stata_2/foreground.png");
    Image trimap = scribble2trimap(im_ref, scribble);

    cout << "processing image foreground" << endl;
    std::string filename = "./Photos/stata_2/foreground";
    Image matte_ref = cfm_with_trimap(im_ref, trimap);
    matte_ref.write(filename + "_matte.png");

    Image fg_ref(im_ref.width(), im_ref.height(), im_ref.channels());
    Image bg_ref(im_ref.width(), im_ref.height(), im_ref.channels());
    compute_fg_bg(im_ref, matte_ref, fg_ref, bg_ref);
    fg_ref.write(filename + "_fg.png");
    bg_ref.write(filename + "_bg.png");

    for(int i=3781; i<=3791; i++) {
        cout << "processing image i=" << i << endl;
        std::ostringstream fname;
        fname << "./Photos/stata_2/DSC_";
        fname << std::setfill('0') << std::setw(4);
        fname << i;

        Image im(fname.str() + ".png");
        Image matte = cfm_with_trimap(im, trimap);
        matte.write(fname.str() + "_matte.png");

        Image fg(im.width(), im.height(), im.channels());
        Image bg(im.width(), im.height(), im.channels());
        compute_fg_bg(im, matte, fg, bg);
        fg.write(fname.str() + "_fg.png");
        bg.write(fname.str() + "_bg.png");
    }
}

void stack_ec()
{
    int start_value = 3659;
    //int end_value   = 3661;
    int end_value   = 3674;
    int reference_value = (end_value+start_value)/2;
    int reference_index = reference_value - start_value;
    std::string filename_prefix = "./Photos/ec/DSC_";
    std::string output_prefix   = "./Photos/ec/";

    // load background images
    std::vector<Image> backgrounds;
    for(int i=start_value; i<=end_value; i++) {
        std::ostringstream fname;
        fname << filename_prefix;
        fname << std::setfill('0') << std::setw(4);
        fname << i;
        backgrounds.push_back(Image(fname.str() + "_bg.png"));
    }

    // load foreground and matte for middle image
    std::ostringstream fname_fg;
    fname_fg << filename_prefix;
    fname_fg << std::setfill('0') << std::setw(4);
    fname_fg << reference_value;
    Image foreground(fname_fg.str() + "_fg.png");
    Image matte(fname_fg.str() + "_matte.png");

    // stack stars without alignment
    cout << "stacking unaligned images" << endl;
    Image stack_unaligned = stack_images(backgrounds, stack_max);
    stack_unaligned.write(output_prefix + "stack_bg_unaligned.png");
    blend_fg_bg(matte, foreground, stack_unaligned).write(output_prefix + "stack_unaligned.png");

    // compute homographies needed to align stars
    cout << "aligning images" << endl;
    auto h_pairwise = compute_pairwise_homographies(backgrounds, reference_index);
    auto h_absolute = compute_absolute_homographies(h_pairwise, reference_index);

    // stack stars with alignment
    cout << "stacking aligned images" << endl;
    Image stack_aligned = stack_images_with_homographies(backgrounds, h_absolute, stack_mean);
    stack_aligned.write(output_prefix + "stack_bg_aligned.png");
    blend_fg_bg(matte, foreground, stack_aligned).write(output_prefix + "stack_aligned.png");
    Image stack_aligned_max = stack_images_with_homographies(backgrounds, h_absolute, stack_max);
    stack_aligned_max.write(output_prefix + "stack_bg_aligned_max.png");
    blend_fg_bg(matte, foreground, stack_aligned_max).write(output_prefix + "stack_aligned_max.png");
}

void stack_stata()
{
    int start_value = 3781;
    int end_value   = 3791;
    int reference_value = (end_value+start_value)/2;
    int reference_index = reference_value - start_value;
    std::string filename_prefix = "./Photos/stata/DSC_";
    std::string output_prefix   = "./Photos/stata/";

    // load background images
    std::vector<Image> backgrounds;
    for(int i=start_value; i<=end_value; i++) {
        std::ostringstream fname;
        fname << filename_prefix;
        fname << std::setfill('0') << std::setw(4);
        fname << i;
        backgrounds.push_back(Image(fname.str() + "_bg.png"));
    }

    // load foreground and matte for middle image
    Image foreground(output_prefix + "foreground_fg.png");
    Image matte(output_prefix + "foreground_matte.png");

    // stack stars without alignment
    cout << "stacking unaligned images" << endl;
    Image stack_unaligned = stack_images(backgrounds, stack_max);
    stack_unaligned.write(output_prefix + "stack_bg_unaligned.png");
    blend_fg_bg(matte, foreground, stack_unaligned).write(output_prefix + "stack_unaligned.png");

    // compute homographies needed to align stars
    cout << "aligning images" << endl;
    auto h_pairwise = compute_pairwise_homographies(backgrounds, reference_index);
    auto h_absolute = compute_absolute_homographies(h_pairwise, reference_index);

    // stack stars with alignment
    cout << "stacking aligned images" << endl;
    Image stack_aligned = stack_images_with_homographies(backgrounds, h_absolute, stack_mean);
    stack_aligned.write(output_prefix + "stack_bg_aligned.png");
    blend_fg_bg(matte, foreground, stack_aligned).write(output_prefix + "stack_aligned.png");
    Image stack_aligned_max = stack_images_with_homographies(backgrounds, h_absolute, stack_max);
    stack_aligned_max.write(output_prefix + "stack_bg_aligned_max.png");
    blend_fg_bg(matte, foreground, stack_aligned_max).write(output_prefix + "stack_aligned_max.png");
}

void stack_stata_one_step()
{
    int start_value = 3781;
    int end_value   = 3791;
    int skip = 3786;
    int reference_value = 3785;
    int reference_index = reference_value - start_value;
    std::string filename_prefix = "./Photos/stata_2/DSC_";
    std::string output_prefix   = "./Photos/stata_2/";

    // load background images
    std::vector<Image> backgrounds;
    for(int i=start_value; i<=end_value; i++) {
        if(i == skip)
            continue;
        std::ostringstream fname;
        fname << filename_prefix;
        fname << std::setfill('0') << std::setw(4);
        fname << i;
        backgrounds.push_back(Image(fname.str() + "_bg.png"));
    }

    // load foreground and matte for middle image
    Image foreground(output_prefix + "foreground_fg.png");
    Image background(output_prefix + "foreground_bg.png");
    Image matte(output_prefix + "foreground_matte.png");

    // stack stars without alignment
    cout << "stacking unaligned images" << endl;
    Image stack_unaligned = stack_images(backgrounds, stack_max);
    stack_unaligned.write(output_prefix + "stack_bg_unaligned.png");
    blend_fg_bg(matte, foreground, stack_unaligned).write(output_prefix + "stack_unaligned.png");

    // compute homographies needed to align stars
    cout << "aligning images" << endl;
    auto h_pairwise = compute_pairwise_homographies(backgrounds, reference_index);
    auto h_absolute = compute_absolute_homographies(h_pairwise, reference_index);

    // stack stars with alignment
    cout << "stacking aligned images" << endl;
    Image stack_aligned = stack_images_with_homographies(backgrounds, h_absolute, stack_mean);
    stack_aligned.write(output_prefix + "stack_bg_aligned.png");
    blend_fg_bg(matte, foreground, stack_aligned).write(output_prefix + "stack_aligned.png");

    Image stack_aligned_max = stack_images_with_homographies(backgrounds, h_absolute, stack_max);
    stack_aligned_max.write(output_prefix + "stack_bg_aligned_max.png");
    blend_fg_bg(matte, foreground, stack_aligned_max).write(output_prefix + "stack_aligned_max.png");

    // stack the pure background parts of the background images
    matte = threshold_matte(matte, 0.1, 0.9);
    matte.write(output_prefix + "foreground_matte_th.png");
    stack_aligned = stack_images_with_homographies_and_matte(backgrounds, h_absolute, matte, stack_mean);

    // copy over unstacked pixels from a single background image
    for(int y=0; y<matte.height(); y++) {
        for(int x=0; x<matte.width(); x++) {
            if(matte(x, y) > 0.0f) {
                for(int c=0; c<stack_aligned.channels(); c++) {
                    stack_aligned(x, y, c) = background(x, y, c);
                }
            }
        }
    }

    stack_aligned.write(output_prefix + "stack_bg_aligned_with_matte.png");
    blend_fg_bg(matte, foreground, stack_aligned).write(output_prefix + "stack_aligned_with_matte.png");
}

void merge_stata_edited()
{
    std::string prefix = "./Photos/stata/";

    Image background(prefix + "stack_bg_aligned_edited.png");
    Image foreground(prefix + "foreground_fg.png");
    Image matte(prefix + "foreground_matte.png");

    blend_fg_bg(matte, foreground, background).write(prefix + "stack_aligned_edited.png");
}

void merge_stata_one_step_edited()
{
    std::string prefix = "./Photos/stata_2/";

    Image background(prefix + "stack_bg_aligned_edited.png");
    Image foreground(prefix + "foreground_fg.png");
    Image matte(prefix + "foreground_matte.png");

    blend_fg_bg(matte, foreground, background).write(prefix + "stack_aligned_edited.png");

    matte = threshold_matte(matte, 0.1f, 0.9f);
    Image background2(prefix + "stack_bg_aligned_with_matte_edited.png");
    blend_fg_bg(matte, foreground, background2).write(prefix + "stack_aligned_with_matte_edited.png");
}

void merge_ec_edited()
{
    std::string prefix = "./Photos/ec/";

    Image background(prefix + "stack_bg_aligned_edited.png");
    Image foreground(prefix + "DSC_3666_fg.png");
    Image matte(prefix + "DSC_3666_matte.png");

    blend_fg_bg(matte, foreground, background).write(prefix + "stack_aligned_edited.png");
}

void test_correspondence_stata2()
{
    Image stars1("./Photos/stata_2/DSC_3785_bg.png");
    Image stars2("./Photos/stata_2/DSC_3787_bg.png");

    //stars1 = threshold_matte(stars1, 0.1f, 1.0f);
    //stars2 = threshold_matte(stars2, 0.1f, 1.0f);

    // try finding stars
    //cout << "test finding stars" << endl;
    //Image stars = find_stars_color(stars1, 3.0f);
    //visualize_stars(stars).write("./Output/stars1-stars.png");

    // compute harris corners
    cout << "computing harris corners" << endl;
    vector<Point> h1 = HarrisCorners_stars(stars1);
    vector<Point> h2 = HarrisCorners_stars(stars2);
    visualizeCorners(stars1, h1).write("./Output/stars1-corners.png");
    visualizeCorners(stars2, h2).write("./Output/stars2-corners.png");

    // compute features
    cout << "computing features" << endl;
    float sigmaBlurDescriptor = 2.0f;
    float radiusDescriptor = 5;
    vector<Feature> f1 = computeFeatures(stars1, h1, sigmaBlurDescriptor, radiusDescriptor);
    vector<Feature> f2 = computeFeatures(stars2, h2, sigmaBlurDescriptor, radiusDescriptor);
    visualizeFeatures(stars1, f1).write("./Output/stars1-features.png");
    visualizeFeatures(stars2, f2).write("./Output/stars2-features.png");

    cout << "computing correspondences" << endl;
    //vector<FeatureCorrespondence> corr = findCorrespondences(f1, f2, 1.2f);
    vector<FeatureCorrespondence> corr = findCorrespondences_stars(f1, f2, 1000.0f, 1000.0f);
    visualizePairs(stars1, stars2, corr).write("./Output/stars-featcorr.png");

    // compute homography
    cout << "computing homography" << endl;
    float epsilon = 1;
    Matrix H = RANSAC(corr, 50000, epsilon);
    vector<bool> ins = inliers(H, corr, epsilon);
    visualizePairsWithInliers(stars1, stars2, corr, ins).write("./Output/stars-RANSAC-featcorr-inliers.png");

    stars2.write("./Output/stars2.png");
    align_with_homography(stars1, H).write("./Output/stars1-aligned.png");
}


int main()
{
    // Test your intermediate functions
    //test_img2eigen();
    //test_window();
    //test_laplacian();
    //test_scribble2trimap();
    //test_sky_ground_trimap();
    //test_cfm_simple();
    //test_cfm();
    //test_cfm_iterative();
    //test_fg_bg();
    //test_ec();
    //test_wallace();
    //test_stata();
    //test_asb();
    //test_find_remove_stars();
    //test_homography_svd();
    //test_correspondence();
    //test_align_stars();
    //test_stack_stars();

    // compute mattes for real
    // these functions are really slow (tens of minutes)
    //matte_ec();
    //matte_stata();
    //matte_stata_one_step();

    // stack for real
    // these functions are really slow (tens of minutes)
    //stack_ec();
    //stack_stata();
    //test_correspondence_stata2();
    //stack_stata_one_step();

    //merge_stata_edited();
    //merge_ec_edited();
    //merge_stata_one_step_edited();

    return EXIT_SUCCESS;
}
