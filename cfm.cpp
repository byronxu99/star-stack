#include "cfm.h"

using std::cout;
using std::endl;

// convert a WxHxC image to a (WH)xC 2D matrix
Matrix img2eigen(const Image &im)
{
    int w = im.width();
    int h = im.height();
    int c = im.channels();
    Matrix I(w*h, c);

    for(int x=0; x<w; x++) {
        for(int y=0; y<h; y++) {
            for(int z=0; z<c; z++) {
                I(x + y*w, z) = im(x, y, z);
            }
        }
    }

    return I;
}

// return a (window area)x(channels) matrix around pixel i,j of a matrix representing
// an image of width
// note: function does not check bounds
Matrix get_window(const Matrix &I, int i, int j, int radius, int image_width)
{
    int window_size = 2*radius + 1;
    int window_area = window_size * window_size;
    int channels = I.cols();
    Matrix W(window_area, channels);

    int current_row = 0;
    for(int dj=-radius; dj<=radius; dj++) {
        for(int di=-radius; di<=radius; di++) {
            for(int c=0; c<channels; c++) {
                W(current_row, c) = I((i+di) + (j+dj)*image_width, c);
            }
            current_row++;
        }
    }

    return W;
}

// compute average for each channel, returning a 1xC matrix
Matrix mu_window(const Matrix &W)
{
    Matrix sums = W.colwise().sum();
    return sums / W.rows();
}

// subtract mu from each row of W
Matrix window_minus_mu(const Matrix &W, const Matrix &mu)
{
    // convert mu to a Vector
    Vector mu_v(mu.cols());
    for(int c=0; c<mu.cols(); c++) {
        mu_v(c) = mu(0, c);
    }

    // convert mu_v to a column vector and
    // subtract it from each row of W
    return W.rowwise() - mu_v.transpose();
}

// compute CxC covariance matrix
Matrix sigma_window(const Matrix &W, const Matrix &mu)
{
    Matrix W_minus_mu = window_minus_mu(W, mu);

    Matrix sigma = W_minus_mu.transpose() * W_minus_mu;
    sigma = sigma / (sigma.rows() * sigma.cols());
    return sigma;
}

// compute the matting laplacian matrix given image and mask of known pixels
SparseMatrix matting_laplacian(const Image &im, const Image &known, float epsilon, int radius)
{
    if(im.width() != known.width() || im.height() != known.height()) {
        throw InvalidArgument();
    }

    int width = im.width();
    int height = im.height();
    int N = width*height; // number of pixels
    int window_size = 2*radius + 1;
    int window_area = window_size * window_size;
    int channels = im.channels();

    // list of sparse matrix entries, where duplicates will be automatically added together
    std::vector<Triplet> triplets;
    // preallocate enough memory
    triplets.reserve(N*window_area);

    // create matrix from image
    Matrix I = img2eigen(im);

    // loop over all possible window centers
    for(int wy=radius; wy<height-radius; wy++) {
        for(int wx=radius; wx<width-radius; wx++) {
            // discard windows that don't contain at least one unknown pixel
            bool has_unknown = false;
            for(int dy=-radius; dy<=radius; dy++)
                for(int dx=-radius; dx<=radius; dx++)
                    if(known(wx+dx, wy+dy) < epsilon)
                        has_unknown = true;
            if(!has_unknown) continue;

            // get pixels in window
            Matrix W = get_window(I, wx, wy, radius, width);

            // compute mean and covariance
            Matrix mu = mu_window(W);
            Matrix sigma = sigma_window(W, mu);

            // compute window minus average
            // this is the equivalent of (I-mu) in equation 12 of the paper
            Matrix W_minus_mu = window_minus_mu(W, mu);

            // compute the inverted term in equation 12 of the paper
            // (sigma + (epsilon/|w|) I3)^-1
            Matrix inv_term = (sigma + (epsilon/(float)window_area) * Matrix::Identity(channels, channels)).inverse();

            // compute the matrix product in equation 12 of the paper
            // (I-mu) * inv_term * (I-mu)
            Matrix mat_product = W_minus_mu * inv_term * W_minus_mu.transpose();

            // compute the entire expression inside the summation
            // the actual summation is done automatically through
            // duplicate entries in the sparse triplet list
            Matrix sum_expr = Matrix::Identity(window_area, window_area) - (1.0f/(float)window_area) * (Matrix::Ones(window_area, window_area) + mat_product);

            // add triplets to list
            // there is an entry in sum_expr for each pair of possible offsets from the window center
            for(int dy=-radius; dy<=radius; dy++) {
                for(int dx=-radius; dx<=radius; dx++) {
                    for(int dy2=-radius; dy2<=radius; dy2++) {
                        for(int dx2=-radius; dx2<=radius; dx2++) {
                            // output row and column in the matting laplacian
                            int L_row = (wx+dx)  + (wy+dy) *width;
                            int L_col = (wx+dx2) + (wy+dy2)*width;
                            // input row and column in sum_expr
                            int s_row = (dx+radius)  + (dy+radius) *window_size;
                            int s_col = (dx2+radius) + (dy2+radius)*window_size;
                            // triplet
                            triplets.emplace_back(L_row, L_col, sum_expr(s_row, s_col));
                        }
                    }
                }
            }

        }
    }

    // create the sparse matrix
    SparseMatrix L(N, N);
    L.setFromTriplets(triplets.begin(), triplets.end());
    return L;
}

// compute matte given image and trimap
Image cfm_with_trimap(const Image &im, const Image &trimap, float lambda, float epsilon, int radius)
{
    if(im.width() != trimap.width() || im.height() != trimap.height())
        throw InvalidArgument();

    int width = im.width();
    int height = im.height();
    int N = width*height; // number of pixels

    // convert trimap to known pixel mask
    Image known(width, height, 1);
    for(int y=0; y<height; y++) {
        for(int x=0; x<width; x++) {
            if(trimap(x, y) < epsilon || trimap(x, y) >= 1-epsilon) {
                known(x, y) = 1.0f;
            } else {
                known(x, y) = 0.0f;
            }
        }
    }

    //cout << "computing laplacian" << endl;

    // compute sparse matting laplacian matrix
    SparseMatrix L = matting_laplacian(im, known, epsilon, radius);

    // compute bS in equation 14 of the paper
    //cout << "computing linear system" << endl;
    Matrix known_mat  = img2eigen(known);
    Matrix trimap_mat = img2eigen(trimap);
    Vector bS = known_mat.cwiseProduct(trimap_mat).col(0);

    // compute lambda*DS in equation 14 of the paper
    // as a sparse matrix
    Vector lambda_DS = (lambda*known_mat).col(0);
    std::vector<Triplet> triplets;
    for(int i=0; i<lambda_DS.size(); i++) {
        triplets.emplace_back(i, i, lambda_DS(i));
    }
    SparseMatrix lambda_DS_sparse(N, N);
    lambda_DS_sparse.setFromTriplets(triplets.begin(), triplets.end());

    // create the linear system lhs*alpha = rhs
    SparseMatrix lhs = L + lambda_DS_sparse;
    lhs.makeCompressed();
    Vector rhs = lambda*bS;

    // solve the sparse linear system
    cout << "solving linear system of size " << lhs.nonZeros() << endl;
    //Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<float>> solver;
    //Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper> solver;
    Eigen::SparseLU<SparseMatrix> solver;
    //solver.setMaxIterations(5000);

    solver.compute(lhs);
    Vector alpha = solver.solve(rhs);
    //Vector alpha = solver.solveWithGuess(rhs, trimap_mat.col(0));

    // clamp values
    //cout << "creating matte" << endl;
    alpha = alpha.array().min(1.0f).max(0.0f);

    // create output
    Image matte(width, height, 1);
    for(int y=0; y<height; y++) {
        for(int x=0; x<width; x++) {
            matte(x, y) = alpha(x + y*width);
        }
    }

    return matte;
}

// convert a set of scribbles to a trimap
Image scribble2trimap(const Image &im, const Image &scribble)
{
    static const float eps = 0.1;
    Image difference = im - scribble;
    Image trimap(im.width(), im.height(), 1);

    for(int y=0; y<im.height(); y++) {
        for(int x=0; x<im.width(); x++) {
            float sum_scribble = 0.0f;
            float sum_difference = 0.0f;
            for(int c=0; c<im.channels(); c++) {
                sum_scribble += scribble(x, y, c);
                sum_difference += abs(difference(x, y, c));
            }

            // if no difference between original and scribble, it's not scribbled over
            if(sum_difference < 1e-6) {
                trimap(x, y) = 0.5f;
                continue;
            }

            // if all channels are near zero, it's a black scribble
            if(sum_scribble < eps)
                trimap(x, y) = 0.0f;
            // if all channels are near one, it's a white scribble
            else if(sum_scribble > im.channels()-eps)
                trimap(x, y) = 1.0f;
            // otherwise it's not a scribble
            else
                trimap(x, y) = 0.5f;
        }
    }

    return trimap;
}

// generate a trimap guess for landscape photos
Image sky_ground_trimap(const Image &im, float frac_sky, float frac_ground)
{
    int sky_limit = frac_sky * im.height();
    int ground_limit = im.height() - (frac_ground * im.height());
    Image trimap(im.width(), im.height(), 1);

    for(int y=0; y<im.height(); y++) {
        for(int x=0; x<im.width(); x++) {
            if(y <= sky_limit)
                trimap(x, y) = 0.0f;
            else if(y >= ground_limit)
                trimap(x, y) = 1.0f;
            else
                trimap(x, y) = 0.5f;
        }
    }

    return trimap;
}

Image threshold_matte(const Image &matte, float lo, float hi)
{
    Image out(matte.width(), matte.height(), 1);
    for(int y=0; y<matte.height(); y++) {
        for(int x=0; x<matte.width(); x++) {
            if(matte(x, y) < lo)
                out(x, y) = 0.0f;
            else if(matte(x, y) > hi)
                out(x, y) = 1.0f;
            else
                out(x, y) = matte(x, y);
        }
    }
    return out;
}

Image iterative_cfm(const Image &im, const Image &trimap, float threshold, float decimation, int max_size, float lambda, float epsilon, int radius)
{
    if(im.width() != trimap.width() || im.height() != trimap.height())
        throw InvalidArgument();

    int width = im.width();
    int height = im.height();
    int N = width*height; // number of pixels

    // for small enough images, directly solve it
    if(N <= max_size) {
        return cfm_with_trimap(im, trimap, lambda, epsilon, radius);
    }

    // otherwise, scale image smaller and solve recursively
    Image im_half = scaleLin(im, decimation);
    Image trimap_half = scaleLin(trimap, decimation);
    Image alpha_half = iterative_cfm(im_half, trimap_half, threshold, decimation, max_size, lambda, epsilon, radius);

    // given solution to image at half scale
    // scale up, threshold known pixels, and solve ambiguous pixels
    Image alpha_half_double = scaleLin(alpha_half, 1.0f/decimation);
    Image new_alpha(im.width(), im.height(), 1);
    for(int y=0; y<im.height(); y++) {
        for(int x=0; x<im.width(); x++) {
            // copy over interpolated alpha values
            new_alpha(x, y) = alpha_half_double.smartAccessor(x, y, 0, true);
            // copy over known alpha values
            if(trimap(x, y) < epsilon || trimap(x, y) >= 1-epsilon) {
                new_alpha(x, y) = trimap(x, y);
            }
        }
    }
    Image new_trimap = threshold_matte(new_alpha, threshold, 1.0f-threshold);
    return cfm_with_trimap(im, new_trimap, lambda, epsilon, radius);
}

// compute foreground and background given matte
// this is a very simple approximation instead of the method described in the paper
void compute_fg_bg(const Image &im, const Image &matte, Image &fg_out, Image &bg_out, float fg_threshold, float bg_threshold)
{
    if(im.width() != matte.width() || im.height() != matte.height())
        throw InvalidArgument();
    if(im.width() != fg_out.width() || im.height() != fg_out.height() || im.channels() != fg_out.channels())
        throw InvalidArgument();
    if(im.width() != bg_out.width() || im.height() != bg_out.height() || im.channels() != bg_out.channels())
        throw InvalidArgument();

    int width = im.width();
    int height = im.height();
    int channels = im.channels();

    Image matte_th = threshold_matte(matte, bg_threshold, fg_threshold);

    // make the background image by alpha-blending
    // between image color and black
    for(int c=0; c<channels; c++) {
        for(int y=0; y<height; y++) {
            for(int x=0; x<width; x++) {
                bg_out(x, y, c) = (1-matte_th(x, y)) * im(x, y, c); // + matte(x, y) * 0
            }
        }
    }

    // make the foreground image by solving the matting equation
    // (equation 1 in the paper) for all foreground pixels
    // F = (I - (1-alpha)B) / alpha
    for(int c=0; c<channels; c++) {
        for(int y=0; y<height; y++) {
            for(int x=0; x<width; x++) {
                // solve for foreground value given matte and background
                if(matte(x, y) > 0.0f) {
                    fg_out(x, y, c) = (im(x, y, c) - (1.0f-matte(x, y)) * bg_out(x, y, c)) / (matte(x, y));
                    fg_out(x, y, c) = std::min(std::max(fg_out(x, y, c), 0.0f), 1.0f);
                } else {
                    fg_out(x, y, c) = 1.0f;
                }
            }
        }
    }

    /*
    // make the foreground image by alpha-blending
    // between image color and white
    for(int c=0; c<channels; c++) {
        for(int y=0; y<height; y++) {
            for(int x=0; x<width; x++) {
                fg_out(x, y, c) = matte_th(x, y) * im(x, y, c) + (1.0f-matte_th(x, y)) * 1.0f;
            }
        }
    }
    */
}

// blend foreground and background given a matte
// this is a simple application of equation 1 in the paper
Image blend_fg_bg(const Image &matte, const Image &fg, const Image &bg)
{
    if(matte.width() != fg.width() || matte.height() != fg.height())
        throw InvalidArgument();
    if(matte.width() != bg.width() || matte.height() != bg.height())
        throw InvalidArgument();
    if(fg.channels() != bg.channels())
        throw InvalidArgument();

    int width = matte.width();
    int height = matte.height();
    int channels = fg.channels();

    Image out(width, height, channels);
    for(int c=0; c<channels; c++) {
        for(int y=0; y<height; y++) {
            for(int x=0; x<width; x++) {
                out(x, y, c) = matte(x, y) * fg(x, y, c) + (1.0f-matte(x, y)) * bg(x, y, c);
            }
        }
    }
    return out;
}

/* average color

    // we assume that the background has a near-uniform color
    // compute average color of background
    float color[channels] = { 0 };
    float count[channels] = { 0 };
    for(int c=0; c<channels; c++) {
        for(int y=0; y<height; y++) {
            for(int x=0; x<width; x++) {
                color[c]  += bg_mask(x, y) * im(x, y, c);
                count[c] += bg_mask(x, y);
            }
        }
        color[c] /= count[c];
    }

    for(int c=0; c<channels; c++)
        color[c] = 0.0f;

    // make the background image by alpha-blending
    // between image color and average color of background pixels
    for(int c=0; c<channels; c++) {
        for(int y=0; y<height; y++) {
            for(int x=0; x<width; x++) {
                if(bg_mask(x, y) > 0.0f) {
                    float amount_bg  = (bg_threshold - matte(x, y)) / bg_threshold;
                    float amount_avg = matte(x, y) / bg_threshold;
                    bg_out(x, y, c) = amount_bg * im(x, y, c) + amount_avg * color[c];
                } else {
                    // fill non-background pixels with average color
                    bg_out(x, y, c) = color[c];
                }
            }
        }
    }
*/
