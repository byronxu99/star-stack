# Closed-Form Matting and Astrophotography
Final project for 6.815 Digital and Computational Photography, Spring 2019.

This is a C++ implementation of the Closed-Form Image Matting algorithm by [Levin, Lischinski, and Weiss](https://sites.fas.harvard.edu/~cs278/papers/matting.pdf), applied to the context of astrophotography image stacking. The code involves an extensive amount of matrix math and sparse linear system solving, done using the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) linear algebra library.

## About
Taking pictures of the night sky is difficult because of low light levels, which leads to photos with high levels of [image noise](https://en.wikipedia.org/wiki/Image_noise). Usually, noise can be reduced by longer exposure times. But in astrophotography, the Earth's rotation causes stars to move relative to the ground, creating star trails when exposure times are too long. The solution involves *stacking*, or aligning and averaging many short exposures.

Good photos of the night sky include both stars and some landscape on the ground, which again gives rise to the issue of moving stars. If we align individual exposures with respect to the ground, stars will again form streaks. If we align stars, the ground becomes blurry. What's necessary here is to first separate ground from sky.

## Closed-Form Matting
Image matting is the problem of separating foreground from background, given small sets of pixels known to be in each category. The goal is to generate a *matte*, a grayscale image representing the opacity of the foreground at each pixel of the original image. Many approaches exist for this problem, usually involving assumptions of local smoothness within the foreground or background, and that foreground and background pixels are drawn from disjoint color distributions. Closed-form matting is unique in that it derives a quadratic cost function from these assumptions. Setting the derivative of this expression to zero results in a sparse system of linear equations, solveable to produce a globally optimal solution.

<p align="center">
  <img src="https://github.com/byronxu99/star-stack/blob/master/TestPhotos/dandelion_clipped.png?raw=true"/>
  &nbsp;
  <img src="https://github.com/byronxu99/star-stack/blob/master/TestPhotos/dandelion_clipped_m.png?raw=true"/>
  &nbsp;
  <img src="https://github.com/byronxu99/star-stack/blob/master/TestPhotos/dandelion_matte.png?raw=true"/>
</p>

<p align="center">
(a) Original image; (b) Image with "scribbles" to indicate areas as foreground and background; (c) Image matte produced by closed-form matting
</p>

An image with `N` pixels results in a `N x N` sparse system of equations, which can be quite large for images where `N` is several million pixels. To make this approach practical for large images, a recursive strategy is used. The input image is repeatedly scaled down in resolution until reaching a manageable size, where a low-resolution matte can be solved for. The matte is then upsampled and thresholded to leave only the most ambiguous values, making the resulting system of equations more sparse. The process is repeated until the full-resolution image matte is solved.

## Results
These are night sky photos taken on the MIT campus, where light pollution makes it difficult to see any stars. Individual exposures look like these:
<p align="center">
  <img src="https://github.com/byronxu99/star-stack/blob/master/Photos/ec/DSC_3659.png?raw=true" width="400"/>
  &nbsp;
  <img src="https://github.com/byronxu99/star-stack/blob/master/Photos/stata_2/DSC_3781.png?raw=true" width="400"/>
</p>

Stacking many exposures without alignment gives star trails:
<p align="center">
  <img src="https://github.com/byronxu99/star-stack/blob/master/Photos/ec/stack_unaligned.png?raw=true" width="400"/>
  &nbsp;
  <img src="https://github.com/byronxu99/star-stack/blob/master/Photos/stata_2/stack_unaligned.png?raw=true" width="400"/>
</p>

Closed-form matting is used to solve for image mattes. Note that matting isn't perfect; for example, the Walker Memorial radio antenna in the first image is not separated correctly from the background.
<p align="center">
  <img src="https://github.com/byronxu99/star-stack/blob/master/Photos/ec/DSC_3659_matte.png?raw=true" width="400"/>
  &nbsp;
  <img src="https://github.com/byronxu99/star-stack/blob/master/Photos/stata_2/DSC_3781_matte.png?raw=true" width="400"/>
</p>

After separating foreground and background, stars can be aligned and stacked independently. Merging the resulting stacked images gives reduced noise (difficult to see in the small thumbnails) without producing blurry stars:
<p align="center">
  <img src="https://github.com/byronxu99/star-stack/blob/master/Photos/ec/stack_aligned.png?raw=true" width="400"/>
  &nbsp;
  <img src="https://github.com/byronxu99/star-stack/blob/master/Photos/stata_2/stack_aligned_with_matte.png?raw=true" width="400"/>
</p>

Enhancing contrast in the stars before merging leads to a more impressive result. Stars are visible in the Boston night sky!
<p align="center">
  <img src="https://github.com/byronxu99/star-stack/blob/master/Photos/ec/stack_aligned_edited2.png?raw=true" width="400"/>
  &nbsp;
  <img src="https://github.com/byronxu99/star-stack/blob/master/Photos/stata_2/stack_aligned_with_matte_edited.png?raw=true" width="400"/>
</p>

## Usage
Edit `main.cpp` to run the code you want. Compile with CMake.

Future work: Build a proper command-line interface. Maybe also switch out the 6.815 PNG loader for an image library that handles more formats.

## More info
[Read my report (PDF) for more information](https://github.com/byronxu99/star-stack/blob/master/6.815_Final_Report.pdf).
