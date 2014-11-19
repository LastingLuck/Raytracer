#include "image.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/**
 * Image
 **/
Image::Image (int width_, int height_)
{
	bmpImg = new BMP();

    assert(width_ > 0);
    assert(height_ > 0);

    width           = width_;
    height          = height_;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;

	if (!bmpImg->SetSize(width, height)){
		printf("Error allocating image.");
		exit(-1);
	}
	
	pixels = new Pixel[num_pixels];
	int c = 0;
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			pixels[c].r = bmpImg->GetPixel(i,j).Red;
			pixels[c].g = bmpImg->GetPixel(i,j).Green;
			pixels[c].b = bmpImg->GetPixel(i,j).Blue;
			pixels[c].a = bmpImg->GetPixel(i,j).Alpha;
			c++;
		}
	}

    assert(pixels != NULL);
}


Image::Image (const Image& src)
{
    bmpImg = new BMP();
	
	width           = src.width;
    height          = src.height;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;

	RangedPixelToPixelCopy( *src.bmpImg, 0, width-1,
                         height-1 , 0, 
                         *bmpImg, 0,0 );

	if (!bmpImg->SetSize(width, height)){
		printf("Error allocating image.");
		exit(-1);
	}
	
	pixels = new Pixel[num_pixels];
	int c = 0;
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			pixels[c].r = bmpImg->GetPixel(i,j).Red;
			pixels[c].g = bmpImg->GetPixel(i,j).Green;
			pixels[c].b = bmpImg->GetPixel(i,j).Blue;
			pixels[c].a = bmpImg->GetPixel(i,j).Alpha;
			c++;
		}
	}

	assert(pixels != NULL);
}

Image::Image (char* fname){
	bmpImg = new BMP();

	if (!bmpImg->ReadFromFile(fname)){
		printf("Error loading image: %s", fname);
		exit(-1);
	}

	width = bmpImg->TellWidth();
	height = bmpImg->TellHeight();
	num_pixels = width * height;
	sampling_method = IMAGE_SAMPLING_POINT;

	pixels = new Pixel[num_pixels];
	int c = 0;
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			pixels[c].r = bmpImg->GetPixel(i,j).Red;
			pixels[c].g = bmpImg->GetPixel(i,j).Green;
			pixels[c].b = bmpImg->GetPixel(i,j).Blue;
			pixels[c].a = bmpImg->GetPixel(i,j).Alpha;
			c++;
		}
	}
}

Image::~Image (){
	delete bmpImg;
    pixels = NULL;
}

void Image::Write(char* fname){
	int c = 0;
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			RGBApixel p;
			p.Red = pixels[c].r;
			p.Green = pixels[c].g;
			p.Blue = pixels[c].b;
			p.Alpha = pixels[c].a;
			bmpImg->SetPixel(i,j,p);
			c++;
		}
	}
	bmpImg->SetBitDepth(24);
	bmpImg->WriteToFile(fname);
}

void Image::AddNoise (double factor)
{
	for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Pixel p = GetPixel(x, y);
            Pixel p_noise = p;
            double ran = ((double)rand() / RAND_MAX);
            if (factor > ran * 2) {
                Pixel randval = PixelRandom();
                p_noise.r += (randval.r * factor);
                p_noise.g += (randval.g * factor);
                p_noise.b += (randval.b * factor);
                p_noise.a += (randval.a * factor);
            }
            GetPixel(x, y) = p_noise;
        }
    }
}

void Image::Brighten (double factor)
{
	int x,y;
	for (x = 0 ; x < Width() ; x++)
	{
		for (y = 0 ; y < Height() ; y++)
		{
			Pixel p = GetPixel(x, y);
			Pixel scaled_p = p*factor;
			GetPixel(x,y) = scaled_p;
		}
	}
}


void Image::ChangeContrast (double factor)
{
    //Calculate the average luminance throughout the picture
	int ave_lum = 0;
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            ave_lum += GetPixel(x, y).Luminance();
        }
    }
    ave_lum /= num_pixels;
    //Interpolate/Extrapolate the pixels from the greyscale of the ave_lum
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Pixel p = GetPixel(x, y);
            Pixel sat_pixel = Pixel();
            sat_pixel.SetClamp(ave_lum+((p.r-ave_lum)*factor), ave_lum+((p.g-ave_lum)*factor), ave_lum+((p.b-ave_lum)*factor), ave_lum+((p.a-ave_lum)*factor));
            GetPixel(x, y) = sat_pixel;
        }
    }
}


void Image::ChangeSaturation(double factor)
{
	for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Pixel p = GetPixel(x, y);
            int lum = p.Luminance();
            Pixel sat_pixel = Pixel();
            sat_pixel.SetClamp(lum+((p.r-lum)*factor), lum+((p.g-lum)*factor), lum+((p.b-lum)*factor), lum+((p.a-lum)*factor));
            GetPixel(x, y) = sat_pixel;
        }
    }
}


Image* Image::Crop(int x, int y, int w, int h)
{
    if (ValidCoord(x, y)) {
        if (!ValidCoord(x+w, y)) {
            w = (w > 0) ? width - x : 0;
        }
        if (!ValidCoord(x, y+h)) {
            h = (h > 0) ? height - y : 0;
        }
        Image* im_crop = new Image(w, h);
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                im_crop->GetPixel(i, j) = GetPixel(i+x, j+y);
            }
        }
        return im_crop;
    }
	return NULL;
}


void Image::ExtractChannel(int channel)
{
	for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            //Does not specify which channel is which so I'm just going down the list
            //0 - red, 1 - green, 2 - blue, 3 - alpha
            Pixel p = GetPixel(x, y);
            switch (channel) {
                case 0:
                    p.g = 0;
                    p.b = 0;
                    p.a = 0;
                    break;
                case 1:
                    p.r = 0;
                    p.b = 0;
                    p.a = 0;
                    break;
                case 2:
                    p.r = 0;
                    p.g = 0;
                    p.a = 0;
                    break;
                case 3:
                    p.r = 0;
                    p.g = 0;
                    p.b = 0;
                    break;
                default:
                    break;
            }
            GetPixel(x, y) = p;
        }
    }
}

void Image::Quantize (int nbits)
{
	for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Pixel p = GetPixel(x, y);
            GetPixel(x, y) = PixelQuant(p, nbits);
        }
    }
}

void Image::RandomDither (int nbits)
{
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Component rand = ComponentRandom();
            Pixel p = GetPixel(x, y);
            Pixel o_dith = Pixel();
            o_dith.SetClamp(p.r+rand, p.g+rand, p.b+rand, p.a+rand);
            GetPixel(x, y) = PixelQuant(o_dith, nbits);
        }
    }
}


static int Bayer4[4][4] =
{
    {15,  7, 13,  5},
    { 3, 11,  1,  9},
    {12,  4, 14,  6},
    { 0,  8,  2, 10}
};


void Image::OrderedDither(int nbits)
{
	for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            int i = x % 4;
            int j = y % 4;
            Pixel p = GetPixel(x, y);
            int b_val = Bayer4[i][j];
            Pixel o_dith = Pixel();
            o_dith.SetClamp(p.r+b_val, p.g+b_val, p.b+b_val, p.a+b_val);
            GetPixel(x, y) = PixelQuant(o_dith, nbits);
        }
    }
}

/* Error-diffusion parameters */
const double
    ALPHA = 7.0 / 16.0,
    BETA  = 3.0 / 16.0,
    GAMMA = 5.0 / 16.0,
    DELTA = 1.0 / 16.0;

void Image::FloydSteinbergDither(int nbits)
{
	for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Pixel p = GetPixel(x, y);
            Pixel fs_dith = Pixel();
            fs_dith = PixelQuant(p, nbits);
            GetPixel(x, y) = fs_dith;
            //Pixel quant_err = Pixel(p.r-fs_dith.r, p.g-fs_dith.g, p.b-fs_dith.b, p.a-fs_dith.a);
            Pixel quant_err = Pixel();
            quant_err.SetClamp(p.r-fs_dith.r, p.g-fs_dith.g, p.b-fs_dith.b, p.a-fs_dith.a);
            int x1 = x+1, x2 = x-1, y1 = y+1;
            if (x1 >= width) {
                x1 = width - 1;
            }
            if (x2 < 0) {
                x2 = 0;
            }
            if (y1 >= height) {
                y1 = height - 1;
            }
            //Alpha
            Pixel alpha = GetPixel(x1, y);
            alpha = alpha + (quant_err * ALPHA);
            GetPixel(x1, y) = alpha;
            //Beta
            Pixel beta = GetPixel(x2, y1);
            beta = beta + (quant_err * BETA);
            GetPixel(x2, y1) = beta;
            //Gamma
            Pixel gamma = GetPixel(x, y1);
            gamma = gamma + (quant_err * GAMMA);
            GetPixel(x, y1) = gamma;
            //Delta
            Pixel delta = GetPixel(x1, y1);
            delta = delta + (quant_err * DELTA);
            GetPixel(x1, y1) = delta;
        }
    }
}

void Image::Blur(int n)
{
    //Not sure how to handle even values of n
	//Blurs with the gaussian in the x direction first, then the y direction
    //Make the kernel
    int size = n * 2 - 1;
    double sigma = n;
    double sum = 0.0;
    int center = ceil(size / 2.0);
    double kernal[size];
    double div = sqrt(2 * M_PI) * sigma;
    for (int i = 1; i <= size; i++) {
        double dist = ((i-center) * (i-center));
        //printf("Dist: %d\n", i-center);
        double gauss = exp(-dist/(2*sigma*sigma)) / div;
        //printf("%f\n", exp(-dist/(2*sigma*sigma)));
        //printf("%f\n", sqrt(2*M_PI)*sigma);
        kernal[i-1] = gauss;
        sum += gauss;
    }
    //Normalize
    for (int i = 0; i < size; i++) {
        kernal[i] /= sum;
    }
    //printf("\nKernal Size %d\n", size);
    //printf("Center: %d\n", center);
    //for (int i = 0; i < size; i++) {
    //    printf("%f, ", kernal[i]);
    //}
    //printf("\n");
    //X
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            Pixel sum = Pixel();
            for (int i = 0; i < size; i++) {
                int xindex = i - (center - 1) + x;
                //Extend pixel past border if applicable
                if (xindex < 0) {
                    xindex = 0;
                }
                if (xindex >= width) {
                    xindex = width - 1;
                }
                Pixel p = GetPixel(xindex, y);
                sum = sum + (p * kernal[i]);
            }
            GetPixel(x, y) = sum;
        }
    }
    //Y
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Pixel sum = Pixel();
            for (int i = 0; i < size; i++) {
                int yindex = i - (center - 1) + y;
                //Extend pixel past border if applicable
                if (yindex < 0) {
                    yindex = 0;
                }
                if (yindex >= height) {
                    yindex = height - 1;
                }
                Pixel p = GetPixel(x, yindex);
                sum = sum + (p * kernal[i]);
            }
            GetPixel(x, y) = sum;
        }
    }
}

void Image::Sharpen(int n)
{
    //Can't figure out how to make it work from the blur
	Image blur = Image(*this);
    blur.Blur(n);
    //Extrapolate from the blurred image
    Pixel sharp, p;
    double amt = .5;
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            p = blur.GetPixel(x, y);
            sharp = GetPixel(x, y);
            //sharp.SetClamp(sharp.r+(p.r-sharp.r)*.2, sharp.g+(p.g-sharp.g)*.2, sharp.b+(p.b-sharp.b)*.2, sharp.a+(p.a-sharp.a)*.2);
            //sharp.SetClamp(sharp.r+.2*(sharp.r-p.r), sharp.g+.2*(sharp.g-p.g), sharp.b+.2*(sharp.b-p.b), sharp.a+.2*(sharp.a-p.a));
            sharp.SetClamp((1.0+amt)*sharp.r-amt*p.r, (1.0+amt)*sharp.g-amt*p.g, (1.0+amt)*sharp.b-amt*p.b, (1.0+amt)*sharp.a-amt*p.a);
            GetPixel(x, y) = sharp;
        }
    }
}

void Image::EdgeDetect()
{
	/*                                  -1 -1 -1
    * Multiply each pixel by the matrix -1  8 -1
    *                                   -1 -1 -1
    */
    Image dest = Image(*this);
    int kernal[][3] = {{-1, -1, -1}, {-1, 8, -1}, {-1, -1, -1}};
    Pixel p;
    for (int x = 1; x < width-1; x++) {
        for (int y = 1; y < height-1; y++) {;
            int sumr = 0, sumg = 0, sumb = 0, suma = 0;
            for (int i = -1; i < 2; i++) {
                for (int j = -1; j < 2; j++) {
                    p = GetPixel(x+i, y+j);
                    sumr += p.r * kernal[i+1][j+1];
                    sumg += p.g * kernal[i+1][j+1];
                    sumb += p.b * kernal[i+1][j+1];
                    suma += p.a * kernal[i+1][j+1];
                }
            }
            p.SetClamp(sumr, sumg, sumb, suma);
            dest.GetPixel(x, y) = p;
        }
    }
    delete pixels;
    pixels = new Pixel[num_pixels];
	int c = 0;
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			pixels[c].r = dest.GetPixel(i,j).r;
			pixels[c].g = dest.GetPixel(i,j).g;
			pixels[c].b = dest.GetPixel(i,j).b;
			pixels[c].a = dest.GetPixel(i,j).a;
			c++;
		}
	}
}

Image* Image::Scale(double sx, double sy)
{
    int sWidth = ceil(width*sx), sHeight = ceil(height*sy);
    //printf("%d / %d\n", height, sy);
    Image* scale = new Image(sWidth, sHeight);
    for (int x = 0; x < sWidth; x++) {
        for (int y = 0; y < sHeight; y++) {
            float u = x / sx;
            float v = y / sy;
            //printf("X: %d\nY: %d\n", x, y);
            scale->GetPixel(x, y) = Sample(u, v);
        }
    }
	return scale;
}

Image* Image::Rotate(double angle)
{
    //Not sure if the angle is in degrees or radians
    //Going for degrees since it's easier to check. (Negative theta because of reverse mapping)
    angle *= (M_PI / 180.0);
    //Assuming image is rotated around the center, find new width/height
    double center_x = width / 2.0;
    double center_y = height / 2.0;
    int rWidth, rHeight;
    //TR, TL, BL, BR
    angle = -angle;
    double c1_x = (width-center_x)*cos(angle) - (height-center_y)*sin(angle);
    double c1_y = (width-center_x)*sin(angle) + (height-center_y)*cos(angle);
    double c2_x = (-center_x)*cos(angle) - (height-center_y)*sin(angle);
    double c2_y = (-center_x)*sin(angle) + (height-center_y)*cos(angle);
    double c3_x = (-center_x)*cos(angle) - (-center_y)*sin(angle);
    double c3_y = (-center_x)*sin(angle) + (-center_y)*cos(angle);
    double c4_x = (width-center_x)*cos(angle) - (-center_y)*sin(angle);
    double c4_y = (width-center_x)*sin(angle) + (-center_y)*cos(angle);
    c1_x += center_x;
    c1_y += center_y;
    c2_x += center_x;
    c2_y += center_y;
    c3_x += center_x;
    c3_y += center_y;
    c4_x += center_x;
    c4_y += center_y;
    //Test each value to find the max
    //X
    double val1, val2, val3, val4;
    if (c1_x > c2_x) {
        val1 = c1_x; //Max Candidate
        val2 = c2_x; //Min Candidate
    }
    else {
        val1 = c2_x;
        val2 = c1_x;
    }
    if (c3_x > c4_x) {
        val3 = c3_x;
        val4 = c4_x;
    }
    else {
        val3 = c4_x;
        val4 = c3_x;
    }
    double max = (val1 > val3) ? ceil(val1) : ceil(val3);
    double min = (val2 < val4) ? floor(val2) : floor(val4);
    rWidth = ceil(max - min);
    //Y
    if (c1_y > c2_y) {
        val1 = c1_y; //Max Candidate
        val2 = c2_y; //Min Candidate
    }
    else {
        val1 = c2_y;
        val2 = c1_y;
    }
    if (c3_y > c4_y) {
        val3 = c3_y;
        val4 = c4_y;
    }
    else {
        val3 = c4_y;
        val4 = c3_y;
    }
    max = (val1 > val3) ? ceil(val1) : ceil(val3);
    min = (val2 < val4) ? floor(val2) : floor(val4);
    rHeight = ceil(max - min);
    //Rotate
    Image* rotate = new Image(rWidth, rHeight);
    double rcenter_x = rWidth / 2.0;
    double rcenter_y = rHeight / 2.0;
    double cent_diff_x = center_x - rcenter_x;
    double cent_diff_y = center_y - rcenter_y;
    for (int x = 0; x < rWidth; x++) {
        for (int y = 0; y < rHeight; y++) {
            //calculate rotation from local space of (0,0) in center of picture
            float u = rcenter_x + ((x - rcenter_x) * cos(-angle) - (y - rcenter_y) * sin(-angle));
            float v = rcenter_y + ((x - rcenter_x) * sin(-angle) + (y - rcenter_y) * cos(-angle));
            u += cent_diff_x;
            v += cent_diff_y;
            if(u < 0.0 || u > width || v < 0.0 || v > height) {
                rotate->GetPixel(x, y) = Pixel(255, 255, 255, 255);
            }
            else {
                rotate->GetPixel(x, y) = Sample(u, v);
            }
        }
    }
	return rotate;
}

void Image::Fun()
{
    //I have no idea what's going on in this filter, but it looks cool
    double theta = 45.0;
    Image idk = Image(width, height);
    theta *= (M_PI / 180.0);
    double center_x = width / 2.0;
    double center_y = height / 2.0;
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            double dist_x = sqrt((x-center_x) * (x-center_x));
            double dist_y = sqrt((y-center_y) * (y-center_y));
            double u = center_x + ((x - center_x) * cos(theta*dist_x) - (y - center_y) * sin(theta*dist_x));
            double v = center_y + ((x - center_x) * sin(theta*dist_y) + (y - center_y) * cos(theta*dist_y));
            idk.GetPixel(x, y) = Sample(u, v);
        }
    }
    delete pixels;
    pixels = new Pixel[num_pixels];
	int c = 0;
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			pixels[c].r = idk.GetPixel(i,j).r;
			pixels[c].g = idk.GetPixel(i,j).g;
			pixels[c].b = idk.GetPixel(i,j).b;
			pixels[c].a = idk.GetPixel(i,j).a;
			c++;
		}
	}
}

/**
 * Image Sample
 **/
void Image::SetSamplingMethod(int method)
{
    assert((method >= 0) && (method < IMAGE_N_SAMPLING_METHODS));
    sampling_method = method;
}


Pixel Image::Sample (double u, double v){
    switch (sampling_method) {
        case IMAGE_SAMPLING_POINT: {
            int x = round(u), y = round(v);
            x = (x < 0) ? 0 : (x >= width) ? width - 1 : x;
            y = (y < 0) ? 0 : (y >= height) ? height - 1 : y;
            return GetPixel(x, y);
        }
        case IMAGE_SAMPLING_BILINEAR: {
            int x1 = floor(u), x2 = ceil(u);
            int y1 = floor(v), y2 = ceil(v);
            if (x2 >= width) {
                x2 = width - 1;
            }
            if (y2 >= height) {
                y2 = height - 1;
            }
            double rval1, rval2;
            if (x2 - x1 == 0) {
                rval1 = 1.0;
                rval2 = 0.0;
            }
            else {
                rval1 = (x2 - u) / (double)(x2 - x1);
                rval2 = (u - x1) / (double)(x2 - x1);
            }
            Pixel q11 = GetPixel(x1, y1);
            Pixel q12 = GetPixel(x1, y2);
            Pixel q21 = GetPixel(x2, y1);
            Pixel q22 = GetPixel(x2, y2);
            Pixel r1 = Pixel(), r2 = Pixel();
            r1.SetClamp(rval1*q11.r+rval2*q21.r, rval1*q11.g+rval2*q21.g, rval1*q11.b+rval2*q21.b, rval1*q11.a+rval2*q21.a);
            r2.SetClamp(rval1*q12.r+rval2*q22.r, rval1*q12.g+rval2*q22.g, rval1*q12.b+rval2*q22.b, rval1*q12.a+rval2*q22.a);
            double pval1, pval2;
            if (y2 - y1 == 0) {
                pval1 = 1.0;
                pval2 = 0.0;
            }
            else {
                pval1 = (y2 - v) / (double)(y2 - y1);
                pval2 = (v - y1) / (double)(y2 - y1);
            }
            Pixel p = Pixel();
            p.SetClamp(pval1*r1.r+pval2*r2.r, pval1*r1.g+pval2*r2.g, pval1*r1.b+pval2*r2.b, pval1*r1.a+pval2*r2.a);
            return p;
        }
        case IMAGE_SAMPLING_GAUSSIAN: {
            //Make the kernel
            //double sigma = n/2;
            double sum = 0.0;
            double div = 4.0 * M_PI;
            double kernal[4][4];
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    int x = floor(u) + (i - 1);
                    int y = floor(v) + (j - 1);
                    double dist = ((u-x)*(u-x)) + ((v-y)*(v-y));
                    double gauss = exp(-dist/4.0) / div;
                    kernal[i][j] = gauss;
                    sum += gauss;
                }
            }
            //Normalize
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    kernal[i][j] /= sum;
                }
            }
            Pixel gsum = Pixel();
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    int xindex = floor(u) + (i - 1);
                    if (xindex < 0) {
                        xindex = 0;
                    }
                    if (xindex >= width) {
                        xindex = width - 1;
                    }
                    int yindex = floor(v) + (j - 1);
                    if (yindex < 0) {
                        yindex = 0;
                    }
                    if (yindex >= height) {
                        yindex = height - 1;
                    }
                    Pixel p = GetPixel(xindex, yindex);
                    gsum = gsum + (p * kernal[i][j]);
                }
            }
            return gsum;
        }
    }
	return Pixel();
}