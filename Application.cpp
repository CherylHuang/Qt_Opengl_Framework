#include "Application.h"
#include "qt_opengl_framework.h"
#include <vector>
#include <stdlib.h>

float* Sorting(float ary[], int size);
int FindMax(int ary[]);
int BilinearInterpolation(int Q11, int Q12, int Q21, int Q22, int x1, int x2, int y1, int y2, double x, double y);

Application::Application()
{

}
Application::~Application()
{

}
//****************************************************************************
//
// * 初始畫面，並顯示Ntust.png圖檔
// 
//============================================================================
void Application::createScene( void )
{
	
	ui_instance = Qt_Opengl_Framework::getInstance();
	
}

//****************************************************************************
//
// * 打開指定圖檔
// 
//============================================================================
void Application::openImage( QString filePath )
{
	mImageSrc.load(filePath);
	mImageDst.load(filePath);

	renew();

	img_data = mImageSrc.bits();
	img_width = mImageSrc.width();
	img_height = mImageSrc.height();

	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
}
//****************************************************************************
//
// * 刷新畫面
// 
//============================================================================
void Application::renew()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageDst));

	std::cout << "Renew" << std::endl;
}

//****************************************************************************
//
// * 畫面初始化
// 
//============================================================================
void Application::reload()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageSrc));
}

//****************************************************************************
//
// * 儲存圖檔
// 
//============================================================================
void Application::saveImage(QString filePath )
{
	mImageDst.save(filePath);
}

//****************************************************************************
//
// * 將圖檔資料轉換為RGB色彩資料
// 
//============================================================================
unsigned char* Application::To_RGB( void )
{
	unsigned char *rgb = new unsigned char[img_width * img_height * 3];
	int i, j;

	if (! img_data )
		return NULL;

	// Divide out the alpha
	for (i = 0; i < img_height; i++)
	{
		int in_offset = i * img_width * 4;
		int out_offset = i * img_width * 3;

		for (j = 0 ; j < img_width ; j++)
		{
			RGBA_To_RGB(img_data + (in_offset + j*4), rgb + (out_offset + j*3));
		}
	}

	return rgb;
}

void Application::RGBA_To_RGB( unsigned char *rgba, unsigned char *rgb )
{
	const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

	unsigned char  alpha = rgba[3];

	if (alpha == 0)
	{
		rgb[0] = BACKGROUND[0];
		rgb[1] = BACKGROUND[1];
		rgb[2] = BACKGROUND[2];
	}
	else
	{
		float	alpha_scale = (float)255 / (float)alpha;
		int	val;
		int	i;

		for (i = 0 ; i < 3 ; i++)
		{
			val = (int)floor(rgba[i] * alpha_scale);
			if (val < 0)
				rgb[i] = 0;
			else if (val > 255)
				rgb[i] = 255;
			else
				rgb[i] = val;
		}
	}
}
//------------------------Color------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Gray()
{
	unsigned char *rgb = To_RGB();

	for (int i=0; i<img_height; i++)
	{
		for (int j=0; j<img_width; j++)
		{
			int offset_rgb = i*img_width*3+j*3;
			int offset_rgba = i*img_width*4+j*4;
			unsigned char gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			for (int k=0; k<3; k++)
				img_data[offset_rgba+k] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Uniform()
{
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			img_data[offset_rgba + rr] = rgb[offset_rgb + rr] / 32 * 32;
			img_data[offset_rgba + gg] = rgb[offset_rgb + gg] / 32 * 32;
			img_data[offset_rgba + bb] = rgb[offset_rgb + bb] / 64 * 64;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Populosity()
{
	unsigned char *rgb = this->To_RGB();
	int ColorCount[32768] = { 0 };
	int PopularColor[256] = { 0 };

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int r, g, b, color;
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			r = (rgb[offset_rgb + rr] >> 3) << 3;	// convert to 5 bit color
			g = (rgb[offset_rgb + gg] >> 3) << 3;
			b = (rgb[offset_rgb + bb] >> 3) << 3;
			color = ((r >> 3) << 5 * rr) + ((g >> 3) << 5 * gg) + ((b >> 3) << 5 * bb);
			ColorCount[color]++;
		}
	}
	for (int i = 0; i < 256; i++)
	{
		int color = FindMax(ColorCount);
		PopularColor[i] = color;
		ColorCount[color] = 0;	// clear
	}
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int r, g, b, minR, minG, minB;
			minR = minG = minB = 0;
			int MinColorDist = 999999999;
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			for (int k = 0; k < 256; k++)
			{
				//r = (PopularColor[k] >> 5 * rr) << 3;
				//g = (PopularColor[k] >> 5 * gg) << 3;
				//b = (PopularColor[k] >> 5 * bb) << 3;
				r = ((PopularColor[k] & (0b11111 << 5 * rr)) >> 5 * rr) << 3;
				g = ((PopularColor[k] & (0b11111 << 5 * gg)) >> 5 * gg) << 3;
				b = ((PopularColor[k] & (0b11111 << 5 * bb)) >> 5 * bb) << 3;
				int ColorDist[3] = { r - rgb[offset_rgb + rr], g - rgb[offset_rgb + gg], b - rgb[offset_rgb + bb] };
				int dist = ColorDist[0] * ColorDist[0] + ColorDist[1] * ColorDist[1] + ColorDist[2] * ColorDist[2];
				if (MinColorDist > dist)
				{
					MinColorDist = dist;
					minR = r; minG = g; minB = b;
				}
			}
			img_data[offset_rgba + rr] = minR;
			img_data[offset_rgba + gg] = minG;
			img_data[offset_rgba + bb] = minB;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Dithering------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Threshold()
{
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			float gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];
			float range = gray / WHITE;
			for (int k = 0; k < 3; k++) {
				if (range < 0.5)img_data[offset_rgba + k] = BLACK;
				else img_data[offset_rgba + k] = WHITE;
			}
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Bright()
{
	unsigned char *rgb = this->To_RGB();
	int count = 0;
	float sum, threshold, average;
	sum = threshold = average = 0;
	float *arySort;
	float **aryGray;
	arySort = new float[img_height * img_width];
	aryGray = new float*[img_height];

	for (int i = 0; i < img_height; i++)
	{
		aryGray[i] = new float[img_width];
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			unsigned char gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];
			float range = (float)gray / WHITE;
			aryGray[i][j] = range;
			arySort[count] = range;
			count++;
			sum += range;
		}
	}
	average = sum / (img_height * img_width);
	count = count * (1 - average);
	arySort = Sorting(arySort, img_height * img_width);
	threshold = arySort[count];
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i * img_width * 4 + j * 4;
			for (int k = 0; k < 3; k++) {
				if (aryGray[i][j] < threshold)img_data[offset_rgba + k] = BLACK;
				//if (aryGray[i][j] < average)img_data[offset_rgba + k] = BLACK;
				else img_data[offset_rgba + k] = WHITE;
			}
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	delete[] arySort;
	for (int i = 0; i < img_height; i++) delete[] aryGray[i];
	delete[] aryGray;
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Random()
{
	srand((unsigned)time(NULL));

	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			float gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];
			float range = gray / WHITE;
			range += -0.2f + (rand() % 5 / 10.f);
			for (int k = 0; k < 3; k++) {
				if (range < 0.5)img_data[offset_rgba + k] = BLACK;
				else img_data[offset_rgba + k] = WHITE;
			}
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Cluster()
{
	unsigned char *rgb = this->To_RGB();

	float ThresholdMatrix[4][4] = {
		{0.7059, 0.3529, 0.5882, 0.2353 },
		{0.0588, 0.9412, 0.8235, 0.4118 },
		{0.4706, 0.7647, 0.8824, 0.1176 },
		{0.1765, 0.5294, 0.2941, 0.6471 },
	};

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			float gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];
			float range = gray / WHITE;
			float threshold = ThresholdMatrix[i % 4][j % 4];
			for (int k = 0; k < 3; k++) {
				if (range < threshold)img_data[offset_rgba + k] = BLACK;
				else img_data[offset_rgba + k] = WHITE;
			}
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_FS()
{
	unsigned char *rgb = this->To_RGB();
	int offset_rgb, offset_rgba;
	bool flag = true;
	float **pixel;
	pixel = new float*[img_height];

	for (int i = 0; i < img_height; i++)
	{
		pixel[i] = new float[img_width];
		for (int j = 0; j < img_width; j++)
		{
			offset_rgb = i * img_width * 3 + j * 3;
			float gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];
			pixel[i][j] = gray;
		}
	}

	for (int i = 0; i < img_height; i++)
	{
		for (int j = flag ? 0 : img_width - 1; j < img_width && j >= 0; j = flag ? j+1 : j-1)
		{
			float oldpixel = pixel[i][j];
			float newpixel = (oldpixel / WHITE) < 0.5f ? BLACK : WHITE;
			float quant_error = oldpixel - newpixel;
			pixel[i][j] = newpixel;
			if (j + 1 < img_width)
				pixel[i][j + 1] += 7.0f / 16.0f * quant_error;
			if (i + 1 < img_height && j - 1 >= 0)
				pixel[i + 1][j - 1] += 3.0f / 16.0f * quant_error;
			if (i + 1 < img_height)
				pixel[i + 1][j] += 5.0f / 16.0f * quant_error;
			if (i + 1 < img_height && j + 1 < img_width)
				pixel[i + 1][j + 1] += 1.0f / 16.0f * quant_error;
		}
		flag = flag ? false : true;
	}

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i * img_width * 4 + j * 4;
			for (int k = 0; k < 3; k++) {
				img_data[offset_rgba + k] = pixel[i][j];
			}
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	for (int i = 0; i < img_height; i++) delete[] pixel[i];
	delete[] pixel;
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
float Application::FindClosestColor(float pixel, int color)
{
	float RG[8] = { 0, 36, 73, 109, 146, 182, 219, 255 };
	float B[4] = { 0, 85, 170, 255 };
	if (color == rr || color == gg) return RG[(int)(pixel / (WHITE / 7.0f))];
	else if (color == bb) return B[(int)(pixel / (WHITE / 3.0f))];
}
void Application::Dither_Color()
{
	unsigned char *rgb = this->To_RGB();
	bool flag = true;
	float **pixel;
	pixel = new float*[img_height];

	for (int i = 0; i < img_height; i++)
	{
		pixel[i] = new float[img_width * 3];
		for (int j = 0; j < img_width; j+=3)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			for (int k = 0; k < 3; k++) pixel[i][j + k] = rgb[offset_rgb + k];
		}
	}

	for (int i = 0; i < img_height; i++)
	{
		for (int j = flag ? 0 : img_width - 3; j < img_width - 3 && j >= 0; j = flag ? j + 3 : j - 3)
		{
			for (int k = 0; k < 3; k++) {
				float oldpixel = pixel[i][j + k];
				float newpixel = FindClosestColor(oldpixel, k);
				float quant_error = oldpixel - newpixel;
				pixel[i][j + k] = newpixel;
				if (j + 3 < img_width)
					pixel[i][j + 3 + k] += 7.0f / 16.0f * quant_error;
				if (i + 1 < img_height && j - 3 >= 0)
					pixel[i + 1][j - 3 + k] += 3.0f / 16.0f * quant_error;
				if (i + 1 < img_height)
					pixel[i + 1][j + k] += 5.0f / 16.0f * quant_error;
				if (i + 1 < img_height && j + 3 < img_width)
					pixel[i + 1][j + 3 + k] += 1.0f / 16.0f * quant_error;
			}
		}
		flag = flag ? false : true;
	}

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j+=3)
		{
			int offset_rgba = i * img_width * 4 + j * 4;
			for (int k = 0; k < 3; k++) {
				img_data[offset_rgba + k] = pixel[i][j + k];
			}
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	for (int i = 0; i < img_height; i++) delete[] pixel[i];
	delete[] pixel;

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Filter------------------------

///////////////////////////////////////////////////////////////////////////////
//
//     Filtering the img_data array by the filter from the parameters
//
///////////////////////////////////////////////////////////////////////////////
int Application::filtering(unsigned char * rgb, int startIndex, float **filter, int width, int height)
{
	int Sum = 0, weight = 0;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int offset_rgb = startIndex + i * this->img_width * 3 + j * 3;
			if (offset_rgb < 0 || offset_rgb > this->img_height * this->img_width * 3) {
				Sum += 0;
			}
			else {
				Sum += rgb[offset_rgb] * filter[i][j];
			}
			weight += filter[i][j];
		}
	}
	return (int)(Sum / weight);
}
void Application::Filter(float FilterData[])
{
	unsigned char *rgb = this->To_RGB();
	int FWidth = 5, FHeight = 5;
	float **FilterMatrix;
	FilterMatrix = new float*[FHeight];
	for (int i = 0; i < FHeight; i++)
	{
		FilterMatrix[i] = new float[FWidth];
		for (int j = 0; j < FWidth; j++)
		{
			FilterMatrix[i][j] = FilterData[i * FWidth + FWidth];
		}
	}
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			int startInx = offset_rgb - 2 * img_width * 3 + 2 * 3;
			for (int k = 0; k < 3; k++) {
				img_data[offset_rgba + k] = filtering(rgb, startInx + k, FilterMatrix, FWidth, FHeight);
			}
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	for (int i = 0; i < FHeight; i++) delete[] FilterMatrix[i];
	delete[] FilterMatrix;
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
void Application::Filter(float **FilterData, int n)
{
	unsigned char *rgb = this->To_RGB();
	int FWidth = n, FHeight = n;
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			int startInx = offset_rgb - ((int)(n / 2)) * img_width * 3 - ((int)(n / 2)) * 3;
			for (int k = 0; k < 3; k++) {
				img_data[offset_rgba + k] = filtering(rgb, startInx + k, FilterData, FWidth, FHeight);
			}
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Box()
{
	float FilterData[25] = {
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1
	};
	this->Filter(FilterData);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Bartlett()
{
	float FilterData[25] = {
		1, 2, 3, 2, 1,
		2, 4, 6, 4, 2,
		3, 6, 9, 6, 3,
		2, 4, 6, 4, 2,
		1, 2, 3, 2, 1
	};
	this->Filter(FilterData);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian()
{
	//float FilterData[25] = {
	//	1, 4,  6,  4,  1,
	//	4, 16, 24, 16, 4,
	//	6, 24, 36, 24, 6,
	//	4, 16, 24, 16, 4,
	//	1, 4,  6,  4,  1,
	//};
	//this->Filter(FilterData);
	this->Filter_Gaussian_N(5);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian_N( int N )
{
	float **iPascal, **GaussianMatrix;
	float *GaussianAry;
	GaussianAry = new float[N];
	iPascal = new float*[N];
	GaussianMatrix = new float*[N];
	for (int i = 0; i < N; i++) {
		iPascal[i] = new float[N];
		GaussianMatrix[i] = new float[N];
	}
	for (int i = 0; i < N; i++) { //Pascal init
		iPascal[0][i] = 1;
		iPascal[i][0] = 1;
	}
	for (int i = 1; i < N - 1; i++) { //complete Pascal
		for (int j = 1; j < N - 1; j++) {
			iPascal[i][j] = iPascal[i][j - 1] + iPascal[i - 1][j];
		}
	}
	for (int i = 0; i < N; i++) GaussianAry[i] = iPascal[N - 1 - i][i];	//get the first line of GaussianMatrix
	for (int i = 0; i < N; i++) {
		GaussianMatrix[0][i] = GaussianAry[i];
		for (int j = 1; j < N; j++) {
			if (i > 0)
				GaussianMatrix[i][j] = GaussianAry[j] * GaussianMatrix[0][i];
		}
	}

	this->Filter(GaussianMatrix, N);	//call filter

	delete[] GaussianAry;
	for (int i = 0; i < N; i++) {
		delete[] iPascal[i];
		delete[] GaussianMatrix[i];
	}
	delete[] iPascal;
	delete[] GaussianMatrix;
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Edge()
{

}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Enhance()
{
	unsigned char *rgb = this->To_RGB();



	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Size------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Half_Size()
{
	this->Resize(0.5f);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Double_Size()
{
	this->Resize(2.0f);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Scale the image dimensions by the given factor.  The given factor is 
//	assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Resize( float scale )
{
	unsigned char * rgb = To_RGB();

	int newHeight = scale * img_height;
	int newWidth = scale * img_width;
	img_data = new unsigned char[newWidth * newHeight * 4];

	for (int i = 0; i < newHeight; i++) {
		for (int j = 0; j < newWidth; j++) {
			float x = (float)img_height / newHeight * i;
			float y = (float)img_width / newWidth * j;
			int x1, x2, y1, y2;
			x1 = (int)x;
			x2 = (x1 >= img_height) ? x1 : (x1 + 1);
			y1 = (int)y;
			y2 = (y1 >= img_width) ? y1 : (y1 + 1);
			int Q11, Q12, Q21, Q22;
			Q11 = x1 * img_width * 3 + y1 * 3;
			Q12 = x1 * img_width * 3 + y2 * 3;
			Q21 = x2 * img_width * 3 + y1 * 3;
			Q22 = x2 * img_width * 3 + y2 * 3;
			int offset_rgba = i * newWidth * 4 + j * 4;
			img_data[offset_rgba + rr] = BilinearInterpolation(rgb[Q11 + rr], rgb[Q12 + rr], rgb[Q21 + rr], rgb[Q22 + rr], x1, x2, y1, y2, x, y);
			img_data[offset_rgba + gg] = BilinearInterpolation(rgb[Q11 + gg], rgb[Q12 + gg], rgb[Q21 + gg], rgb[Q22 + gg], x1, x2, y1, y2, x, y);
			img_data[offset_rgba + bb] = BilinearInterpolation(rgb[Q11 + bb], rgb[Q12 + bb], rgb[Q21 + bb], rgb[Q22 + bb], x1, x2, y1, y2, x, y);
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	this->img_height = newHeight;
	this->img_width = newWidth;
	this->mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
int BilinearInterpolation(int Q11, int Q12, int Q21, int Q22, int x1, int x2, int y1, int y2, double x, double y)
{
	float a0 = (float)x - x1;
	float b0 = (float)y - y1;
	float a1 = (float)x2 - x1;
	float b1 = (float)y2 - y1;
	return (int)((a1 - a0)*(b1 - b0)*Q11 + a0*(b1 - b0)*Q21 + b0*(a1 - a0)*Q12 + a0*b0*Q22) / (a1 * b1);
}
///////////////////////////////////////////////////////////////////////////////
//
//  resample_src for resize and rotate
//
///////////////////////////////////////////////////////////////////////////////
void Application::resample_src(int u, int v, float ww, unsigned char* rgba)
{

}

//////////////////////////////////////////////////////////////////////////////
//
//  Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Rotate( float angleDegrees )
{

	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Composing------------------------


void Application::loadSecondaryImge( QString filePath )
{
	mImageSrcSecond.load(filePath);

	renew();

	img_data2 = mImageSrcSecond.bits();
	img_width2 = mImageSrcSecond.width();
	img_height2 = mImageSrcSecond.height();
}

//////////////////////////////////////////////////////////////////////////
//
//	Composite the image A and image B by Over, In, Out, Xor and Atom. 
//
//////////////////////////////////////////////////////////////////////////
void Application::Comp_image( int tMethod )
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Over()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_In()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Out()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Atop()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Xor()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

//------------------------NPR------------------------

///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::NPR_Paint()
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

void Application::NPR_Paint_Layer( unsigned char *tCanvas, unsigned char *tReferenceImage, int tBrushSize )
{

}

///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void Application::Paint_Stroke( const Stroke& s )
{
	int radius_squared = (int)s.radius * (int)s.radius;
	for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) 
	{
		for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) 
		{
			int x_loc = (int)s.x + x_off;
			int y_loc = (int)s.y + y_off;

			// are we inside the circle, and inside the image?
			if ((x_loc >= 0 && x_loc < img_width && y_loc >= 0 && y_loc < img_height)) 
			{
				int dist_squared = x_off * x_off + y_off * y_off;
				int offset_rgba = (y_loc * img_width + x_loc) * 4;

				if (dist_squared <= radius_squared) 
				{
					img_data[offset_rgba + rr] = s.r;
					img_data[offset_rgba + gg] = s.g;
					img_data[offset_rgba + bb] = s.b;
					img_data[offset_rgba + aa] = s.a;
				} 
				else if (dist_squared == radius_squared + 1) 
				{
					img_data[offset_rgba + rr] = (img_data[offset_rgba + rr] + s.r) / 2;
					img_data[offset_rgba + gg] = (img_data[offset_rgba + gg] + s.g) / 2;
					img_data[offset_rgba + bb] = (img_data[offset_rgba + bb] + s.b) / 2;
					img_data[offset_rgba + aa] = (img_data[offset_rgba + aa] + s.a) / 2;
				}
			}
		}
	}
}





///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
	unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}


float* Sorting(float ary[], int size)
{
	int iSwitch;
	float temp;
	do {
		iSwitch = 0;
		for (int i = 0; i < size - 1; i++) {
			if (ary[i] > ary[i + 1]) {
				temp = ary[i]; ary[i] = ary[i + 1]; ary[i + 1] = temp;
				iSwitch = 1;
			}
		}
	} while (iSwitch);
	return ary;
}

int FindMax(int ary[])
{
	int size, max = 0;
	size = 32 * 32 * 32;
	for (int i = 0; i < size; i++) {
		if (ary[max] < ary[i]) max = i;
	}
	return max;
}