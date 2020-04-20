#include <string>
#include <vector>
#include <iostream>
#include <Magick++.h>
#include "imageio.h"

using namespace std;
using namespace Magick;

//
//	 imageio_init() - 	// Any initialisation that may be required
//
void imageio_init()
{
	InitializeMagick("C:\Program Files\ImageMagick-7.0.8-Q16");
}

//
//	readImage(std::string, int *, int *) -  Reads an image from disk, returns data as a float array
//
float *readImage(std::string fileName, int *w, int *h)
{
	Image mainImage(fileName);
	Pixels mainPixelCache(mainImage);
	*w = mainImage.columns();
	*h = mainImage.rows();
	Quantum *pixels = mainPixelCache.get(0, 0, *w, *h);

	float *buf = (float *)(malloc(*w * *h * 3 * sizeof(float)));
	for (int i = 0; i < *w; i++)
	{
		for (int j = 0; j < *h; j++)
		{
			for (int c = 0; c < 3; c++) buf[(j * *w + i) * 3 + c] = (float)(pixels[(j * *w + i) * 3 + c]) / (float)(QuantumRange);
		}
	}

	return buf;
}

//
//	writeImage(float *, std::string, int, int) - Writes an image in float array format to disk
//
void writeImage(float *buf, std::string fileName, int w, int h)
{
	Image outImage(Geometry(w, h), Color(QuantumRange, QuantumRange, QuantumRange, 0));
	Pixels outPixelCache(outImage);
	Quantum *outPixels = outPixelCache.get(0, 0, w, h);

	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			for (int c = 0; c < 3; c++)  outPixels[((h - 1 - j) * w + i) * 4 + c] = buf[(j * w + i) * 4 + c] * QuantumRange;
			outPixels[((h - 1 - j) * w + i) * 4 + 3] = buf[(j * w + i) * 4 + 3] * QuantumRange;
		}
	}

	outImage.write(fileName);
}
