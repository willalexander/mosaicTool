#define PI 3.14159																						// Handy to have
#define SUBSAMPLES 16																					// How many tile subsections (per X/Y) to use for matching original image
#define JC_VORONOI_IMPLEMENTATION																		// For Voronoi library
#define MIN(a,b) ((a < b)? a : b)																		// Simple min function
#define MAX(a,b) ((a > b)? a : b)																		// Simple max function
#define CLAMP(v, a, b) (MIN(MAX(v, a), b))																// Simple clamp function
#define OUTSIDE_IMAGE_BOUNDS(i, j, w, h) ((i<0)?true:((i>=w)?true:((j<0)?true:((j>=h)?true:false))))	// Whether or not the i,j pixel position is outside the image bounds:

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <filesystem>
#include <algorithm>
#include "types.h"
#include "imageio.h"
#include "jc_voronoi.h"
#include "mosaic.h"


using namespace std::experimental::filesystem;

void generateTilePattern(int, int);																		// Setup & record tile pattern properties 
void getTileProperties(int, pixVal *, pixVal *, float *, pixVal *, float *);							// Return geometric properties of a particular image tile
void placeTile(int, float *, pixVal, pixVal, pixVal, pixVal, float *);									// Draw a tile to the output image

void generateRegularTilePattern(int, int);																// Specific generation functions for each pattern type
void generateVoronoiPattern(int, int);
void generateVaryingVoronoiPattern(int, int);

void getTilePropertiesRegular(int, pixVal *, pixVal *, float *, pixVal *, float *);						// Specific properties functions for each pattern type
void getTilePropertiesVoronoi(int, pixVal *, pixVal *, float *, pixVal *, float *);
void getTilePropertiesVaryingVoronoi(int, pixVal *, pixVal *, float *, pixVal *, float *);

void placeTileRegular(int, float *, pixVal, pixVal, pixVal, pixVal, float *);							// Specific placement functions for each pattern type
void placeTileVoronoi(int, float *, pixVal, pixVal, pixVal, pixVal, float *);
void placeTileVaryingVoronoi(int, float *, pixVal, pixVal, pixVal, pixVal, float *);

void outputTilePixel(float *, float *, pixVal, pixVal, pixVal, pixVal);									// Draws a particular pixel the the output image

std::vector<std::string> findLibraryImages(std::string);												// Returns a vector of all images in the lib directory
void loadLibraryImages(std::string);																	// Reads into memory all iamges in the library directory
void loadLibraryImage(std::string, int);																// Reads a particular library image into memory 
void xform(int, int, float, float, int, int, int *, int *);												// Performs rotation & scaling of images
pixVal computeAverage(float *, int, int, int, int, int, int);											// Computes the averge rgb value of a rectangular section of an image buffer
pixVal computeContrast(float *, pixVal, int, int, int, int, int, int);									// Computes the contrast of a rectangular section of an image buffer

void generateVoronoiDiagram();																			// Uses the jc_voronoi library to create a voronoi diagram of tiles
bboxf voronoiCellBoundingBoxf(const jcv_site *);														// Returns the axis aligned bounding box of a voronoi cell tile
bool pixelIsInCell(const jcv_site *, int, int);															// Returns whether or nor a pixel is inside a voronoi cell tile

void generateDensityMap();																				// Creates a density map for the input image, based on variance
float densityMapLookup(int, int);																		// Returns the value of the density map at a point
float computeIntegral();																				// Integrates the density map over the image (for use in normalisation)
void generateSampler();																					// Creates an object to sample the image proportionally to the density map
void cleanUpSampler();																					// Frees memory etc.
void generateSample(int *, int *);																		// Uses PDF sampler object to generate a random sample in the image plane

void cleanup();																							// Free memory etc.


float subtileFactor;																					// The size of the subtiles relative to tiles
float *inImg, *outImg;																					// Input & output image floating point pixel buffers
int inW, inH;																							// The diemsions of the input image
int outW, outH;																							// The dimensions of the output image
float inAspectRatio;																					// Input image aspect ratio
int tileWidth, tileHeight;																				// Pixel width & height of a tile
int numTiles;																							// The number of tile in the moasic
int patternType;																						// The type of tile pattern 
float tileAspectRatio;																					// Aspect ratio of tile images after loading from the library
float cheat;																							// How much to adjust the rgb & contrast of tiles to fit the in image				
int numLibImageVariants;																				// Number of scale/rotation/colour variants of each library image to create 	
float libImageVariation;																				// How much scale/colour variation to apply to each library image
point tileCentre;																						// Centre point (in output image space) of the current tile

std::vector<float *>images;																				// Pointers to image data for each library image
std::vector<pixVal> imageAverages;																		// Average colour value of each library image
std::vector<pixVal *> imageSubAverages;																	// Subtile average colours for each library image
std::vector<pixVal> imageContrasts;																		// Contrast value for each library image
std::vector<int> imageIDs;																				// ID value for each library image (indicates which image file it came from )
float libImageScaleMin, libImageScaleMax;																// The range of scale variation to apply to library images
float libImageRotationMax;																				// The range of rotation variation to apply to library images
float libImageColMin, libImageColMax;																	// The range of colour variation to apply to library images

// For use with the regular tile pattern type:
int numLandscapeTilesWidth, numPortraitTilesWidth;														// The number of tiles in each row type 
int th, tw;																								// Tile position in row/column indices
int numRowTiles;																						// The number of tiles in the current row
int tileWidth1, tileHeight1;																			// The dimensions of a particular tile 
int numRows;																							// The number of tile rows in the image
int rowBase;																							// Pixel position of the bottom of the current tile's row
bool portrait;																							// Whether the current tile is portrait or not

// For use with the voronoi tile pattern:
int numTilesWidth, numTilesHeight;																		// The number of tiles in each dimension
jcv_point *points;																						// For storing voronoi points
std::vector<point> voronoiSites;																		// Also for voronoi points
jcv_diagram diagram;																					// Stores the voronoi diagram
int tileci, tilecj;																						// The central pixel of a tile
bboxf cellBB;																							// The bounding box of a tile
float densityMapSampleRate;																				// Resolution of the density map
int densityMapW, densityMapH;																			// Num rows/colums of the density map
int numSamples;																							// Number of samples for used for a local density measure
float *densityMap;																						// Pointer to density map data
float *rowProb;																							// For numerical integration
float *rowCDF;																							
float **pointProb;
float **pointCDF;

//
//	run() -	called from outside, contains the mosaic creation process from start to finish:
//
void run(std::string inImgFilePath, std::string libDir, int numTiles0, int outW0, int patternType0, float cheat0, int numLibImageVariants0, float libImageVariation0)
{
	// Record key settigs chosen by user:
	patternType = patternType0;
	cheat = cheat0;
	numLibImageVariants = numLibImageVariants0;
	libImageVariation = libImageVariation0;

	// Run initial setup required by the image IO system:
	imageio_init();

	// Load the main image and determine its aspect ratio:
	inImg = readImage(inImgFilePath, &inW, &inH);
	inAspectRatio = (float)(inW) / (float)(inH);

	// The size of the subtiles relative to the tile:
	subtileFactor = 1.0 / (float)(SUBSAMPLES);

	// Generate the pattern of tiles for the output image:
	generateTilePattern(numTiles0, outW0);

	// Create a map of all tiles to record which library image has been chosen for each:
	tileInfo *IDmap = (tileInfo *)(malloc(numTiles * sizeof(tileInfo)));
	for (int t = 0; t < numTiles; t++) IDmap[t].id = -1;

	// Load up the library images:
	loadLibraryImages(libDir);

	// Allocate memory for the output image:
	float *outImg = (float *)(malloc(outW * outH * 4 * sizeof(float)));


	pixVal av;																							// Average colour of a tile in the input image
	pixVal *avs = (pixVal *)malloc(SUBSAMPLES*SUBSAMPLES * sizeof(pixVal));								// Average colour of a tile subtile in the input image
	float *avs_lum = (float *)malloc(SUBSAMPLES*SUBSAMPLES * sizeof(float));							// Average brightness of a subtile in the input image									
	float *weights = (float *)(malloc(SUBSAMPLES*SUBSAMPLES * sizeof(float)));							// Subtile weight (proportional to how much of its area is visible)
	pixVal contrast;																					// Contrast of a tile in the input image

	// For each tile, find the library image which best matches the input image, and place it in the output image:
	for(int t = 0; t < numTiles; t++)
	{
		std::cout << "Creating tile " << t << " of " << numTiles << std::endl;
		
		// From the tile index t, determine the geometry of this tile & colour properties of the tile's area in the main image:
		getTileProperties(t, &av, avs, avs_lum, &contrast, weights);

		// Loop over all the library images, determine a score for how similar each is to this tile, and choose the image with the best (lowest) score:
		float lowestScore = 1000.0;
		int chosenImg;

		for (int i = 0; i < images.size(); i++)
		{
			float score = 0.0;

			// Only consider this library image if it or a variant it of it has *not* been used nearby (within 10% of the image width):
			bool fail = false;
			for (int tt = 0; tt < numTiles; tt++)
			{
				if (IDmap[tt].id == -1) continue;

				if (IDmap[tt].id == imageIDs.at(i))
				{
					float dist = sqrt(pow(tileCentre.x - IDmap[tt].centre.x, 2) + pow(tileCentre.y - IDmap[tt].centre.y, 2));
					if (dist < 0.1*(float)(outW)) fail = true;
				}
			}
			if (fail) continue;

			// Compute a score based on the closeness of the colour of each library image subtile to its equivalent in the main image tile:
			for (int s = 0; s < SUBSAMPLES*SUBSAMPLES; s++)
			{
				if(avs[s].c[0] == -1) continue;
				pixVal imgAv = imageSubAverages.at(i)[s];
				float imgAv_lum = (imgAv.c[0] + imgAv.c[1] + imgAv.c[2]) / 3.0;

				// The score is based 20% on colour differance, and 80% on luminance difference:
				for (int c = 0; c < 3; c++) score += weights[s] * 1.0 * (pow(avs[s].c[c] - imgAv.c[c], 2) / (3.0 * SUBSAMPLES*SUBSAMPLES));
				score += weights[s] * 0.0 * pow(avs_lum[s] - imgAv_lum, 2) / (SUBSAMPLES*SUBSAMPLES);
			}

			// f this score is the lowest score so far, keep track of this library image:
			if(score < lowestScore)
			{
				chosenImg = i;
				lowestScore = score;
			}
		}

		// Record which image has been chosen:
		IDmap[t].id = imageIDs.at(chosenImg);
		IDmap[t].centre = tileCentre;
		
		// Get the buffer for the chosen image:
		float *tile_buf = images.at(chosenImg);

		// Determine, depending on the cheat value, how much to adjust the image's hue and contrast to match the original image:
		pixVal hueMult;
		pixVal tileAv0 = imageAverages.at(chosenImg);
		pixVal tileContrast0 = imageContrasts.at(chosenImg);
		pixVal tileAv1, tileContrast1;
		for (int c = 0; c < 3; c++)
		{
			hueMult.c[c] = (1.0 - cheat) + cheat * av.c[c] / imageAverages.at(chosenImg).c[c];
			tileAv1.c[c] = tileAv0.c[c] * hueMult.c[c];
			tileContrast1.c[c] = tileContrast0.c[c] * hueMult.c[c];
		}

		// Draw this tile to the output mosaic image:
		placeTile(t, tile_buf, hueMult, tileAv1, contrast, tileContrast1, outImg);
	}

	// Output the final image:
	writeImage(outImg, "out.tif", outW, outH);

	// Clean up:
	free(IDmap);
	free(avs);
	free(avs_lum);
	free(weights);
	free(inImg);
	free(outImg);
	if(patternType == VARYING_VORONOI_PATTERN_TYPE) cleanUpSampler();

	cleanup();
}

//
//	generateTilePattern() - Setup & record tile pattern properties 
//
void generateTilePattern(int numTiles0, int outW0)
{
	if(patternType == REGULAR_PATTERN_TYPE) generateRegularTilePattern(numTiles0, outW0);
	if(patternType == VORONOI_PATTERN_TYPE) generateVoronoiPattern(numTiles0, outW0);
	if(patternType == VARYING_VORONOI_PATTERN_TYPE) generateVaryingVoronoiPattern(numTiles0, outW0);
}

// 
//	 getTileProperties() - Return geometric properties of a particular image tile
//
void getTileProperties(int t, pixVal *av, pixVal *avs, float *avs_lum, pixVal *contrast, float *weights)
{
	if(patternType == REGULAR_PATTERN_TYPE) getTilePropertiesRegular(t, av, avs, avs_lum, contrast, weights);
	if(patternType == VORONOI_PATTERN_TYPE) getTilePropertiesVoronoi(t, av, avs, avs_lum, contrast, weights);
	if(patternType == VARYING_VORONOI_PATTERN_TYPE) getTilePropertiesVaryingVoronoi(t, av, avs, avs_lum, contrast, weights);
}

//
//	placeTile() - Draw a tile to the output image
//
void placeTile(int t, float *tile_buf, pixVal hueMult, pixVal tileAv1, pixVal contrast, pixVal tileContrast1, float *outImg)
{
	if(patternType == REGULAR_PATTERN_TYPE) placeTileRegular(t, tile_buf, hueMult, tileAv1, contrast, tileContrast1, outImg);
	if(patternType == VORONOI_PATTERN_TYPE) placeTileVoronoi(t, tile_buf, hueMult, tileAv1, contrast, tileContrast1, outImg);
	if(patternType == VARYING_VORONOI_PATTERN_TYPE) placeTileVaryingVoronoi(t, tile_buf, hueMult, tileAv1, contrast, tileContrast1, outImg);
}


//
//	generateRegularTilePattern() - sets up the regular tile pattern
//
void generateRegularTilePattern(int numTiles0, int outW0)
{
	// Use 4/3 tiles:
	tileAspectRatio = 1.33333333;

	// Determine the number of landscape tiles in a row in the input image, based un the user's chosen total number of tiles:
	numLandscapeTilesWidth = (int)((float)(inW) / sqrt(((float)(inW) * (float)(inH) / (float)(numTiles0)) * tileAspectRatio));

	// Determine the width of the landscape tiles in the output image: 
	tileWidth = (int)((float)(outW0) / (float)(numLandscapeTilesWidth));

	// Hence determine the tile height:
	tileHeight = (int)((float)(tileWidth) / tileAspectRatio);

	// The output image width must be a multiple of the tile width:
	outW = tileWidth * numLandscapeTilesWidth;

	// The exact output image height is determined by the output image width and the input image aspect ratio:
	outH = (int)((float)(outW) / inAspectRatio);

	// Determine the number of portrait tiles across a row of the output image (it will almost certainly overshoot the right edge)
	numPortraitTilesWidth = (int)((float)(outW) / (float)(tileHeight)) + 1;

	// The number of rows is the smallest integer that covers the height of the output image:
	numRows = 2 * (int)(((float)(outH) / (float)(tileWidth + tileHeight)) + 1.0);

	// Determine exact the total number of tiles:
	numTiles = (numRows / 2) * (numLandscapeTilesWidth + numPortraitTilesWidth);
}


//
//	 generateVoronoiPattern() - sets up the voronoi tile pattern
//
void generateVoronoiPattern(int numTiles0, int outW0)
{
	// Compute the image dimensions based on user's setting:
	outW = outW0;
	outH = (int)((float)(outW) / inAspectRatio);

	// Generate uniformly spaced random voronoi sites over the output image area:
	numTilesWidth = (int)(sqrt((float)(numTiles0) * (float)(inW) / (float)(inH)));
	numTilesHeight = (int)((float)(numTilesWidth) * (float)(outH) / (float)(outW));
	numTiles = numTilesWidth * numTilesHeight;

	generateVoronoiDiagram();
	memset(&diagram, 0, sizeof(jcv_diagram));
	_jcv_rect imgBbox = { { 0.0, 0.0 },{ (float)(outW), (float)(outH) } };
	jcv_diagram_generate(numTiles, points, &imgBbox, 0, &diagram);

	// The voronoi generator may have pruned some points:
	numTiles = diagram.numsites;

	// Use 100x100 for all library images for now:
	tileWidth = tileHeight = 100;
	tileAspectRatio = 1.0;
	tileci = (int)(0.5 * (float)(tileWidth));
	tilecj = (int)(0.5 * (float)(tileHeight));
}


//
//	 generateVaryingVoronoiPattern() - sets up the voronoi tile pattern
//
void generateVaryingVoronoiPattern(int numTiles0, int outW0)
{
	// Out image dimensions based on user's setting:
	outW = outW0;
	outH = (int)((float)(outW) / inAspectRatio);

	generateDensityMap();

	// Generate uniformly spaced random voronoi sites over the output image area:
	numTilesWidth = (int)(sqrt((float)(numTiles0) * (float)(inW) / (float)(inH)));
	numTilesHeight = (int)((float)(numTilesWidth) * (float)(outH) / (float)(outW));
	numTiles = numTilesWidth * numTilesHeight;

	generateVoronoiDiagram();
	memset(&diagram, 0, sizeof(jcv_diagram));
	_jcv_rect imgBbox = { { 0.0, 0.0 },{ (float)(outW), (float)(outH) } };
	jcv_diagram_generate(numTiles, points, &imgBbox, 0, &diagram);

	// The voronoi generator may have pruned some points:
	numTiles = diagram.numsites;

	// Use 100x100 for all library images for now:
	tileWidth = tileHeight = 100;
	tileAspectRatio = 1.0;
	tileci = (int)(0.5 * (float)(tileWidth));
	tilecj = (int)(0.5 * (float)(tileHeight));
}

//
//	getTilePropertiesRegular() - given a tile index, computes geometric properties of that tile & colour properties of the tile's contents in the input image
//
void getTilePropertiesRegular(int t, pixVal *av, pixVal *avs, float *avs_lum, pixVal *contrast, float *weights)
{
	// First determine which landscape/portrait row pair this tile is on:
	int rowPair = (int)((float)(t) / (float)(numLandscapeTilesWidth + numPortraitTilesWidth));

	// Hence count tiles from start of landscape row to determine whether it's in the landscape row (above) or the portrait row (below)
	tw = t - rowPair * (numLandscapeTilesWidth + numPortraitTilesWidth);

	// Depending on whether tile is on a landscape or portrait row, compute info about tile position & dimensions:
	if(tw < numLandscapeTilesWidth)
	{
		th = rowPair * 2;

		tileWidth1 = tileWidth;
		tileHeight1 = tileHeight;
		rowBase = th * (tileWidth + tileHeight) / 2;
		numRowTiles = numLandscapeTilesWidth;
		portrait = false;
	}

	// if portrait:
	else
	{
		th = rowPair * 2 + 1;
		tw -= numLandscapeTilesWidth;

		tileWidth1 = tileHeight;
		tileHeight1 = tileWidth;
		rowBase = (int)(th / 2.0) * (tileWidth + tileHeight) + tileHeight;
		numRowTiles = numPortraitTilesWidth;
		portrait = true;
	}

	// Compute the centre point (in output image space) of this tile:
	tileCentre = { (float)(((float)(tw)+0.5) * (float)(tileWidth1)), (float)((float)(rowBase)+0.5 * (float)(tileHeight1)) };

	// Determine the range of pixels in the original image which occupy this tile:
	int i0 = (int)((float)(inW) * (float)(tw * tileWidth1) / (float)(outW));
	int i1 = (int)((float)(inW) * (float)((tw + 1) * tileWidth1) / (float)(outW)) - 1;
	int j0 = (int)((float)(inH) * (float)(rowBase) / (float)(outH));
	int j1 = (int)((float)(inH) * (float)(rowBase + tileHeight1) / (float)(outH)) - 1;

	// Compute average colour of the tile:
	*av = computeAverage(inImg, inW, inH, i0, i1 + 1, j0, j1 + 1);

	// Compute the average colour value of each subsample:
	for (int si = 0; si < SUBSAMPLES; si++)
	{
		for (int sj = 0; sj < SUBSAMPLES; sj++)
		{
			// Determine the set of pixels in this subtile:
			int si0 = (int)((float)(inW) * (float)(tw * tileWidth1 + (float)(si)* subtileFactor * (float)(tileWidth)) / (float)(outW));
			int si1 = (int)((float)(inW) * (float)(tw * tileWidth1 + (float)(si + 1)* subtileFactor * (float)(tileWidth)) / (float)(outW)) - 1;
			int sj0 = (int)((float)(inH) * (float)(rowBase + (float)(sj)* subtileFactor * (float)(tileHeight)) / (float)(outH));
			int sj1 = (int)((float)(inH) * (float)(rowBase + (float)(sj + 1)* subtileFactor * (float)(tileHeight)) / (float)(outH)) - 1;

			if (portrait)
			{
				si0 = (int)((float)(inW) * (float)(tw * tileWidth1 + (float)(sj)* subtileFactor * (float)(tileHeight)) / (float)(outW));
				si1 = (int)((float)(inW) * (float)(tw * tileWidth1 + (float)(sj + 1)* subtileFactor * (float)(tileHeight)) / (float)(outW)) - 1;
				sj0 = (int)((float)(inH) * (float)(rowBase + (float)(si)* subtileFactor * (float)(tileWidth)) / (float)(outH));
				sj1 = (int)((float)(inH) * (float)(rowBase + (float)(si + 1)* subtileFactor * (float)(tileWidth)) / (float)(outH)) - 1;
			}

			avs[sj * SUBSAMPLES + si] = computeAverage(inImg, inW, inH, si0, si1 + 1, sj0, sj1 + 1);

			// Compute the average luminance value of this sub tile:
			avs_lum[sj * SUBSAMPLES + si] = (avs[sj * SUBSAMPLES + si].c[0] + avs[sj * SUBSAMPLES + si].c[1] + avs[sj * SUBSAMPLES + si].c[2]) / 3.0;

			// For a regular tile pattern, each subtile is the same size, so carries the same weight:
			weights[sj * SUBSAMPLES + si] = 1.0;
		}
	}

	// Compute the contrast:
	*contrast = computeContrast(inImg, *av, inW, inH, i0, i1 + 1, j0, j1 + 1);
}

//
//	getTilePropertiesVoronoi() - given a tile index, computes geometric properties of that tile & colour properties of the tile's contents in the input image
//
void getTilePropertiesVoronoi(int t, pixVal *av, pixVal *avs, float *avs_lum, pixVal *contrast, float *weights)
{
	// Get reference to the voronoi site for this tile:
	const jcv_site *sites = jcv_diagram_get_sites(&diagram);

	const jcv_site *site = &sites[t];

	// Generate the cell's bounding box:
	cellBB = voronoiCellBoundingBoxf(site);

	tileCentre = {(float)(site->p.x), (float)(site->p.y)};

	int count;							
	
	// Convert the cell bounding box to main image space:
	float i0 = cellBB.i0 * (float)(inW) / (float)(outW);
	float i1 = cellBB.i1 * (float)(inW) / (float)(outW);
	float j0 = cellBB.j0 * (float)(inH) / (float)(outH);
	float j1 = cellBB.j1 * (float)(inH) / (float)(outH);

	// Compute average:
	for (int c = 0; c < 3; c++) av->c[c] = 0.0;
	count = 0;
	for (float fj = (int)(j0)+0.5; fj < j1; fj++)
	{
		for (float fi = (int)(i0)+0.5; fi < i1; fi++)
		{
			if ((fj < j0) || (fi < i0)) continue;

			// To check whether this main image pixel is within the voronoi cell, convert it back to out image space:
			int itest = (int)(fi * (float)(outW) / (float)(inW));
			int jtest = (int)(fj * (float)(outH) / (float)(inH));

			if (pixelIsInCell(site, itest, jtest))
			{
				int i = (int)(fi);
				int j = (int)(fj);

				if ((i < 0) || (i >= inW) || (j < 0) || (j >= inH)) continue;

				for (int c = 0; c < 3; c++) av->c[c] += inImg[(j * inW + i) * 3 + c];
				count++;
			}
		}
	}
	for (int c = 0; c < 3; c++) av->c[c] /= (float)(count);

	// Compute the average value of each subsample in the main image:
	float maxPixels = pow((float)(inW) * ((cellBB.i1 - cellBB.i0) / (float)(SUBSAMPLES)) / (float)(outW), 2);

	// The cell's bounding box is not square, so make sure to compensate when converting from bbox space to image space:
	float wmult = 1.0;
	float hmult = 1.0;
	if ((cellBB.i1 - cellBB.i0) >(cellBB.j1 - cellBB.j0))	hmult = (cellBB.i1 - cellBB.i0) / (cellBB.j1 - cellBB.j0);
	else wmult = (cellBB.j1 - cellBB.j0) / (cellBB.i1 - cellBB.i0);


	float ci = (int)(cellBB.i0 + 0.5 * (float)(cellBB.i1 - cellBB.i0));
	float cj = (int)(cellBB.j0 + 0.5 * (float)(cellBB.j1 - cellBB.j0));

	for (int si = 0; si < SUBSAMPLES; si++)
	{
		for (int sj = 0; sj < SUBSAMPLES; sj++)
		{
			for (int c = 0; c < 3; c++) avs[sj * SUBSAMPLES + si].c[c] = 0.0;

			// Compute the subcell bounding box in out image space;
			float i0_0 = (cellBB.i0 + si * ((cellBB.i1 - cellBB.i0) / (float)(SUBSAMPLES)));
			float i1_0 = (cellBB.i0 + (si + 1) * ((cellBB.i1 - cellBB.i0) / (float)(SUBSAMPLES)));
			float j0_0 = (cellBB.j0 + sj * ((cellBB.j1 - cellBB.j0) / (float)(SUBSAMPLES)));
			float j1_0 = (cellBB.j0 + (sj + 1) * ((cellBB.j1 - cellBB.j0) / (float)(SUBSAMPLES)));

			// Compensate for non-square cell shape:
			i0_0 = (i0_0 - ci) * wmult + ci;
			i1_0 = (i1_0 - ci) * wmult + ci;
			j0_0 = (j0_0 - cj) * hmult + cj;
			j1_0 = (j1_0 - cj) * hmult + cj;

			// Convert to main image space:
			int i0 = (int)((float)(inW)* i0_0 / (float)(outW));
			int i1 = (int)((float)(inW)* i1_0 / (float)(outW)) - 1;
			int j0 = (int)((float)(inH)* j0_0 / (float)(outH));
			int j1 = (int)((float)(inH)* j1_0 / (float)(outH)) - 1;

			count = 0;
			for (int i = i0; i <= i1; i++)
			{
				for (int j = j0; j <= j1; j++)
				{
					if ((i < 0) || (i >= inW) || (j < 0) || (j >= inH)) continue;

					// Only consider pixels that are inside the cell:
					int testi = (int)((float)(outW) * (i + 0.5) / (float)(inW));
					int testj = (int)((float)(outH) * (j + 0.5) / (float)(inH));

					if (pixelIsInCell(site, testi, testj) == false) continue;

					for (int c = 0; c < 3; c++) avs[sj * SUBSAMPLES + si].c[c] += inImg[(j * inW + i) * 3 + c];
					count++;
				}
			}
			for (int c = 0; c < 3; c++) avs[sj * SUBSAMPLES + si].c[c] /= (float)(count);
			if (count == 0)
			{
				for (int c = 0; c < 3; c++) avs[sj * SUBSAMPLES + si].c[c] = -1;
			}

			weights[sj * SUBSAMPLES + si] = (float)(count) / maxPixels;
			if (weights[sj * SUBSAMPLES + si] > 1.0) weights[sj * SUBSAMPLES + si] = 1.0;

			// Compute the average luminance value of this sub tile:
			avs_lum[sj * SUBSAMPLES + si] = (avs[sj * SUBSAMPLES + si].c[0] + avs[sj * SUBSAMPLES + si].c[1] + avs[sj * SUBSAMPLES + si].c[2]) / 3.0;
		}
	}
	
	// Compute the contrast:
	for (int c = 0; c < 3; c++) contrast->c[c] = 0.0;
	count = 0;
	for (float fj = (int)(j0)+0.5; fj < j1; fj++)
	{
		for (float fi = (int)(i0)+0.5; fi < i1; fi++)
		{
			if ((fj < j0) || (fi < i0)) continue;

			// To check whether this main image pixel is within the voronoi cell, convert it back to out image space:
			int itest = (int)(fi * (float)(outW) / (float)(inW));
			int jtest = (int)(fj * (float)(outH) / (float)(inH));

			if (pixelIsInCell(site, itest, jtest))
			{
				int i = (int)(fi);
				int j = (int)(fj);

				if ((i < 0) || (i >= inW) || (j < 0) || (j >= inH)) continue;

				for (int c = 0; c < 3; c++) contrast->c[c] += fabs(inImg[(j * inW + i) * 3 + c] - av->c[c]);
				count++;
			}
		}
	}
	for (int c = 0; c < 3; c++) contrast->c[c] /= (float)(count);
}

//
//	getTilePropertiesVaryingVoronoi() - given a tile index, computes geometric properties of that tile & colour properties of the tile's contents in the input image
//
void getTilePropertiesVaryingVoronoi(int t, pixVal *av, pixVal *avs, float *avs_lum, pixVal *contrast, float *weight)
{
	getTilePropertiesVoronoi(t, av, avs, avs_lum, contrast, weight);
}

//
//	placeTileRegular() - draws a given library image into the tile's position
//
void placeTileRegular(int t, float *tile_buf, pixVal hueMult, pixVal tileAv1, pixVal contrast, pixVal tileContrast1, float *outImg)
{
	// For each pixel in the output image, determine the source pixel from the library image:
	for (int j = rowBase; j < (rowBase + tileHeight1); j++)
	{
		for (int i = tw * tileWidth1; i < (tw + 1) * tileWidth1; i++)
		{
			if ((i < 0) || (i >= outW) || (j < 0) || (j >= outH)) continue;
			int i0 = i - tw * tileWidth1;
			int j0 = j - rowBase;
			int ti = i0;
			int tj = j0;
			if (portrait)
			{
				ti = j0;
				tj = i0;
			}
			if ((ti < 0) || (ti >= tileWidth) || (tj < 0) || (tj >= tileHeight)) continue;

			// Adjust the lib image's hue and contrast to more closely match the main image:
			outputTilePixel(tile_buf + (tj * tileWidth + ti) * 3, outImg + ((outH - j - 1) * outW + i) * 4, hueMult, tileAv1, contrast, tileContrast1);
		}
	}
}

//
//	 placeTileVoronoi() - draws a given library image into the tile's position
//
void placeTileVoronoi(int t, float *tile_buf, pixVal hueMult, pixVal tileAv1, pixVal contrast, pixVal tileContrast1, float *outImg)
{
	// Get reference to the voronoi site for this tile:
	const jcv_site *sites = jcv_diagram_get_sites(&diagram);
	const jcv_site *site = &sites[t];

	int ci = (int)(cellBB.i0 + 0.5 * (float)(cellBB.i1 - cellBB.i0));
	int cj = (int)(cellBB.j0 + 0.5 * (float)(cellBB.j1 - cellBB.j0));

	float vcellw = (float)(cellBB.i1 - cellBB.i0);
	float vcellh = (float)(cellBB.j1 - cellBB.j0);

	float wmult = 1.0;
	float hmult = 1.0;

	if(vcellw > vcellh)	hmult = vcellh / vcellw;
	if(vcellh > vcellw)	wmult = vcellw / vcellh;


	for (int i = (int)(cellBB.i0); i <= (int)(cellBB.i1 + 1.0); i++)
	{
		for (int j = (int)(cellBB.j0); j <= (int)(cellBB.j1 + 1.0); j++)
		{
			if (pixelIsInCell(site, i, j))
			{
				if ((i < 0) || (i >= outW) || (j < 0) || (j >= outH)) continue;

				// 'relative' i & j, i.e. relative to the bottom left corner of the voronoi cell's bounding box: (scale to compensate for fact that the cell's bounding box is not square)
				float irel = (float)(i - ci + 0.5) * wmult + (ci - (int)(cellBB.i0));
				float jrel = (float)(j - cj + 0.5) * hmult + (cj - (int)(cellBB.j0));

				int ti = CLAMP((int)((float)(tileWidth) * irel / (float)(cellBB.i1 - cellBB.i0)), 0, tileWidth);
				int tj = CLAMP((int)((float)(tileHeight)* jrel / (float)(cellBB.j1 - cellBB.j0)), 0, tileHeight);

				// Adjust the lib image's hue and contrast to more closely match the main image:
				outputTilePixel(tile_buf + (tj * tileWidth + ti) * 3, outImg + ((outH - j - 1) * outW + i) * 4, hueMult, tileAv1, contrast, tileContrast1);
			}
		}
	}
}

//
//	 placeTileVaryingVoronoi() - draws a given library image into the tile's position
//
void placeTileVaryingVoronoi(int t, float *tile_buf, pixVal hueMult, pixVal tileAv1, pixVal contrast, pixVal tileContrast1, float *outImg)
{
	placeTileVoronoi(t, tile_buf, hueMult, tileAv1, contrast, tileContrast1, outImg);
}

//
//	findLibrayImages() - Returns a vector of all images in the lib directory
//
std::vector<std::string> findLibraryImages(std::string libDir)
{
	std::vector<std::string> fileNames;

	// Find any jpg images in this current directory:
	for (directory_iterator next(path(libDir + "\\")), end; next != end; ++next)
	{
		path  myPath = next->path();
		std::string fileName(myPath.generic_string());
		fileNames.push_back(fileName);
	}

	return fileNames;
}

//
//	loadLibraryImage() - Reads into memory all iamges in the library directory
//
void loadLibraryImages(std::string libDir)
{
	// Get a list of all possible library image files in this directory:
	std::vector<std::string> fileNames = findLibraryImages(libDir);

	// If the user has specified more than 0 library image variants, set the range of scale and colour adjustment here based on the variation setting:
	if (numLibImageVariants > 0)
	{
		libImageScaleMin = 1.0 - libImageVariation;
		libImageScaleMax = 1.0;

		libImageRotationMax = 2.0 * PI;

		libImageColMin = 1.0 - 0.5 * libImageVariation;
		libImageColMax = 1.0 + 0.5 * libImageVariation;
	}

	// If no variants are requested, then use each library image without any variation of scale, rotation or colour:
	else
	{
		libImageScaleMin = libImageScaleMax = 1.0;
		libImageRotationMax = 0.0;
		libImageColMin = libImageColMax = 1.0;
	}

	// Load each of these images, downsize them to the tile size, store them, and record the average rgb value of their pixels:
	for (int f = 0; f < fileNames.size(); f++)
	{
		std::cout << "Loading " << fileNames.at(f) << "(" << (f + 1) << " of " << fileNames.size() << ")" << std::endl;
		loadLibraryImage(fileNames.at(f), f);
	}
}

//
//	loadLibraryImage() - Reads a particular library image into memory 
//
void loadLibraryImage(std::string fileName, int ID)
{
	int w, h;
	float scaleFactor;
	int i0, i1, j0, j1;

	// Load the image pixel data:
	float *img_buf = readImage(fileName, &w, &h);

	// Scale the image to fit the tile shape. If tile shape is wider, scale by width. If taller, scale by height:
	if (((float)(w)/(float)(h)) > tileAspectRatio)
	{
		scaleFactor = (float)(tileHeight) / (float)(h);
		i0 = (int)(0.5 * (w - (float)(tileWidth) / scaleFactor));
		i1 = w - 1 - i0;
		j0 = 0;
		j1 = h - 1;
	}
	else
	{
		scaleFactor = (float)(tileWidth) / (float)(w);
		i0 = 0;
		i1 = w - 1;
		j0 = (int)(0.5 * (h - (float)(tileHeight) / scaleFactor));
		j1 = h - 1 - j0;
	}

	// Create several tile-sized versions of the image, scaled and rotated & colour adjusted for more variation:
	for(int v = 0; v < MAX(numLibImageVariants, 1); v++)
	{
		// For each variation, choose a random scale, rotation and tint:
		float scale = libImageScaleMin + ((float)(rand())/(float)(RAND_MAX)) * (libImageScaleMax - libImageScaleMin);
		float theta = libImageRotationMax * (float)(rand()) / (float)(RAND_MAX);
		pixVal colRand;
		for(int c = 0; c < 3; c++) colRand.c[c] = libImageColMin + ((float)(rand()) / (float)(RAND_MAX)) * (libImageColMax - libImageColMin);

		// Allocate memory for this version:
		float *tileBuf = (float *)(malloc(tileWidth * tileHeight * 3 * sizeof(float)));

		// Fill the tile pixels from the image buffer:
		for(int ti = 0; ti < tileWidth; ti++)
		{
			for(int tj = 0; tj < tileHeight; tj++)
			{
				// Determine which pixel in the main image corresponds to this pixel in the tile image:
				int i, j;
				xform(ti, tj, scale, theta, i1 + 1 - i0, j1 + 1 - j0, &i, &j);

				// Copy pixel values:
				for (int c = 0; c < 3; c++)
				{
					// Copy the pixels:
					tileBuf[(tj * tileWidth + ti) * 3 + c] = CLAMP(img_buf[((j + j0) * w + (i + i0)) * 3 + c] * colRand.c[c], 0.0f, 1.0f);
				}
			}
		}

		// Compute the average colour of this new tile image:
		pixVal av = computeAverage(tileBuf, tileWidth, tileHeight, 0, tileWidth, 0, tileHeight);

		// Compute the average colour of each subtile of the tile image:
		pixVal *avs = (pixVal *)(malloc(SUBSAMPLES * SUBSAMPLES * sizeof(pixVal)));

		for (int si = 0; si < SUBSAMPLES; si++)
		{
			for (int sj = 0; sj < SUBSAMPLES; sj++)
			{
				int i0 = (int)(si * (float)(tileWidth)* subtileFactor);
				int i1 = (int)((si + 1) * (float)(tileWidth)* subtileFactor) - 1;
				int j0 = (int)(sj * (float)(tileHeight)* subtileFactor);
				int j1 = (int)((sj + 1) * (float)(tileHeight)* subtileFactor) - 1;

				avs[sj * SUBSAMPLES + si] = computeAverage(tileBuf, tileWidth, tileHeight, i0, i1+1, j0, j1+1);
			}
		}

		// Compute the image's contrast:
		pixVal contrast = computeContrast(tileBuf, av, tileWidth, tileHeight, 0, tileWidth, 0, tileHeight);

		// Add this tile image & info to the global list:
		images.push_back(tileBuf);
		imageAverages.push_back(av);
		imageSubAverages.push_back(avs);
		imageContrasts.push_back(contrast);
		imageIDs.push_back(ID);
	}

	free(img_buf);
}

//
//	xform() - Performs rotation & scaling of images
//
void xform(int ti, int tj, float scale, float theta, int w, int h, int *i, int *j)
{
	// Convert the tile image coordinates into scalar coordinates with origin in the centre:
	float x = ti + 0.5 - 0.5*(float)(tileWidth);
	float y = tj + 0.5 - 0.5*(float)(tileHeight);

	// Apply scale factor:
	x /= scale;
	y /= scale;

	// Now convert this xformed coordinate to the size of the input image:
	x *= ((float)(w) / (float)(tileWidth));
	y *= ((float)(h) / (float)(tileHeight));

	// Rotate:
	float r = sqrt(x*x + y * y);
	float angle0 = atan(x / y);
	if (y < 0.0) angle0 += PI;

	float angle1 = angle0 + theta;
	float x1 = r * sin(angle1);
	float y1 = r * cos(angle1);

	// Finally convert these coordinates back to integer pixel coordinates:
	*i = (int)(x1 + 0.5*(float)(w));
	*j = (int)(y1 + 0.5*(float)(h));

	// extrapolate pixel values that are outside the bounds of the original image:
	if (*i < 0)	*i = 0;
	if (*i >= w) *i = w - 1;
	if (*j < 0) *j = 0;
	if (*j >= h) *j = h - 1;
}

//
//	computeAverage() - Computes the averge rgb value of a rectangular section of an image buffer
//
pixVal computeAverage(float *buf, int w, int h, int i0, int i1, int j0, int j1)
{
	pixVal av = { 0.0, 0.0, 0.0 };

	for(int c = 0; c < 3; c++)
	{
		int count = 0;

		for(int i = i0; i < i1; i++)
		{
			for(int j = j0; j < j1; j++)
			{
				// Skip if this pixel is outside the image bounds:
				if(OUTSIDE_IMAGE_BOUNDS(i, j, w, h)) continue;

				av.c[c] += buf[(j * w + i) * 3 + c];
				count++;
			}
		}

		av.c[c] /= (float)(count);
	}

	return av;
}

//
//	computeContrast() - Computes the contrast of a rectangular section of an image buffer
//
pixVal computeContrast(float *buf, pixVal av, int w, int h, int i0, int i1, int j0, int j1)
{
	pixVal contrast = {0.0, 0.0, 0.0};

	for(int c = 0; c < 3; c++)
	{
		int count = 0;

		for(int i = i0; i < i1; i++)
		{
			for(int j = j0; j < j1; j++)
			{
				// Skip if this pixel is outside the image bounds:
				if(OUTSIDE_IMAGE_BOUNDS(i, j, w, h)) continue;

				contrast.c[c] += fabs(av.c[c] - buf[(j * w + i) * 3 + c]);
				count++;
			}
		}

		contrast.c[c] /= (float)(count);
	}

	return contrast;
}

//
//	generateVoronoiDiagram() -  Uses the jc_voronoi library to create a voronoi diagram of tiles
//
void generateVoronoiDiagram()
{
	points = (_jcv_point *)(malloc(numTilesWidth * numTilesHeight * sizeof(_jcv_point)));

	if(patternType == VORONOI_PATTERN_TYPE)
	{
		for (int i = 0; i < numTilesWidth; i++)
		{
			for (int j = 0; j < numTilesHeight; j++)
			{
				points[j*numTilesWidth + i] = { (jcv_real)(((float)(i)+(float)(rand()) / (float)(RAND_MAX)) * ((float)(outW) / (float)(numTilesWidth))), (jcv_real)(((float)(j)+(float)(rand()) / (float)(RAND_MAX)) * ((float)(outH) / (float)(numTilesHeight))) };
			}
		}
	}

	if(patternType == VARYING_VORONOI_PATTERN_TYPE)
	{
		generateSampler();
		for (int i = 0; i < (numTilesWidth * numTilesHeight); i++)
		{
			int pi, pj;

			generateSample(&pi, &pj);
			points[i] = { (jcv_real)(pi), (jcv_real)(pj) };
		}
	}
}

//
//	voronoiCellBoundingBoxf() - Returns the axis aligned bounding box of a voronoi cell tile
//
bboxf voronoiCellBoundingBoxf(const jcv_site *site)
{
	const jcv_graphedge *e = site->edges;

	float i0, i1, j0, j1;
	i0 = i1 = e->pos[0].x;
	j0 = j1 = e->pos[0].y;

	while (e)
	{
		for (int v = 0; v < 2; v++)
		{
			if (e->pos[v].x < i0) i0 = e->pos[v].x;
			if (e->pos[v].x > i1) i1 = e->pos[v].x;
			if (e->pos[v].y < j0) j0 = e->pos[v].y;
			if (e->pos[v].y > j1) j1 = e->pos[v].y;
		}

		e = e->next;
	}

	return { i0, i1, j0, j1 };
}

//
//	pixelIsInCell() - Returns whether or nor a pixel is inside a voronoi cell tile
//
bool pixelIsInCell(const jcv_site *site, int i, int j)
{
	// exact pos is centre of pixel:
	float x = (float)(i)+0.5;
	float y = (float)(j)+0.5;

	const jcv_graphedge *e = site->edges;

	// Loop over each cell edge, and check whether the pixel is on the 'correct' side of the edge
	while (e)
	{
		// Compute the unit normal vector of the edge:
		float normx = e->pos[0].y - e->pos[1].y;
		float normy = e->pos[1].x - e->pos[0].x;
		float len = sqrt(pow(normx, 2) + pow(normy, 2));
		normx /= len;
		normy /= len;

		// Compute the signed normal distance from the edge, for the cell centre and the pixel. (allow an extra tolerance of 1px for the pixel distance) 
		float centreDot = (site->p.x - e->pos[0].x) * normx + (site->p.y - e->pos[0].y) * normy;
		float pixelDot = (x - e->pos[0].x + normx) * normx + (y - e->pos[0].y + normy) * normy;

		// If the pixel is on the wrong side of the edge, then it is outside the cell:
		if ((centreDot < 0.0) && (pixelDot > 0.0)) return false;
		if ((centreDot > 0.0) && (pixelDot < 0.0)) return false;

		e = e->next;
	}

	return true;
}

//
//	generateDensityMap() - Creates a density map for the input image, based on variance
//
void generateDensityMap()
{
	numSamples = 64;
	densityMapSampleRate = 1.0 / (float)(numSamples);
	densityMapW = (float)(inW)* densityMapSampleRate;
	densityMapH = (float)(inH)* densityMapSampleRate;
	densityMap = (float *)(malloc(numSamples * numSamples * sizeof(float)));

	int radius = (int)(0.1 * (float)(inW));

	for (int mi = 0; mi < numSamples; mi++)
	{
		for (int mj = 0; mj < numSamples; mj++)
		{
			int ci = (int)((mi + 0.5) * ((float)(inW) / (float)(numSamples)));
			int cj = (int)((mj + 0.5) * ((float)(inH) / (float)(numSamples)));

			// Determine the average colour in this section of the original image:
			pixVal av = { 0.0, 0.0, 0.0 };
			int num = 0;
			for (int i = ci - radius; i < (ci + radius); i++)
			{
				for (int j = cj - radius; j < (cj + radius); j++)
				{
					if ((i < 0) || (i >= inW) || (j < 0) || (j >= inH)) continue;

					if (sqrt(pow(i - ci, 2) + pow(j - cj, 2)) > radius) continue;

					for (int c = 0; c < 3; c++) av.c[c] += inImg[(j * inW + i) * 3 + c];
					num++;
				}
			}
			for (int c = 0; c < 3; c++) av.c[c] /= (float)(num);

			// Determine the variance in this section:
			float variance = 0.0;
			num = 0;

			for (int i = ci - radius; i < (ci + radius); i++)
			{
				for (int j = cj - radius; j < (cj + radius); j++)
				{
					if ((i < 0) || (i >= inW) || (j < 0) || (j >= inH)) continue;
					if (sqrt(pow(i - ci, 2) + pow(j - cj, 2)) > radius) continue;

					for (int c = 0; c < 3; c++) variance += fabs(inImg[(j * inW + i) * 3 + c] - av.c[c]);
					num++;
				}
			}
			for (int c = 0; c < 3; c++) av.c[c] /= (float)(num);

			densityMap[mj * numSamples + mi] = variance;
		}
	}

	// Normalize:
	float min = 1.0, max = 0.0;
	for (int i = 0; i < numSamples*numSamples; i++)
	{
		if (densityMap[i] < min) min = densityMap[i];
		if (densityMap[i] > max) max = densityMap[i];
	}

	for (int i = 0; i < numSamples*numSamples; i++)
	{
		densityMap[i] = (densityMap[i] - min) / (max - min);
		if (densityMap[i] < 0.5) densityMap[i] = 0.5*densityMap[i];
		else densityMap[i] = 0.5 * 1.0 + 0.5 * densityMap[i];
	}
}

//
//	densityMapLookup() - Returns the value of the density map at a point
//
float densityMapLookup(int i, int j)
{
	// Determine the 4 density map data points that this pixel is closest to:
	int mi0 = (int)(((float)(i) / ((float)(outW)* densityMapSampleRate)) - 0.5);
	int mi1 = (int)(((float)(i) / ((float)(outW)* densityMapSampleRate)) + 0.5);
	int mj0 = (int)(((float)(j) / ((float)(outH)* densityMapSampleRate)) - 0.5);
	int mj1 = (int)(((float)(j) / ((float)(outH)* densityMapSampleRate)) + 0.5);

	if (mi1 >= numSamples) mi1 = numSamples - 1;
	if (mj1 >= numSamples) mj1 = numSamples - 1;

	// Determine the parametric distances to the data points:
	float ti = (i - ((mi0 + 0.5)*((float)(outW)* densityMapSampleRate))) / ((float)(outW)* densityMapSampleRate);
	float tj = (j - ((mj0 + 0.5)*((float)(outH)* densityMapSampleRate))) / ((float)(outH)* densityMapSampleRate);

	// Determine the interpolated map value at this point:
	float valj0 = (1.0 - ti)*(densityMap[mj0 * numSamples + mi0]) + (ti)*(densityMap[mj0 * numSamples + mi1]);
	float valj1 = (1.0 - ti)*(densityMap[mj1 * numSamples + mi0]) + (ti)*(densityMap[mj1 * numSamples + mi1]);

	return (1.0 - tj)*(valj0)+(tj)*(valj1);
}

//
//	computeIntegral() - Integrates the density map over the image (for use in normalisation)
//																						
float computeIntegral()
{
	float result = 0.0;

	for (int x = 0; x < numSamples*numSamples; x++) result += densityMap[x] * ((float)(outW*outH) / (float)(numSamples*numSamples));

	return result;
}

//
//	generateSampler() - Creates an object to sample the image proportionally to the density map
//																				
void generateSampler()
{
	float normFactor = 1.0 / computeIntegral();

	float integral = 0.0;
	for (int x = 0; x < numSamples*numSamples; x++) integral += normFactor * densityMap[x] * ((float)(outW*outH) / (float)(numSamples*numSamples));
	
	rowProb = (float *)(malloc(numSamples * sizeof(float)));
	pointProb = (float **)(malloc(numSamples * sizeof(float *)));
	pointCDF = (float **)(malloc(numSamples * sizeof(float *)));
	for (int r = 0; r < numSamples; r++)
	{
		rowProb[r] = 0.0;
		for (int x = 0; x < numSamples; x++) rowProb[r] += densityMap[r * numSamples + x] * normFactor * ((float)(outW*outH) / (float)(numSamples*numSamples));

		pointProb[r] = (float *)(malloc(numSamples * sizeof(float)));
		for (int x = 0; x < numSamples; x++) pointProb[r][x] = densityMap[r * numSamples + x] * normFactor * ((float)(outW*outH) / (float)(numSamples*numSamples)) * (1.0 / rowProb[r]);

		pointCDF[r] = (float *)(malloc((numSamples + 1) * sizeof(float)));
		pointCDF[r][0] = 0.0;
		for (int x = 0; x < numSamples; x++) pointCDF[r][x + 1] = pointCDF[r][x] + pointProb[r][x];

	}
	rowCDF = (float *)(malloc((numSamples + 1) * sizeof(float)));
	rowCDF[0] = 0.0;
	for (int i = 0; i < numSamples; i++)
	{
		rowCDF[i + 1] = rowCDF[i] + rowProb[i];
	}
}

//
//	cleanUpSampler() - Frees memory etc.
//
void cleanUpSampler()
{
	for (int i = 0; i < numSamples; i++) free(pointProb[i]);
	for (int i = 0; i <= numSamples; i++) free(pointCDF[i]);

	free(rowProb);
	free(rowCDF);
	free(pointProb);
	free(pointCDF);
	free(densityMap);
}

//
//	generateSample(int *, int *) - Uses PDF sampler object to generate a random sample in the image plane
//
void generateSample(int *i, int *j)
{
	float x, y;

	// Choose a row:
	float sample = ((float)(rand())) / (float)(RAND_MAX);
	int row = 0;
	for (int i = 1; i <= numSamples; i++)
	{
		if (rowCDF[i] > sample)
		{
			row = i - 1;
			break;
		}
	}

	// Choose a column:
	int column = 0;
	sample = ((float)(rand())) / (float)(RAND_MAX);
	for (int j = 1; j <= numSamples; j++)
	{
		if (pointCDF[row][j] > sample)
		{
			column = j - 1;
			break;
		}
	}

	// Generata a random sample point in this tile:
	*j = (int)((row + (float)(rand()) / (float)(RAND_MAX)) * ((float)(outH) / (float)(numSamples)));
	*i = (int)((column + (float)(rand()) / (float)(RAND_MAX)) * ((float)(outW) / (float)(numSamples)));

}

//
//	outputTilePixel() - Draws a particular pixel the the output image
//
void outputTilePixel(float *tilePtr, float *outPtr, pixVal hueMult, pixVal tileAv1, pixVal contrast, pixVal tileContrast1)
{
	for (int c = 0; c < 3; c++)
	{
		// Get the pixel compoenent value from the library image, and adjust its average hue to match this main image tile:
		float cval = tilePtr[c] * hueMult.c[c];

		// The contrast multiplier scales the pixel values relative to the mean, to change the contrast of the image:
		float contrastMultiplier = contrast.c[c] / tileContrast1.c[c];

		// The tile contrast could be zero, in which case there is no point having a contrast multiplier:
		if (std::isinf(contrastMultiplier)) contrastMultiplier = 1.0;

		// Adjust the contrast to be nearer the contrast of the main image tile:
		cval = tileAv1.c[c] + (cval - tileAv1.c[c]) * ((1.0 - 0.5*cheat) + 0.5*cheat * contrastMultiplier);

		if (cval < 0.0) cval = 0.0;
		if (cval >= 1.0) cval = 1.0;

		outPtr[c] = cval;
	}
	outPtr[3] = 1.0;
}

//
//	cleanup() -  Free memory etc.
//
void cleanup()
{
	for (int i = 0; i < images.size(); i++)
	{
		free(images.at(i));
		free(imageSubAverages.at(i));
	}
}