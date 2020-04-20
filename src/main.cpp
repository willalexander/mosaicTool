#include "mosaic.h"
#include <iostream>
#include <map>
#include <iostream>
#include <string>

int main(int argc, char **argv)
{
	// Declare variable for input arguments, and set defaults:
	std::string inImgFilePath = "None";
	std::string libDir = "None";
	int numtiles = 2000;
	int mosaicwidth = 1024;
	int pattern = REGULAR_PATTERN_TYPE;
	float cheat = 0.5;
	int numLibImageVariants = 32;
	float libImageVariation = 0.5;

	// Process command line arguments:
	for(int a = 0; a < argc; a++)
	{
		if(argv[a] == std::string("-image")) inImgFilePath = std::string(argv[a + 1]);
		if(argv[a] == std::string("-library")) libDir = std::string(argv[a + 1]);
		if(argv[a] == std::string("-numtiles")) numtiles = atoi(argv[a + 1]);
		if(argv[a] == std::string("-mosaicwidth")) mosaicwidth = atoi(argv[a + 1]);
		if(argv[a] == std::string("-pattern"))
		{
			if (argv[a + 1] == std::string( "regular")) pattern = REGULAR_PATTERN_TYPE;
			if(argv[a + 1] == std::string("voronoi")) pattern = VORONOI_PATTERN_TYPE;
			if(argv[a + 1] == std::string("varyingvoronoi")) pattern = VARYING_VORONOI_PATTERN_TYPE;
		}
		if (argv[a] == std::string("-cheat")) cheat = atof(argv[a + 1]);
		if (argv[a] == std::string("-libVariants")) numLibImageVariants = atof(argv[a + 1]);
		if (argv[a] == std::string("-libVariation")) libImageVariation = atof(argv[a + 1]);
	}

	// If either the input image or library images directory is not given, fail:
	if((inImgFilePath == std::string("None")) || (libDir == std::string("None")))
	{
		if(inImgFilePath == std::string("None")) std::cout << "ERROR. No image provided." << std::endl;
		if(libDir == std::string("None")) std::cout << "ERROR. No library image directory provided." << std::endl;

		std::cout << "\nUsage: mosaicTool -image <main image> -library <directory of library images> -numtiles <number of tiles in mosaic> -mosaicwidth <pixel width of output mosaic image> -pattern <tile pattern type> -cheat <cheat value> -libVariants <number> -libVariation <amount>" << std::endl;

		return 0;
	}

	// Generate the mosaic:
	run(inImgFilePath, libDir, numtiles, mosaicwidth, pattern, cheat, numLibImageVariants, libImageVariation);
}

