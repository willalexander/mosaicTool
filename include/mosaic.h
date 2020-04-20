#include <string>

#define REGULAR_PATTERN_TYPE 0																				// Rectangular tiles
#define VORONOI_PATTERN_TYPE 1																				// Cellular tiles, uniformly distributed sites
#define VARYING_VORONOI_PATTERN_TYPE 2																		// Cellular tiles, sites distributed proportionally to variance

void run(std::string, std::string, int, int, int, float, int, float);										// Runs the mosaic generation process