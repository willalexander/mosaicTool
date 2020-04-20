#ifndef __TYPES__
#define __TYPES__

struct pixVal																		// rgb pixel colour value
{
	float c[3];
};

struct point																		// simple 2D point
{
	float x, y;
};
	
struct bboxf																		// Rectangular bounding box
{
	float i0, i1;
	float j0, j1;
};

struct tileInfo																		// Contains info about a library tile image
{
	int id;
	point centre;
};

#endif