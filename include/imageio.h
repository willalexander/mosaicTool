void imageio_init();																		// Any initialisation that may be required
float *readImage(std::string, int *, int *);												// Reads an image from disk, returns data as a float array
void writeImage(float *, std::string, int, int);											// Writes an image in float array format to disk
