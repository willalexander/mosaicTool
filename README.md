# mosaicTool
Creates photomosaics. Recreates an image by building it up out of many smaller images

<img src="/images/mosaicA.bmp" width=50%/><img src="/images/mosaicB.bmp" width=50%/>
<img src="/images/mosaicC.bmp" width=50%/>

### Supports:
* Area lights
* Sphere lights
* OBJ polygon geometry
* Ljcambert diffuse surfaces
* Specular Surfaces
* Homogeneous, isotropic participating media
* Homogoenous, isotropic subsurface scattering

### Dependencies:
* ImageMagick
* jc_voronoi (https://github.com/JCash/voronoi)

### Usage:
```
mosaicTool.exe -image main.jpg -library library -numtiles 2000 -mosaicwidth 1024 -pattern regular -cheat 0.5 -libVariants 8 -libVariation 0.5
```

Parameter | Meaning
----------|----------
image | The main mage for the mosaic
library | A directory of images that will be used as small tiles in the mosaic
numtiles | The (approx) number of tiles which will form the mosaic
mosaicwidth | The width (in pixels) of the mosaic
pattern | The pattern layout of tiles ('regular', 'voronoi' or 'varyingvoronoi'
cheat | How much to tint the tiles to make them match the colour of the main image
libVariants | How many variants of each library image to create (each variant is scaled, rotated & tinted differently)
libVariation | How much scale & tint adjustment to apply to each library variant
