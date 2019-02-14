# Images for cartopy

Note that the image [!50-natural-earth-1-downsampled.png ""](50-natural-earth-1-downsampled.png) was copied from cartopy.

As noted in the [images.json](images.json) file, this comes from [natural earth data](http://www.naturalearthdata.com/downloads/50m-raster-data/50m-natural-earth-1/).

I created the grayscale version of this using imagemagik:

```bash
convert 50-natural-earth-1-downsampled.png -set colorspace Gray -separate -average natural-earth-grayscale.png
```

This gives me the [!grayscale image](natural-earth-grayscale.png) that we are trying in the map, and then as [Nikolay Koldunov described](http://earthpy.org/tag/cartopy.html), I append information about this to the  [images.json](images.json) file so I can import it into my cartopy code.

```json

{"__comment__": "JSON file specifying the image to use for a given type/name and resolution. Read in by cartopy.mpl.geoaxes.read_user_background_images.",
  "ne_shaded": {
    "__comment__": "Natural Earth shaded relief. This is the image from cartopy",
    "__source__": "http://www.naturalearthdata.com/downloads/50m-raster-data/50m-natural-earth-1/",
    "__projection__": "PlateCarree",
    "low": "50-natural-earth-1-downsampled.png"  }

  "grayscale_shaded": {
    "__comment__": "Greyscale shaded relief",
    "__source__": "Converted by Rob Edwards using image magik from http://www.naturalearthdata.com/downloads/50m-raster-data/50m-natural-earth-1/",
    "__projection__": "PlateCarree",
    "low": "natural-earth-grayscale.png" }

}
```

Note that I left the original image in there.

Now I can edit my cartopy code, and add the background image.

```python

import os
os.environ["CARTOPY_USER_BACKGROUNDS"] = "/home/redwards/GitHubs/crAssphage/bin/map_drawing/crassphage_maps/images"

import matplotlib.pyplot as plt

ax = plt.axes(projection=ccrs.Robinson())
ax.background_img(name='grayscale_shaded', resolution='low')
```


