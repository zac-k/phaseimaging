### *class* phaseimaging.Image

| Attributes |  |  
|---|---|
| `resolution` | Dimensions of the image in pixels. |
| `width` | Dimensions of the image in units of length. |
| `image` | Numpy array containing the image data. |

| Methods |  |  
|---|---|
| `plot([limits])` | Generates an on-screen image using matplotlib. |
| `apodise([rad_sup])` | Apodises the image with a circular window function of radius rad_sup (as a fraction of the image width). |
