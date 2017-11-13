### *class* phaseimaging.Wavefield
    
| Attributes |  |  
|---|---|
| `resolution` | Dimensions of the image in pixels. |
| `width` | Dimensions of the image in units of length. |
| `image` | Numpy array containing the wavefield. |
| `defocus` | Distance corresponding to the defocus of the wavefield. Typically zero. |

| Methods |  |  
|---|---|
| `plot([limits])` | Generates an on-screen image using matplotlib. |
| `apodise([rad_sup])` | Apodises the image with a circular window function of radius rad_sup (as a fraction of the image width). |

[Back to Classes](classes.md)