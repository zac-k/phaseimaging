### *class* phaseimaging.Specimen
   
| Attributes |  |  
|---|---|
| `resolution` | Dimensions of the specimen array in pixels. |
| `width` | Dimensions of the specimen array in units of length. |
| `image` | Numpy array containing the specimen. |
| `mean_inner_potential` | Mean inner potential of the specimen. Complex values can be used to simulate attenuation. |
| `magnetisation` | Volume magnetisation. |
| `mhat` | Vector pointing in the magnetisation direction. The length of this vector has no effect. |

| Methods |  |  
|---|---|
| `plot([limits])` | Generates an on-screen image, from a central slice, using matplotlib. |
| `apodise([rad_sup])` | Apodises the with a cylindrical window function of radius rad_sup (as a fraction of the array width). |

[Back to Classes](classes.md)