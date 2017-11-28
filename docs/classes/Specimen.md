### *class* phaseimaging.Specimen
   
| Attributes |  |  
|---|---|
| `resolution` | Dimensions of the specimen array in pixels. |
| `width` | Dimensions of the specimen array in units of length. |
| `image` | Numpy array containing the specimen. |
| `mean_inner_potential` | Mean inner potential of the specimen. Complex values can be used to simulate attenuation. |
| `magnetisation` | Volume magnetisation. |
| `mhat` | Vector pointing in the magnetisation direction. The length of this vector has no effect. |
| `moment` | Three dimensional array of vectors representing the magnetisation direction at all points in space. This is used for arbitrary magnetisation configurations. `mhat` will be ignored if `moment` is not `None`. The magnitude of the vectors is in units of `specimen.magnetisation`.

| Methods |  |  
|---|---|
| `plot([limits])` | Generates an on-screen image, from a central slice, using matplotlib. |
| `apodise([rad_sup])` | Apodises the array with a cylindrical window function of radius rad_sup (as a fraction of the array width). |
| `mask_moment()` | Sets the magnetisation direction vector to `(0, 0, 0)` outside the specimen. This is performed automatically during initialisation.

[Back to Classes](classes.md)