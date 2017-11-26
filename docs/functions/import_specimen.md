### phaseimaging.import_specimen
Generates a 3D numpy array from a specimen file. Input format can be either of two formats: a text file with each block of rows/columns representing a 2D slice of the binary specimen mask, with each block separated by an empty line; a numpy array generated with `numpy.save()`. Any file with a `.npy` extension will be assumed to be in the latter format. Any other extension (or none) will be assumed to be in the former.

Specimen files can be generated using my [random-specimen-generator](https://github.com/zac-k/random-specimen-generator) repository. One specimen file is included in the present repository as an example.

    phaseimaging.import_specimen(specimen_file)

|  |  |  |
|---|---|---|
| Parameters: | **specimen_file**:*string* |  |
|  |  | Path of specimen file |
| Returns: | **specimen**:*ndarray* |  |
|  |  | Three-dimensional binary specimen mask. A value of `1` indicates a voxel where the specimen exists, and a value of `0` indicates a voxel where it does not. |
|  |  |  |

Note: Does not support rectangular (non-cubic) arrays when using the text file format.

[Back to Functions](functions.md)