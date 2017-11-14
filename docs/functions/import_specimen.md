### phaseimaging.import_specimen
Generates a 3D numpy array from a specimen file. Input format is a text file (sans the file extension) with each block of rows/columns representing a 2D slice of the binary specimen mask. Each block is separated by an empty line. Specimen files can be generated using my [random-specimen-generator](https://github.com/zac-k/random-specimen-generator) repository. One specimen file is included in the present repository as an example.

    phaseimaging.import_specimen(specimen_file)

|  |  |  |
|---|---|---|
| Parameters: | **specimen_file**:*string* |  |
|  |  | Path of specimen file |
| Returns: | **specimen**:*ndarray* |  |
|  |  | Three-dimensional binary specimen mask. A value of `1` indicates a voxel where the specimen exists, and a value of `0` indicates a voxel where it does not. |
|  |  |  |

Note: Does not yet support rectangular (non-cubic) arrays.

[Back to Functions](functions.md)