# phaseimaging





This library provides a variety of functions related to phase imaging. Additionally, classes and their associated methods are included for if an object-oriented approach is preferred. Some of the included functions and methods can be used to process experimental images, but I have also included functions to simulate these images. For experimental data, pre-processing will be required for meaningful results. In the future, I will add image registration and other preprocessing functions to enable perfoming all steps after image acquisition using only this library.

Note that although I have coded most of the functions to allow for rectangular images and specimen arrays, most of my testing has used square and cubic arrays, so there is a potential for bugs if you are using rectangular arrays. Please let me know if you experience any issues with any of the functions.

## Important

This library uses the [unwrap](https://pypi.python.org/pypi/unwrap) package for phase unwrapping, to convert the exit-surface wavefield---computed using the multislice method---into a relative phase. This package seems not to play nicely with a lot of python versions, but it---and everything else required for this library---is working in python 3.6.0.

See the [documentation](https://zac-k.github.io/phaseimaging) for details on usage.
