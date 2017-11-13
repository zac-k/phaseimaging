### phaseimaging.project_electrostatic_phase

Generates a phase map by projecting the electrostatic potential of the specimen in the z-direction.

    phaseimaging.project_electrostatic_phase(specimen, accel_volt, mean_inner_potential, image_width)

|  |  |  |
|---|---|---|
| Parameters: | **specimen** : *array_like* |  |
|  |  | Three-dimensional binary specimen mask. |
|  | **accel_volt** : *float* |  |
|  |  | Electron accelerating voltage in volts. |
|  | **mean_inner_potential** : *complex float* |  |
|  |  | Constant mean inner electrostatic potential of the specimen. An imaginary component can be added for simulation of absorption. |
|  | **image_width** : *tuple*, *list* |  |
|  |  | Three element tuple or list containing the width of the specimen array in metres, in the x-, y-, and z-direction, respectively. |
| Returns: | **phase** : *ndarray* |  |
|  |  | Real, two-dimensional numpy array. The computed phase map. |