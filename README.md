# JMultiSlice

JMultiSlice is C++ library code providing routines for fast multislice calculations
on multiple CPU and single GPU. Routine interfaces are delivered by classes.

A multislice algorithm [1] approximates the solution of the coupled system of
differential equations describing the dynamic electron diffraction (multiple
scattering) in a crystal for fast electrons. The approximation is based on
making discrete steps of the forward scattering along the main propagation
direction of the electrons. For each step, the 3d crystal scattering potential
is projected in a thin slice, assuming that the effect of each slice is a small
deflection of the electron. This allows to separate the scattering from the
propagation to the next slice for faster numerical calculation.


[1] J.M. Cowley and A.F. Moodie, Acta. Cryst. 10 (1957) p. 609.
    https://doi.org/10.1107/S0365110X57002194


Detailed reference to the approaches used in the code is given in [2].

[2] J. Barthel, Ultramic. 193 (2018) p. 1.
    https://doi.org/10.1016/j.ultramic.2018.06.003


## Usage

In your code, declare an object of CJMultiSlice and use its member functions to
setup calculation parameters.

Example, see 
[test/JMSBench1](https://github.com/ju-bar/JMultiSliceLib/tree/master/test/JMSBench1)


## Testing

The program JMSBench1 is available as project in the test sub folder. More tests have been
made with the class interface CJMultiSlice in the 
[Dr. Probe graphical user interface](http://www.er-c.org/barthel/drprobe/).


## Further development

Current version: 0.45

See [Changes.txt](https://github.com/ju-bar/JMultiSliceLib/blob/master/Changes.txt) for reading on recent modifications.

You are welcome to propose further changes. A list with possible future steps is maintained
in the file [ToDo.txt](https://github.com/ju-bar/JMultiSliceLib/blob/master/ToDo.txt).
