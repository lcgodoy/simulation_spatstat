Before compiling, you need to make sure you have a c and c++ compiler and the FFTW3 library installed on your computer. Then, go in the folder bin and then type:
cmake ../programs/

This will generate a makefile in the folder. You just need to type
make
to generate a binary file to run GcoPS.

To see what options are available, just type:
./GcoPS -h

For eaxample, if you just want to compute the p-value corresponding to the colocalization test between two segmented images input1.tif and input2.tif, just type:
./GcoPS -i1 input1.tif -i2 input2.tif
