BEAM Support files Readme 3/26/2018

The files in this directory are tables of planar flux albedo (Apf) 
generated for the Mineart scattering law with k = 1.0, 1.05, and 1.15. 

For more background see:

"Rough surfaces: Is the dark stuff just shadow?"
Jeffrey N. Cuzzi, Lindsey B. Chambers, Amanda R. Hendrix
Volume 289, p. 281-294 (2017)

In the tables, each row is for a specific value of mu0 (from 0.0 to 1.0) 
and each column is a specific shadowing parameter value from 0.0 to 2.0.
We assume a planar flux albedo (also known as normal refectance) equal 
to 1.0.

The tables assume a Hapke-like shadowing function and were generated 
using the program planarflux_simpsons.f90 found in the source directory.
If you are interested simulating shadowing based on a simple exponential
(f(alpha) = exp(-S*alpha)) you will need to modify the program (line 112)
and regenerate the tables.

The files are:

apf_table_k115.dat - k=1.15
apf_table_k105.dat - k=1.05
apf_table_k1.dat - k=1.00




