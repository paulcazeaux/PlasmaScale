Some newer versions of MinGW-g++ compiler have problem 
working with *.o codes(and hence static libraries) compiled in older versions of MinGW (3.x). If you are facing this problem try MINGW4.5 compiled static libraries or build wavelet2d from the source code with your own compiler. You'll need to statically build fftw3 library first if you are going the source code way. The DLLs should be working fine in any case so that's always another option.


---------------------------------------------------
1. The package includes release and debug versions.

2. Both files statically link to FFTW-3.3 library.

3. Usage is straightforward. Make sure that wavelet2s.h is included in your program and library/headers are properly
linked.