#pragma once
#ifdef FFTW
#include <fftw3.h>
#include <iostream>
namespace convolution_helpers
{
	void ComplexScale_fftw(fftw_complex, double);
	void ComplexMul_fftw(fftw_complex, const fftw_complex);
	void ComplexPointwiseMulAndScale_fftw(fftw_complex *, const fftw_complex *, int, double);
}
#endif

