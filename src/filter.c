#include <R.h>
#include <Rdefines.h>

// filterMA
//
// Filter [time-]series by replacing all values by the moving average of values
// centered around current one. Border values are averaged with available data.
//
// @param M_ A real matrix of size LxD
// @param w_ The (odd) number of values to average
//
// @return The filtered matrix, of same size as the input
SEXP filterMA(SEXP M_, SEXP w_)
{
	int L = INTEGER(getAttrib(M_, R_DimSymbol))[0],
		D = INTEGER(getAttrib(M_, R_DimSymbol))[1],
		w = INTEGER_VALUE(w_),
		half_w = (w-1)/2,
		i,
		nt; //number of terms in the current sum (max: w)
	double *M = REAL(M_),
		cs, //current sum (max: w terms)
		left; //leftward term in the current moving sum

	SEXP fM_; //the filtered result
	PROTECT(fM_ = allocMatrix(REALSXP, L, D));
	double* fM = REAL(fM_); //pointer to the encapsulated vector

	// NOTE: unused loop parameter: shifting at the end of the loop is more efficient
	for (int col=D-1; col>=0; col--)
	{
		// Initialization
		nt = half_w + 1;
		left = M[0];
		cs = 0.;
		for (i=half_w; i>=0; i--)
			cs += M[i];

		// Left border
		for (i=0; i<half_w; i++)
		{
			fM[i] = cs / nt; //(partial) moving average at current index i
			cs += M[i+half_w+1];
			nt++;
		}

		// Middle: now nt == w, i == half_w
		for (; i<L-half_w-1; i++)
		{
			fM[i] = cs / w; //moving average at current index i
			cs = cs - left + M[i+half_w+1]; //remove oldest items, add next
			left = M[i-half_w+1]; //update first value for next iteration
		}

		// (Last "fully averaged" index +) Right border
		for (; i<L; i++)
		{
			fM[i] = cs / nt; //(partial) moving average at current index i
			cs -= M[i-half_w];
			nt--;
		}

		// Shift by L == increment column index by 1 (storage per column)
		M = M + L;
		fM = fM + L;
	}

	UNPROTECT(1);
	return fM_;
}
