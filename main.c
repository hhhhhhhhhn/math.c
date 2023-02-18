#ifndef MATHC
#define MATHC

#ifdef MATHC_RUN
#include <stdio.h>
#include <stdlib.h>
#endif

#ifdef MATHC_NO_MATRIX
#define MATHC_MALLOC(x) NULL
#define MATHC_FREE(x) NULL
#endif

#if defined(_STDLIB_H) && !defined(MATHC_MALLOC)
#define MATHC_MALLOC malloc
#define MATHC_FREE free
#endif
#if !defined(MATHC_MALLOC) || !defined(MATHC_FREE)
#error "You must include stdlib.h before math.c, set a custom malloc with MATHC_MALLOC and MATHC_FREE, or define MATHC_NO_MATRIX"
#endif

const unsigned int SIGN_MASK = 0b10000000000000000000000000000000;
const unsigned int EXP_MASK  = 0b01111111100000000000000000000000;
const unsigned int MAN_MASK  = 0b00000000011111111111111111111111;
const float PI = 3.14159265358979323;
const float EULER = 2.718281828459045;

// Get exponent part
static inline int math_fexp(float x) {
	unsigned int* bits = (unsigned int*) &x;
	int exponent = (*bits & EXP_MASK) >> 23;
	return exponent - 127;
}

static inline void math_setfexp(float* x, unsigned int exponent) {
	exponent = exponent + 127;
	unsigned int* bits = (unsigned int*) x;
	*bits = (*bits & ~EXP_MASK) | (exponent << 23);
}

static inline void math_setfmantissa(float* x, unsigned int mantissa) {
	unsigned int* bits = (unsigned int*) x;
	*bits = (*bits & ~MAN_MASK) | mantissa;
}

// f = (1 + M/2^23)*2^(E-127)
// bin f = E*2^23 + M
// 
// log(x+ε) ≈ ε + µ
// µ = correction term, 0.05
// 
// log2 f = log2(1 + M/2^23) + E - 127
//        ≈ M/2^23 + E + µ - 127
//        ≈ (M + 2^23*E)/2^23 - 126.95
//        ≈ (bin f)/2^23 - 126.95
// 
// log2(sqrt(x)) = 1/2 log2(x)
// log2(s) = 1/2 log2(x)
// (bin s)/(2^23) - 126.95 = (bin x)/(2^24) - 126.95/2
// (bin s)/(2^23) = (bin x)/(2^24) - 126.95/2 + 126.95
// (bin s)/(2^23) = (bin x)/(2^24) + 126.95/2
// bin s = (bin x)/2 + 2^22*126.95
// bin s = (bin x)/2 + 532466892.8
// bin s = (bin x)/2 + 532466893
static inline float math_sqrt_initial(float x) {
	unsigned int* bits = (unsigned int*) &x;
	*bits = (*bits / 2) + 532466893;
	return x;
}

float math_sqrt(float x) {
	float initial = math_sqrt_initial(x);
	float result = 0.5*(initial + x/initial); // Newton's method
	result = 0.5*(result + x/result); // Newton's method
	return result;
}

static inline float math_abs(float x) {
	unsigned int* bits = (unsigned int*) &x;
	*bits = *bits & ~SIGN_MASK;
	return x;
}

float math_log2(float x) {
	float log = math_fexp(x);
	float s = x;
	math_setfexp(&s, 0);
	s = s - 1.5;
	float mantissa_factor = 0.58496+0.96179*s-0.32059*s*s+0.14248*s*s*s;
	return log + mantissa_factor;
}

int math_floor(float x) {
	return (int)x - (x < 0);
}

float math_exp2(float x) {
	int whole = math_floor(x);
	float decimal = x - whole;
	float s = decimal - 0.5;
	float decimal_factor = 1.41421 + 0.980258*s + 0.339732*s*s + 0.0784947*s*s*s;

	int exp = math_fexp(decimal_factor);
	exp += whole;
	math_setfexp(&decimal_factor, exp);

	return decimal_factor;
}

float math_pow(float base, float exponent) {
	return math_exp2(math_log2(base)*exponent);
}

float math_exp(float exponent) {
	return math_exp2(1.44269504089 * exponent);
}

// NOTE: Slower, less accurate, but funner than normal div
// f = (1 + M/2^23)*2^(E-127)
// bin f = E*2^23 + M
// 
// log(x+ε) ≈ ε + µ
// µ = correction term, 0.05
// 
// log2 f = log2(1 + M/2^23) + E - 127
//        ≈ M/2^23 + E + µ - 127
//        ≈ (M + 2^23*E)/2^23 - 126.95
//        ≈ (bin f)/2^23 - 126.95
// 
// log2(x/y) = log2(x) - log2(y)
// (bin s)/2^23 - 126.95 = (bin x)/2^23 - 126.95 - ((bin y)/2^23 - 126.95)
// (bin s)/2^23 - 126.95 = (bin x)/2^23 - (bin y)/2^23
// (bin s)/2^23 - 126.95 = (bin x - bin y)/2^23
// (bin s)/2^23 = (bin x - bin y)/2^23 + 126.95
// (bin s) = bin x - bin y + 126.95*2^23
// (bin s) = bin x - bin y + 126.95*2^23
// (bin s) = bin x - bin y + 1064933785.6
// (bin s) = bin x - bin y + 1064933786
float math_slow_div(float x, float y) {
	unsigned int* xbits = (unsigned int*) &x;
	*xbits = *xbits - (*(unsigned int*)&y) + 1064933786;
	return x;
}

float math_rem(float dividend, float divisor) {
	return dividend - divisor*((int)(dividend / divisor));
}

float math_mod(float dividend, float divisor) {
	float rem = math_rem(dividend, divisor);
	return rem + (rem < 0)*divisor;
}

// Bhaskara I
float math_sin(float x) {
	x = math_mod(x, 2*PI);
	int negate = x > PI;
	x = x - (x > PI)*PI;
	return (!negate*2 - 1) * 16*x*(PI - x)/(5*PI*PI - 4*x*(PI - x));
}

float math_cos(float x) {
	return math_sin(PI/2 - x);
}

float math_tan(float x) {
	return math_sin(x)/math_sin(PI/2 - x);
}

typedef struct math_matrix {
	int rows;
	int cols;
	float* data;
} math_matrix;

static inline int MATRIX_INDEX(math_matrix mat, int row, int col) {
	return row*mat.cols + col;
}
static inline float math_matrix_at(math_matrix mat, int row, int col) {
	return mat.data[MATRIX_INDEX(mat, row, col)];
}

math_matrix math_create_matrix(int rows, int cols, float* data) {
	math_matrix m = {rows, cols, data};
	return m;
}

// NOTE: Resulting matrix MUST BE FREED
math_matrix math_matrix_add(math_matrix a, math_matrix b) {
	float* data = MATHC_MALLOC(a.rows*a.cols*sizeof(float));
	math_matrix result = math_create_matrix(a.rows, a.cols, data);
	for (int row = 0; row < a.rows; row++) {
		for (int col = 0; col < a.cols; col++) {
			data[MATRIX_INDEX(result, row, col)] = math_matrix_at(a, row, col) + math_matrix_at(b, row, col);
		}
	}
	return result;
}

// NOTE: Resulting matrix MUST BE FREED
math_matrix math_matrix_scale(math_matrix mat, float f) {
	float* data = MATHC_MALLOC(mat.rows*mat.cols*sizeof(float));
	math_matrix result = math_create_matrix(mat.rows, mat.cols, data);
	for (int row = 0; row < mat.rows; row++) {
		for (int col = 0; col < mat.cols; col++) {
			data[MATRIX_INDEX(result, row, col)] = math_matrix_at(mat, row, col)*f;
		}
	}
	return result;
}

// NOTE: Resulting matrix MUST BE FREED
math_matrix math_matrix_transpose(math_matrix mat) {
	float* data = MATHC_MALLOC(mat.rows*mat.cols*sizeof(float));
	math_matrix result = math_create_matrix(mat.cols, mat.rows, data);
	for (int row = 0; row < mat.rows; row++) {
		for (int col = 0; col < mat.cols; col++) {
			data[MATRIX_INDEX(result, col, row)] = math_matrix_at(mat, row, col);
		}
	}
	return result;
}

// NOTE: Resulting matrix MUST BE FREED
math_matrix math_matrix_multiply(math_matrix a, math_matrix b) {
	int n = a.cols; // Common side
	float* data = MATHC_MALLOC(a.rows*b.cols*sizeof(float));

	math_matrix result = math_create_matrix(a.rows, b.cols, data);
	for (int row = 0; row < a.rows; row++) {
		for (int col = 0; col < b.cols; col++) {
			float sum = 0;
			for(int i = 0; i < n; i++) {
				sum += math_matrix_at(a, row, i)*math_matrix_at(b, i, col);
			}
			data[MATRIX_INDEX(result, row, col)] = sum;
		}
	}
	return result;
}

// NOTE: Resulting matrix MUST BE FREED
math_matrix math_submatrix(math_matrix mat, int removedRow, int removedCol) {
	float* data = MATHC_MALLOC((mat.rows-1)*(mat.cols-1)*sizeof(float));
	math_matrix result = math_create_matrix(mat.rows-1, mat.cols-1, data);

	for (int row = 0; row < result.rows; row++) {
		for (int col = 0; col < result.cols; col++) {
			int mat_col = col + (col >= removedCol);
			int mat_row = row + (row >= removedRow);
			data[MATRIX_INDEX(result, row, col)] = math_matrix_at(mat, mat_row, mat_col);
		}
	}
	return result;
}

float math_matrix_determinant(math_matrix mat) {
	if (mat.rows == 1) {
		return math_matrix_at(mat, 0, 0);
	}
	if (mat.rows == 2) {
		return math_matrix_at(mat, 0, 0)*math_matrix_at(mat, 1, 1)
			 - math_matrix_at(mat, 0, 1)*math_matrix_at(mat, 1, 0);
	}
	float sign = 1;
	float determinant = 0;
	for (int col = 0; col < mat.cols; col++) {
		math_matrix submat = math_submatrix(mat, 0, col);
		determinant += sign*math_matrix_at(mat, 0, col)*math_matrix_determinant(submat);
		MATHC_FREE(submat.data);
		sign = -sign;
	}
	return determinant;
}

float math_dot(math_matrix a, math_matrix b) {
	float dot = 0;
	for (int i = 0; i < a.rows*a.cols; i++) {
		dot += a.data[i]*b.data[i];
	}
	return dot;
}

// NOTE: Resulting matrix MUST BE FREED
// Returns a column vector
math_matrix math_cross(math_matrix a, math_matrix b) {
	float* data = MATHC_MALLOC(3*sizeof(float));
	math_matrix result = math_create_matrix(3, 1, data);
	result.data[0] = a.data[1]*b.data[2] - a.data[2]*b.data[1];
	result.data[1] = a.data[2]*b.data[0] - a.data[0]*b.data[2];
	result.data[2] = a.data[0]*b.data[1] - a.data[1]*b.data[0];
	return result;
}

// NOTE: Resulting matrix MUST BE FREED
math_matrix math_matrix_clone(math_matrix mat) {
	float* data = MATHC_MALLOC(mat.rows*mat.cols*sizeof(float));
	math_matrix result = math_create_matrix(mat.rows, mat.cols, data);
	for (int i = 0; i < mat.rows*mat.cols; i++) {
		data[i] = mat.data[i];
	}
	return result;
}

// NOTE: Done in place
void math_matrix_swap_rows(math_matrix mat, int row1, int row2) {
	if (row1 == row2) return;
	for (int col = 0; col < mat.cols; col++) {
		float temp = math_matrix_at(mat, row1, col);
		mat.data[MATRIX_INDEX(mat, row1, col)] = math_matrix_at(mat, row2, col);
		mat.data[MATRIX_INDEX(mat, row2, col)] = temp;
	}
}

// NOTE: Done in place
// Adds the (source row)*factor to the dest row
void math_matrix_add_row_with_factor(math_matrix mat, int src, int dest, float factor) {
	for (int col = 0; col < mat.cols; col++) {
		mat.data[MATRIX_INDEX(mat, dest, col)] += math_matrix_at(mat, src, col)*factor;
	}
}

// NOTE: Done in place
void math_matrix_scale_row(math_matrix mat, int row, float factor) {
	for (int col = 0; col < mat.cols; col++) {
		mat.data[MATRIX_INDEX(mat, row, col)] *= factor;
	}
}

// NOTE: Resulting matrix MUST BE FREED
math_matrix math_gauss_jordan(math_matrix mat) {
	mat = math_matrix_clone(mat);

	int pcol = 0; // Column of the pivot
	int prow = 0; // Row of the pivot
	while(prow < mat.rows && pcol < mat.cols) {
		// First, make sure the row has the earliest pivot
		for(int swap = prow; swap < mat.rows; swap++) {
			if (math_abs(math_matrix_at(mat, prow, pcol)) > 0.00001) {
				math_matrix_swap_rows(mat, prow, swap);
				break;
			}
		}
		float pivot = math_matrix_at(mat, prow, pcol);
		// No pivot found at this column, go to next
		if (math_abs(pivot) <= 0.00001) {
			pcol++;
			continue;
		}

		// Scale so the pivot is 1
		math_matrix_scale_row(mat, prow, 1/math_matrix_at(mat, prow, pcol));
		
		// Make every other element of the column a 0
		for (int target = 0; target < mat.rows; target++) {
			if (target == prow) continue;
			float factor = -math_matrix_at(mat, target, pcol); // divided by pivot == 1
			math_matrix_add_row_with_factor(mat, prow, target, factor);
		}
		prow++;
		pcol++;
	}
	return mat;
}


#ifdef _STDIO_H
int format_matrix(char* dest, int max_length, math_matrix matrix) {
	int i = 0;
	for (int row = 0; row < matrix.rows; row++) {
		for (int col = 0; col < matrix.cols; col++) {
			 i += snprintf(dest+i, max_length-i, "%10.1f", math_matrix_at(matrix, row, col));
		}
		i += snprintf(dest+i, max_length-i, "\n");
	}
	return i;
}
#endif

#endif

//////////////////// RUN SECTION ////////////////////////////////////
#ifdef MATHC_RUN

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
void keep_alive(volatile float result) {}

void bench() {
	int iters = 100000000;
	long start = clock();
	for (int i = 1; i < iters; i++) {
		float result = (float)1000000000/(float)i;
		keep_alive(result);
	}
	long end = clock();
	printf("DIV SPEED: %e iters per second\n", iters*CLOCKS_PER_SEC/(float)(end - start));

	start = clock();
	for (int i = 1; i < iters; i++) {
		float result = math_slow_div((float)1000000000, (float)i);
		keep_alive(result);
	}
	end = clock();
	printf("MATH_DIV SPEED: %e iters per second\n", iters*CLOCKS_PER_SEC/(float)(end - start));

	start = clock();
	for (int i = 1; i < iters; i++) {
		float result = sqrtf((float)i);
		keep_alive(result);
	}
	end = clock();
	printf("SQRT SPEED: %e iters per second\n", iters*CLOCKS_PER_SEC/(float)(end - start));
	
	start = clock();
	for (int i = 1; i < iters; i++) {
		float result = math_sqrt((float)i);
		keep_alive(result);
	}
	end = clock();
	printf("MATH_SQRT SPEED: %e iters per second\n", iters*CLOCKS_PER_SEC/(float)(end - start));
}

float sample_values[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

int main() {
	printf("22/7     = %f\n", math_slow_div(22, 7));
	printf("abs 2    = %f\n", math_abs(2));
	printf("abs -3   = %f\n", math_abs(-3));
	printf("sqrt 2   = %f\n", math_sqrt(2));
	printf("7 mod 3  = %f\n", math_mod(7, 3));
	printf("-1 mod 3 = %f\n", math_mod(-1, 3));
	printf("sin 4    = %f\n", math_sin(4));
	printf("sin 7    = %f\n", math_sin(7));
	printf("cos 9    = %f\n", math_cos(9));
	printf("tan 9    = %f\n", math_tan(9));
	printf("log2 243 = %f\n", math_log2(243));
	printf("⌊24.5⌋   = %i\n", math_floor(24.5));
	printf("⌊-24.5⌋  = %i\n", math_floor(-24.5));
	printf("2^10.7   = %f\n", math_exp2(10.7));
	printf("2^-10.7  = %f\n", math_exp2(-10.7));
	printf("34^9     = %e\n", math_pow(34, 9));
	printf("34^-9    = %e\n", math_pow(34, -9));
	printf("e^5.6    = %e\n", math_exp(5.6));
	printf("\n");
	char* str = MATHC_MALLOC(1000 * sizeof(char));
	math_matrix matrix = math_create_matrix(4, 3, sample_values);
	math_matrix small  = math_create_matrix(3, 2, sample_values);
	math_matrix square = math_create_matrix(3, 3, sample_values);
	math_matrix column = math_create_matrix(4, 1, sample_values);
	math_matrix row    = math_create_matrix(1, 4, sample_values+4);

	format_matrix(str, 1000, matrix);
	printf("mat:\n%s\n", str);

	format_matrix(str, 1000, small);
	printf("small:\n%s\n", str);

	format_matrix(str, 1000, square);
	printf("square:\n%s\n", str);

	format_matrix(str, 1000, column);
	printf("column:\n%s\n", str);

	format_matrix(str, 1000, row);
	printf("row:\n%s\n", str);

	math_matrix byTwo = math_matrix_add(matrix, matrix);
	format_matrix(str, 1000, byTwo);
	printf("mat+mat:\n%s\n", str);

	math_matrix half = math_matrix_scale(matrix, 0.5);
	format_matrix(str, 1000, half);
	printf("0.5*mat:\n%s\n", str);

	math_matrix transposed = math_matrix_transpose(matrix);
	format_matrix(str, 1000, transposed);
	printf("t(mat):\n%s\n", str);

	math_matrix timest = math_matrix_multiply(matrix, transposed);
	format_matrix(str, 1000, timest);
	printf("mat*t(mat):\n%s\n", str);

	math_matrix times_small = math_matrix_multiply(matrix, small);
	format_matrix(str, 1000, times_small);
	printf("mat*small:\n%s\n", str);

	math_matrix sub22 = math_submatrix(matrix, 2, 2);
	format_matrix(str, 1000, sub22);
	printf("sub(mat, 2, 2):\n%s\n", str);

	math_matrix cross = math_cross(column, row);
	format_matrix(str, 1000, cross);
	printf("column ⨯ row:\n%s\n", str);

	math_matrix clone = math_matrix_clone(matrix);
	format_matrix(str, 1000, clone);
	printf("clone(mat):\n%s\n", str);

	math_matrix_swap_rows(clone, 1, 3);
	format_matrix(str, 1000, clone);
	printf("swap(mat, 1, 3):\n%s\n", str);

	math_matrix gj = math_gauss_jordan(transposed);
	format_matrix(str, 1000, gj);
	printf("gaussjordan(t(mat)):\n%s\n", str);

	printf("det(square) = %f\n", math_matrix_determinant(square));

	printf("row·column = %f\n", math_dot(row, column));

	printf("\n");

	free(byTwo.data);
	free(half.data);
	free(transposed.data);
	free(timest.data);
	free(gj.data);
	free(sub22.data);
	free(times_small.data);
	free(clone.data);
	free(cross.data);
	free(str);

	bench();
	return 0;
}

#endif
