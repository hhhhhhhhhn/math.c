#include <stdio.h>
#ifndef MATHC
#define MATHC

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

float math_abs(float x) {
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

#endif

#ifdef MATH_RUN

#include <time.h>
#include <math.h>
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
	bench();
	return 0;
}

#endif
