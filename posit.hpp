#include <ap_int.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <hls_math.h>
#include <vector>
#include "hls_stream.h"
using namespace std;

#define N 14
#define ES 1
#define TERMS 2
#define IN_SIZE 1024
#define APPR_TAILOR 1
#define FFT 1
#define IFFT 1
#define EXACT 0
#if (IN_SIZE ==32)
#define LOG_IN_SIZE_SF 5
#endif
#if (IN_SIZE ==64)
#define LOG_IN_SIZE_SF 6
#endif
#if (IN_SIZE ==128)
#define LOG_IN_SIZE_SF 7

#endif
#if (IN_SIZE ==256)
#define LOG_IN_SIZE_SF 8

#endif
#if (IN_SIZE ==512)
#define LOG_IN_SIZE_SF 9

#endif
#if (IN_SIZE ==1024)
#define LOG_IN_SIZE_SF 10

#endif
#if (IN_SIZE ==2048)
#define LOG_IN_SIZE_SF 11

#endif
#if (IN_SIZE ==4096)
#define LOG_IN_SIZE_SF 12

#endif
typedef ap_uint<LOG_IN_SIZE_SF> log_sf_t;
#define FRAC_LEN (N-(ES+2)+1)
#define MUL_LEN 2*FRAC_LEN

#if (ES==0)
#define REG_SHIFT 1 //(1<<ES)
#define USEED 2 //(1<<REG_SHIFT)
#endif
#if (ES==1)
#define REG_SHIFT 2 //(1<<ES)
#define USEED 4 //(1<<REG_SHIFT)
#endif
#if (ES==2)
#define REG_SHIFT 4 //(1<<ES)
#define USEED 16 //(1<<REG_SHIFT)
#endif
#if (ES==3)
#define REG_SHIFT 8 //(1<<ES)
#define USEED 256 //(1<<REG_SHIFT)
#endif
#if (N==13) 
#define REG_LEN 5 //log2(N)+1
#endif
#if (N==14) 
#define REG_LEN 5 //log2(N)+1
#endif
#if (N==15) 
#define REG_LEN 5 //log2(N)+1
#endif
#if (N==16) 
#define REG_LEN 5 //log2(N)+1
#endif
#if (N==20) 
#define REG_LEN 6 //log2(N)+1
#endif
#if (N==24) 
#define REG_LEN 6 //log2(N)+1
#endif
#if (N==28) 
#define REG_LEN 6 //log2(N)+1
#endif
#if (N==30) 
#define REG_LEN 6 //log2(N)+1
#endif
#if (N==32) 
#define REG_LEN 6 //log2(N)+1
#endif

#if (N==64) 
#define REG_LEN 7 //log2(N)+1
#endif
typedef ap_int<REG_LEN> regime_t;
typedef ap_uint<REG_LEN> ml_t;
typedef ap_uint<FRAC_LEN> mantissa_t;
typedef ap_uint<FRAC_LEN+REG_LEN+ES+1> mantissa_sf_t;
#if ES!=0
typedef ap_uint<ES> exponent_t;
#else
typedef ap_uint<1> exponent_t;
#endif
typedef ap_uint<FRAC_LEN+2> mant_add_t;
typedef ap_uint<FRAC_LEN+1> m_add_t;
typedef ap_uint<2> ovf_t;
typedef ap_uint<MUL_LEN> mul_t;
typedef ap_uint<MUL_LEN+2> dv_t;
typedef ap_int<REG_LEN+ES+1> sf_t;

typedef ap_uint<N-1> reg_t;
typedef ap_uint<N> posit_t;
typedef ap_uint<1> bool_t;


typedef struct POSIT{bool sign=0;bool isZero=1;bool isInf=0;regime_t regime=0;exponent_t exponent=0;mantissa_t mantissa=0;}POSIT;
typedef POSIT ps_t;

inline bool operator==(const POSIT& lhs, const POSIT& rhs) {
    return lhs.sign     == rhs.sign     &&
           lhs.isZero   == rhs.isZero   &&
           lhs.isInf    == rhs.isInf    &&
           lhs.regime   == rhs.regime   &&
           lhs.exponent == rhs.exponent &&
           lhs.mantissa == rhs.mantissa;
}



#if (N==32)
//#define PI_MANTISSA 843251712
#define PI_MANTISSA 1686110208
//#define PI_MANTISSA 838860800
#endif
#if (N==30)
//#define PI_MANTISSA 843251712
#define PI_MANTISSA 421527552
//#define PI_MANTISSA 838860800
#endif
#if (N==28)
//#define PI_MANTISSA 52703232
#define PI_MANTISSA 105381888
//#define PI_MANTISSA 52428800
#endif
#if (N==24)
//#define PI_MANTISSA 3293952
#define PI_MANTISSA 6586368
//#define PI_MANTISSA 3276800
#endif
#if (N==20)
//#define PI_MANTISSA 3293952
#define PI_MANTISSA 411648
//#define PI_MANTISSA 3276800
#endif
#if (N==16)
//#define PI_MANTISSA 3293952
#define PI_MANTISSA 25728
//#define PI_MANTISSA 3276800
#endif
#if (N==15)
//#define PI_MANTISSA 3293952
#define PI_MANTISSA 12864
//#define PI_MANTISSA 3276800
#endif
#if (N==14)
//#define PI_MANTISSA 3293952
#define PI_MANTISSA 6432
//#define PI_MANTISSA 3276800
#endif
#if (N==13)
//#define PI_MANTISSA 3293952
#define PI_MANTISSA 3216
//#define PI_MANTISSA 3276800
#endif
//#define PI 3.141357421875
#define PI 3.140625
//#define PI 3.125
#if (ES==0)

const ps_t POSIT_PI_OVER2 =     {0, false, false, 0, 0, PI_MANTISSA};
const ps_t POSIT_PI =           {0, false, false, 1, 0, PI_MANTISSA};
const ps_t POSIT_2PI =          {0, false, false, 2, 0, PI_MANTISSA}; 
const ps_t POSIT_M_PI_OVER2 =   {1, false, false, 0, 0, PI_MANTISSA};
const ps_t POSIT_M_PI =         {1, false, false, 1, 0, PI_MANTISSA};
const ps_t POSIT_M_2PI =        {1, false, false, 2, 0, PI_MANTISSA}; 
#endif 
#if (ES==1)
const ps_t POSIT_PI_OVER2 =     {0, false, false, 0, 0, PI_MANTISSA/2};
const ps_t POSIT_PI =           {0, false, false, 0, 1, PI_MANTISSA/2};
const ps_t POSIT_2PI =          {0, false, false, 1, 0, PI_MANTISSA/2}; 
const ps_t POSIT_M_PI_OVER2 =   {1, false, false, 0, 0, PI_MANTISSA/2};
const ps_t POSIT_M_PI =         {1, false, false, 0, 1, PI_MANTISSA/2};
const ps_t POSIT_M_2PI =        {1, false, false, 1, 0, PI_MANTISSA/2}; 
#endif 
#if  (ES==2)
//#define PI 3.14159265161
const ps_t POSIT_PI_OVER2 =     {0, false, false, 0, 0, PI_MANTISSA/4};
const ps_t POSIT_PI =           {0, false, false, 0, 1, PI_MANTISSA/4};
const ps_t POSIT_2PI =          {0, false, false, 0, 2, PI_MANTISSA/4};
const ps_t POSIT_M_PI_OVER2 =   {1, false, false, 0, 0, PI_MANTISSA/4};
const ps_t POSIT_M_PI =         {1, false, false, 0, 1, PI_MANTISSA/4};
const ps_t POSIT_M_2PI =        {1, false, false, 0, 2, PI_MANTISSA/4};
#endif

const ps_t ONE = {0, false, false, 0, 0, 1<<(FRAC_LEN-1)};
const ps_t mONE = {1, false, false, 0, 0, 1<<(FRAC_LEN-1)};
const ps_t ZERO = {0, true, false, 0, 0, 1<<(FRAC_LEN-1)}; 
ps_t calculateKFactor(int k);
regime_t LOD(reg_t reg);
int LOD_ADD(mant_add_t in);
ps_t  decode(posit_t posit);
posit_t encode(ps_t X);
posit_t posit_reverse(posit_t value);
ps_t posit_negate(ps_t posit);
ps_t int2posit(int X);
float posit2float(ps_t pos);
double posit2double(ps_t pos);
ps_t float2posit(float in);
ps_t double2posit(double in);
ps_t positMod(ps_t x, ps_t y);
ps_t pReduceAngle(ps_t angle, bool &negate);
ps_t positAdd(ps_t x,ps_t y);
ps_t positMul(ps_t x,ps_t y);
ps_t positDiv(ps_t x,ps_t y);
ps_t positDiv2p(ps_t in,int i);
ps_t positSub(ps_t x,ps_t y);
ps_t positCos(ps_t x);
ps_t positSin(ps_t x);
ps_t positMod(ps_t x, ps_t y);
ps_t normalize_angle(ps_t x) ;
double dTailorSin(double in);
float fTailorSin(float in);
half hTailorSin(half in);
double dTailorCos(double in);
float fTailorCos(float in);
half hTailorCos(half in);
void pEuler(ps_t angle, ps_t *result_real, ps_t *result_imag);
void dEuler(double angle, double *result_real, double *result_imag);
void fEuler(float angle, float *result_real, float *result_imag);
void hEuler(half angle, half *result_real, half *result_imag);

struct pFFTResult {
    ps_t real[IN_SIZE];  // Array for real values
    ps_t imag[IN_SIZE];  // Array for imaginary values

    // Constructor (initialize arrays)
    pFFTResult(size_t size = IN_SIZE) {
        // Arrays are statically sized, no need to resize
        std::fill(std::begin(real), std::end(real), ps_t{ZERO});  // Initialize with zero
        std::fill(std::begin(imag), std::end(imag), ps_t{ZERO});  // Initialize with zero
    }
};

struct dFFTResult {
    double real[IN_SIZE];  // Array for real values
    double imag[IN_SIZE];  // Array for imaginary values

    // Constructor (initialize arrays)
    dFFTResult(size_t size = IN_SIZE) {
        // Arrays are statically sized, no need to resize
        std::fill(std::begin(real), std::end(real), 0.0);  // Initialize with zero
        std::fill(std::begin(imag), std::end(imag), 0.0);  // Initialize with zero
    }
};

struct fFFTResult {
    float real[IN_SIZE];  // Array for real values
    float imag[IN_SIZE];  // Array for imaginary values

    // Constructor (initialize arrays)
    fFFTResult(size_t size = IN_SIZE) {
        // Arrays are statically sized, no need to resize
        std::fill(std::begin(real), std::end(real), 0.0f);  // Initialize with zero
        std::fill(std::begin(imag), std::end(imag), 0.0f);  // Initialize with zero
    }
};
struct hFFTResult {
    half real[IN_SIZE];  // Array for real values
    half imag[IN_SIZE];  // Array for imaginary values

    // Constructor (initialize arrays)
    hFFTResult(size_t size = IN_SIZE) {
        // Arrays are statically sized, no need to resize
        std::fill(std::begin(real), std::end(real), 0.0f);  // Initialize with zero
        std::fill(std::begin(imag), std::end(imag), 0.0f);  // Initialize with zero
    }
};
void fFFT(float signal[], fFFTResult& result);
void hFFT(half signal[], hFFTResult& result);
void dFFT(double signal[], dFFTResult& result);
void pFFT(ps_t signal[], pFFTResult& result);

// IFFT functions for different signal types
void dIFFT(dFFTResult& result, double signal[], int sampleCount);
void fIFFT(fFFTResult& result, float signal[], int sampleCount);
void hIFFT(hFFTResult& result, half signal[], int sampleCount);
void pIFFT(pFFTResult& result, ps_t signal[], int sampleCount);


#ifndef GLOBALS_H
#define GLOBALS_H

#include <cstddef>  

const size_t MAX_SIZE = IN_SIZE*IN_SIZE;
const size_t DTHETA_SIZE = IN_SIZE;
extern double d_Angle_Array[MAX_SIZE];   
extern size_t d_Angle_Counter;           
extern float f_Angle_Array[MAX_SIZE];   
extern size_t f_Angle_Counter;           
extern double p_Angle_Array[MAX_SIZE];   
extern size_t p_Angle_Counter;           

void addTod_Angle_Array(double value);
void f_addTof_Angle_Array(float value);
void p_addTop_Angle_Array(double value);
///RP
extern double d_RP_Array[MAX_SIZE];   
extern size_t d_RP_Counter;           
extern float f_RP_Array[MAX_SIZE];   
extern size_t f_RP_Counter;           
extern double p_RP_Array[MAX_SIZE];   
extern size_t p_RP_Counter; 

void addTod_RP_Array(double value);
void f_addTof_RP_Array(float value);
void p_addTop_RP_Array(double value);
///IMG
extern double d_IMG_Array[MAX_SIZE];   
extern size_t d_IMG_Counter;           
extern float f_IMG_Array[MAX_SIZE];   
extern size_t f_IMG_Counter;           
extern double p_IMG_Array[MAX_SIZE];   
extern size_t p_IMG_Counter; 

void addTod_IMG_Array(double value);
void f_addTof_IMG_Array(float value);
void p_addTop_IMG_Array(double value);
///DTheta
extern double d_DTH_Array[DTHETA_SIZE];   
extern size_t d_DTH_Counter;           
extern float f_DTH_Array[DTHETA_SIZE];   
extern size_t f_DTH_Counter;           
extern double p_DTH_Array[DTHETA_SIZE];   
extern size_t p_DTH_Counter; 

void addTod_DTH_Array(double value);
void f_addTof_DTH_Array(float value);
void p_addTop_DTH_Array(double value);
#endif