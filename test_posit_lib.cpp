#include "posit.hpp"
#include <climits>
#include <etc/hls_sqrt_apfixed.h>
#include <hls_math.h>
#include <ostream>
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <filesystem>


void fillRandomDoublesArray(double* arr, int size, double minVal, double maxVal, unsigned int seed) {
    std::mt19937 gen(seed);  
    std::uniform_real_distribution<double> dist(minVal, maxVal); 

    for (int i = 0; i < size; ++i) {
        arr[i] = dist(gen);
    }
}
void convertPositArrayToDoubleArray(const ps_t* positArray, double* doubleArray, size_t size) {
    for (size_t i = 0; i < size; ++i) {
        ps_t positValue = positArray[i]; // Access array element
        double reconstructed = posit2double(positValue);  // Convert posit to double
        doubleArray[i] = reconstructed; // Store in the output array
    }
}

void convertDoubleArrayToPositArray(ps_t positArray[], double doubleArray[], size_t size) {
    for (size_t i = 0; i < size; ++i) {
        double doubleValue = doubleArray[i]; // Access array element
        
        ps_t converted = double2posit(doubleValue);  // Convert posit to double
        positArray[i] = converted; // Store in the output array
        if(i==1)
        std::cout<<doubleValue<<"--"<<posit2double(converted)<<std::endl;
    }
}
void convertDoubleArrayToFloatArray(float floatArray[], double doubleArray[], size_t size) {
    float converted=0;
    for (size_t i = 0; i < size; ++i) {
        double doubleValue = doubleArray[i]; 
        converted = (float)doubleValue;  
        floatArray[i] = converted; 
    }
}
void convertDoubleArrayToHalfArray(half halfArray[], double doubleArray[], size_t size) {
    half converted=0;
    for (size_t i = 0; i < size; ++i) {
        double doubleValue = doubleArray[i]; 
        converted = (half)doubleValue;  
        halfArray[i] = converted; 
    }
}
void calculateSinArrays(
    const float floatSinCosIn[], float SinFloat[],
    const double doubleSinCosIn[], double SinDouble[],
    const ps_t positSinCosIn[], ps_t SinPosit[],
    int size)
{
    for (int i = 0; i < size; ++i) {
        // Calculate sine for float
        SinFloat[i] = fTailorSin(floatSinCosIn[i]);

        // Calculate sine for double
        SinDouble[i] = dTailorSin(doubleSinCosIn[i]);

        // Calculate sine for posit
        SinPosit[i] = positSin(positSinCosIn[i]);  
    }
}
void calculateCosArrays(
    const float floatSinCosIn[], float CosFloat[],
    const double doubleSinCosIn[], double CosDouble[],
    const ps_t positSinCosIn[], ps_t CosPosit[],
    int size)
{
    for (int i = 0; i < size; ++i) {
        // Calculate sine for float
        CosFloat[i] = fTailorCos(floatSinCosIn[i]);

        // Calculate sine for double
        CosDouble[i] = dTailorCos(doubleSinCosIn[i]);

        // Calculate sine for posit
        CosPosit[i] = positCos(positSinCosIn[i]);  
    }
}
// Compute SNR between original and reconstructed array
double computeSNR(const double original[], const double reconstructed[], size_t size) {
    double signalPower = 0.0, noisePower = 0.0;

    for (size_t i = 0; i < size; ++i) {
        signalPower += original[i] * original[i];
        double error = original[i] - reconstructed[i];
        noisePower += error * error;
    }

    if (noisePower == 0.0) return std::numeric_limits<double>::infinity(); // Perfect reconstruction

    return 10.0 * log10(signalPower / noisePower);
}

// Compute SNR between original double array and float array
double computeSNRFloat(const double original[], const float reconstructed[], size_t size) {
    float signalPower = 0.0, noisePower = 0.0;

    for (size_t i = 0; i < size; ++i) {
        signalPower += original[i] * original[i];
        float error = original[i] - reconstructed[i];
        noisePower += error * error;
    }

    if (noisePower == 0.0) return std::numeric_limits<float>::infinity(); // Perfect reconstruction

    return 10.0 * log10(signalPower / noisePower);
}
double computeSNRHalf(const double original[], const half reconstructed[], size_t size) {
    double signalPower = 0.0, noisePower = 0.0;

    for (size_t i = 0; i < size; ++i) {
        signalPower += original[i] * original[i];
        double error = original[i] - (double)reconstructed[i];
        noisePower += error * error;
    }

    if (noisePower == 0.0) return std::numeric_limits<half>::infinity(); // Perfect reconstruction

    return 10.0 * log10(signalPower / noisePower);
}
void writeToFileFloat(const std::string& filename, const float data_array[], size_t size) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }

    // Loop through the array and write each value to the file
    for (size_t i = 0; i < size; ++i) {
        outFile << data_array[i] << "\n";
    }

    outFile.close();
    std::cout << "Saved " << size << " values to " << filename << std::endl;
}
void writeToFileHalf(const std::string& filename, const half data_array[], size_t size) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }

    // Loop through the array and write each value to the file
    for (size_t i = 0; i < size; ++i) {
        outFile << data_array[i] << "\n";
    }

    outFile.close();
    std::cout << "Saved " << size << " values to " << filename << std::endl;
}


void writeToFile(const std::string& filename, const double data_array[], size_t size) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }

    // Loop through the array and write each value to the file
    for (size_t i = 0; i < size; ++i) {
        outFile <<data_array[i] << "\n";
    }

    outFile.close();
    std::cout << "Saved " << size << " values to " << filename << std::endl;
}
void addDoubleArrays(const double* arr1, const double* arr2, double* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] + arr2[i];
    }
}
void addFloatArrays(const float* arr1, const float* arr2, float* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] + arr2[i];
    }
}
void addHalfArrays(const half* arr1, const half* arr2, half* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] + arr2[i];
    }
}
void addPositArrays(const ps_t* arr1, const ps_t* arr2, ps_t* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = positAdd(arr1[i] , arr2[i]);
    }
}
void subDoubleArrays(const double* arr1, const double* arr2, double* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] - arr2[i];
    }
}
void subFloatArrays(const float* arr1, const float* arr2, float* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] - arr2[i];
    }
}
void subHalfArrays(const half* arr1, const half* arr2, half* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] - arr2[i];
    }
}
void subPositArrays(const ps_t* arr1, const ps_t* arr2, ps_t* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = positSub(arr1[i] , arr2[i]);
    }
}
void mulDoubleArrays(const double* arr1, const double* arr2, double* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] * arr2[i];
    }
}
void mulFloatArrays(const float* arr1, const float* arr2, float* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] * arr2[i];
    }
}
void mulHalfArrays(const half* arr1, const half* arr2, half* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] * arr2[i];
    }
}
void mulPositArrays(const ps_t* arr1, const ps_t* arr2, ps_t* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = positMul(arr1[i] , arr2[i]);
    }
}
void divDoubleArrays(const double* arr1, const double* arr2, double* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] / arr2[i];
    }
}
void divFloatArrays(const float* arr1, const float* arr2, float* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] / arr2[i];
    }
}
void divHalfArrays(const half* arr1, const half* arr2, half* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] / arr2[i];
    }
}
void divPositArrays(const ps_t* arr1, const ps_t* arr2, ps_t* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = positDiv(arr1[i] , arr2[i]);
    }
}
void div2DoubleArrays(const double* arr1, int C, double* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] / C;
    }
}
void div2FloatArrays(const float* arr1, int C, float* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] / C;
    }
}
void div2HalfArrays(const half* arr1, int C, half* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = arr1[i] / C;
    }
}
void div2PositArrays(const ps_t* arr1, int C, ps_t* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = positDiv2p(arr1[i] , C);
    }
}
void modDoubleArrays(const double* arr1, double* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = fmod(arr1[i] , 2*PI);
    }
}
void modFloatArrays(const float* arr1, float* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = fmod(arr1[i] , 2*PI);
    }
}
void modHalfArrays(const half* arr1, half* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = hls::half_fmod(arr1[i] , 2*PI);
    }
}
void modPositArrays(const ps_t* arr1, ps_t* result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = positMod(arr1[i] ,POSIT_2PI);
    }
}
void printPosit(ps_t x){
    std::cout<<"sign: "<<x.sign<<" ";
    std::cout<<"isZero: "<<x.isZero<<" ";
    std::cout<<"isInf: "<<x.isInf<<" ";
    std::cout<<"exponent: "<<x.exponent<<" ";
    std::cout<<"regime: "<<x.regime<<" ";
    std::cout<<"mantissa: "<<x.mantissa<<" "<<std::endl;
}

int main() {
        // Construct the base path
    std::string basePath = "C:/Users/Burak/Desktop/PFFT/PFFT" ;
	std::string d_inFile = basePath + "/d_in.txt";  
    std::string f_inFile = basePath + "/f_in.txt";
	std::string h_inFile = basePath + "/h_in.txt";    
    std::string p_inFile = basePath + "/p_in.txt";
	std::string d_in2File = basePath + "/d_in2.txt";  
    std::string f_in2File = basePath + "/f_in2.txt"; 
	std::string h_in2File = basePath + "/h_in2.txt"; 
    std::string p_in2File = basePath + "/p_in2.txt";
	std::string d_SumFile = basePath + "/d_sum.txt";  
    std::string f_SumFile = basePath + "/f_sum.txt";
	std::string h_SumFile = basePath + "/h_sum.txt";  
    std::string p_SumFile = basePath + "/p_sum.txt";
	std::string d_SubFile = basePath + "/d_sub.txt";  
    std::string f_SubFile = basePath + "/f_sub.txt";
	std::string h_SubFile = basePath + "/h_sub.txt";  
    std::string p_SubFile = basePath + "/p_sub.txt";
    std::string d_ProdFile = basePath + "/d_prod.txt";  
    std::string f_ProdFile = basePath + "/f_prod.txt"; 
	std::string h_ProdFile = basePath + "/h_prod.txt";  
    std::string p_ProdFile = basePath + "/p_prod.txt";
    std::string d_DivFile = basePath + "/d_div.txt";  
    std::string f_DivFile = basePath + "/f_div.txt"; 
	std::string h_DivFile = basePath + "/h_div.txt";  
    std::string p_DivFile = basePath + "/p_div.txt";
    std::string outSNR = basePath + "/N_" + std::to_string(N) 
                            + "_ES_" + std::to_string(ES) +"/result.txt";
	std::ofstream outFile(outSNR);
	
    int TEST_SIZE = 1000;
    double double_in[TEST_SIZE];
    double double_in2[TEST_SIZE];
    float float_in[TEST_SIZE];
    float float_in2[TEST_SIZE];
    half half_in[TEST_SIZE];
    half half_in2[TEST_SIZE];
    ps_t posit_in[TEST_SIZE];
    ps_t posit_in2[TEST_SIZE];
    double d_posit_in[TEST_SIZE];
    double d_posit_in2[TEST_SIZE];


    ///////////////////////////////////////////////INPUT//////////////////////////////////////////
    fillRandomDoublesArray(double_in,TEST_SIZE,-1,1,26);
    fillRandomDoublesArray(double_in2,TEST_SIZE,-1,1,27);
    convertDoubleArrayToPositArray(posit_in,double_in,TEST_SIZE);
    convertPositArrayToDoubleArray(posit_in, d_posit_in, TEST_SIZE);
    convertDoubleArrayToPositArray(posit_in2,double_in2,TEST_SIZE);
    convertPositArrayToDoubleArray(posit_in2, d_posit_in2, TEST_SIZE);

    convertDoubleArrayToFloatArray(float_in,double_in,TEST_SIZE);
    convertDoubleArrayToFloatArray(float_in2,double_in2,TEST_SIZE);
    convertDoubleArrayToHalfArray(half_in,double_in,TEST_SIZE);
    convertDoubleArrayToHalfArray(half_in2,double_in2,TEST_SIZE);
    writeToFile(d_inFile, double_in, TEST_SIZE);
    writeToFile(d_in2File, double_in2, TEST_SIZE);
    writeToFile(p_inFile, d_posit_in, TEST_SIZE);
    writeToFile(p_in2File, d_posit_in2, TEST_SIZE);
    writeToFileFloat(f_inFile, float_in, TEST_SIZE);
    writeToFileFloat(f_in2File, float_in2, TEST_SIZE);
	writeToFileHalf(h_inFile, half_in, TEST_SIZE);
    writeToFileHalf(h_in2File, half_in2, TEST_SIZE);
	//////////////////////////INPUT//////////////////////////
	double SNR_IN_P = computeSNR(double_in,d_posit_in,TEST_SIZE);
    double SNR_IN_F = computeSNRFloat(double_in, float_in, TEST_SIZE);
	double SNR_IN_H = computeSNRHalf(double_in, half_in, TEST_SIZE);
    std::cout<<"SNR IN POSIT VS DOUBLE :"<<SNR_IN_P<<std::endl;
    std::cout<<"SNR IN FLOAT VS DOUBLE :"<<SNR_IN_F<<std::endl;
    std::cout<<"SNR IN HALF VS DOUBLE :"<<SNR_IN_H<<std::endl;		
    //////////////////////////TEST ADD//////////////////////////
    std::cout<<"--------------------------------ADD------------------------"<<std::endl;
    double d_sum[TEST_SIZE];
    float f_sum[TEST_SIZE];
	half h_sum[TEST_SIZE];
    ps_t p_sum[TEST_SIZE];
    double d_p_sum[TEST_SIZE];
    addDoubleArrays(double_in, double_in2, d_sum, TEST_SIZE);
    addFloatArrays(float_in, float_in2, f_sum, TEST_SIZE);
	addHalfArrays(half_in, half_in2, h_sum, TEST_SIZE);
    addPositArrays(posit_in, posit_in2, p_sum, TEST_SIZE);
    convertPositArrayToDoubleArray(p_sum, d_p_sum, TEST_SIZE);
    double SNR_ADD_P = computeSNR(d_sum,d_p_sum,TEST_SIZE);
    double SNR_ADD_F = computeSNRFloat(d_sum, f_sum, TEST_SIZE);
	double SNR_ADD_H = computeSNRHalf(d_sum, h_sum, TEST_SIZE);
    std::cout<<"SNR ADD POSIT VS DOUBLE :"<<SNR_ADD_P<<std::endl;
    std::cout<<"SNR ADD FLOAT VS DOUBLE :"<<SNR_ADD_F<<std::endl;
	std::cout<<"SNR ADD HALF VS DOUBLE :"<<SNR_ADD_H<<std::endl;
  /*  writeToFile(d_SumFile, d_sum, TEST_SIZE);
    writeToFile(p_SumFile, d_p_sum, TEST_SIZE);
    writeToFileFloat(f_SumFile, f_sum, TEST_SIZE);*/
    //////////////////////////TEST SUB//////////////////////////
    std::cout<<"--------------------------------SUB------------------------"<<std::endl;
    double d_sub[TEST_SIZE];
    float f_sub[TEST_SIZE];
	half h_sub[TEST_SIZE];
    ps_t p_sub[TEST_SIZE];
    double d_p_sub[TEST_SIZE];
    subDoubleArrays(double_in, double_in2, d_sub, TEST_SIZE);
    subFloatArrays(float_in, float_in2, f_sub, TEST_SIZE);
	subHalfArrays(half_in, half_in2, h_sub, TEST_SIZE);
    subPositArrays(posit_in, posit_in2, p_sub, TEST_SIZE);
    convertPositArrayToDoubleArray(p_sub, d_p_sub, TEST_SIZE);
    double SNR_SUB_P = computeSNR(d_sub,d_p_sub,TEST_SIZE);
    double SNR_SUB_F = computeSNRFloat(d_sub, f_sub, TEST_SIZE);
	double SNR_SUB_H = computeSNRHalf(d_sub, h_sub, TEST_SIZE);
    std::cout<<"SNR SUB POSIT VS DOUBLE :"<<SNR_SUB_P<<std::endl;
    std::cout<<"SNR SUB FLOAT VS DOUBLE :"<<SNR_SUB_F<<std::endl;
	std::cout<<"SNR SUB HALF VS DOUBLE :"<<SNR_SUB_H<<std::endl;
   /* writeToFile(d_SubFile, d_sub, TEST_SIZE);
    writeToFile(p_SubFile, d_p_sub, TEST_SIZE);
    writeToFileFloat(f_SubFile, f_sub, TEST_SIZE);*/
    //////////////////////////TEST MUL//////////////////////////
    std::cout<<"--------------------------------MUL------------------------"<<std::endl;
	double d_prod[TEST_SIZE];
    float f_prod[TEST_SIZE];
	half h_prod[TEST_SIZE];
    ps_t p_prod[TEST_SIZE];
    double d_p_prod[TEST_SIZE];
    mulDoubleArrays(double_in, double_in2, d_prod, TEST_SIZE);
    mulFloatArrays(float_in, float_in2, f_prod, TEST_SIZE);
	mulHalfArrays(half_in, half_in2, h_prod, TEST_SIZE);
    mulPositArrays(posit_in, posit_in2, p_prod, TEST_SIZE);
    convertPositArrayToDoubleArray(p_prod, d_p_prod, TEST_SIZE);
    double SNR_MUL_P = computeSNR(d_prod,d_p_prod,TEST_SIZE);
    double SNR_MUL_F = computeSNRFloat(d_prod, f_prod, TEST_SIZE);
	double SNR_MUL_H = computeSNRHalf(d_prod, h_prod, TEST_SIZE);
    std::cout<<"SNR MUL POSIT VS DOUBLE :"<<SNR_MUL_P<<std::endl;
    std::cout<<"SNR MUL FLOAT VS DOUBLE :"<<SNR_MUL_F<<std::endl;
	std::cout<<"SNR MUL HALF VS DOUBLE :"<<SNR_MUL_H<<std::endl;
 /*   writeToFile(d_ProdFile, d_prod, TEST_SIZE);
    writeToFile(p_ProdFile, d_p_prod, TEST_SIZE);
    writeToFileFloat(f_ProdFile, f_prod, TEST_SIZE);*/
	//////////////////////////TEST DIV//////////////////////////
	std::cout<<"--------------------------------DIV------------------------"<<std::endl;
    double d_div[TEST_SIZE];
    float f_div[TEST_SIZE];
	half h_div[TEST_SIZE];
    ps_t p_div[TEST_SIZE];
    double d_p_div[TEST_SIZE];
    divDoubleArrays(double_in, double_in2, d_div, TEST_SIZE);
    divFloatArrays(float_in, float_in2, f_div, TEST_SIZE);
	divHalfArrays(half_in, half_in2, h_div, TEST_SIZE);
    divPositArrays(posit_in, posit_in2, p_div, TEST_SIZE);
    convertPositArrayToDoubleArray(p_div, d_p_div, TEST_SIZE);
    double SNR_DIV_P = computeSNR(d_div,d_p_div,TEST_SIZE);
    double SNR_DIV_F = computeSNRFloat(d_div, f_div, TEST_SIZE);
	double SNR_DIV_H = computeSNRHalf(d_div, h_div, TEST_SIZE);
    std::cout<<"SNR DIV POSIT VS DOUBLE :"<<SNR_DIV_P<<std::endl;
    std::cout<<"SNR DIV FLOAT VS DOUBLE :"<<SNR_DIV_F<<std::endl;
	std::cout<<"SNR DIV HALF VS DOUBLE :"<<SNR_DIV_H<<std::endl;
  /*  writeToFile(d_DivFile, d_div, TEST_SIZE);
    writeToFile(p_DivFile, d_p_div, TEST_SIZE);
    writeToFileFloat(f_DivFile, f_div, TEST_SIZE);*/
		//////////////////////////TEST DIV8//////////////////////////
	std::cout<<"--------------------------------DIV8------------------------"<<std::endl;
    double d_div8[TEST_SIZE];
    float f_div8[TEST_SIZE];
	half h_div8[TEST_SIZE];
    ps_t p_div8[TEST_SIZE];
    double d_p_div8[TEST_SIZE];
    div2DoubleArrays(double_in, 8, d_div8, TEST_SIZE);
    div2FloatArrays(float_in, 8, f_div8, TEST_SIZE);
	div2HalfArrays(half_in, 8, h_div8, TEST_SIZE);
    div2PositArrays(posit_in, -3, p_div8, TEST_SIZE);
    convertPositArrayToDoubleArray(p_div8, d_p_div8, TEST_SIZE);
    double SNR_DIV8_P = computeSNR(d_div8,d_p_div8,TEST_SIZE);
    double SNR_DIV8_F = computeSNRFloat(d_div8, f_div8, TEST_SIZE);
	double SNR_DIV8_H = computeSNRHalf(d_div8, h_div8, TEST_SIZE);
    std::cout<<"SNR DIV8 POSIT VS DOUBLE :"<<SNR_DIV8_P<<std::endl;
    std::cout<<"SNR DIV8 FLOAT VS DOUBLE :"<<SNR_DIV8_F<<std::endl;
	std::cout<<"SNR DIV8 HALF VS DOUBLE :"<<SNR_DIV8_H<<std::endl;
		//////////////////////////TEST DIV32//////////////////////////
	std::cout<<"--------------------------------DIV32------------------------"<<std::endl;
    double d_div32[TEST_SIZE];
    float f_div32[TEST_SIZE];
	half h_div32[TEST_SIZE];
    ps_t p_div32[TEST_SIZE];
    double d_p_div32[TEST_SIZE];
    div2DoubleArrays(double_in, 32, d_div32, TEST_SIZE);
    div2FloatArrays(float_in, 32, f_div32, TEST_SIZE);
	div2HalfArrays(half_in, 32, h_div32, TEST_SIZE);
    div2PositArrays(posit_in, -5, p_div32, TEST_SIZE);
    convertPositArrayToDoubleArray(p_div32, d_p_div32, TEST_SIZE);
    double SNR_DIV32_P = computeSNR(d_div32,d_p_div32,TEST_SIZE);
    double SNR_DIV32_F = computeSNRFloat(d_div32, f_div32, TEST_SIZE);
	double SNR_DIV32_H = computeSNRHalf(d_div32, h_div32, TEST_SIZE);
    std::cout<<"SNR DIV32 POSIT VS DOUBLE :"<<SNR_DIV32_P<<std::endl;
    std::cout<<"SNR DIV32 FLOAT VS DOUBLE :"<<SNR_DIV32_F<<std::endl;
	std::cout<<"SNR DIV32 HALF VS DOUBLE :"<<SNR_DIV32_H<<std::endl;
		//////////////////////////TEST DIV128//////////////////////////
	std::cout<<"--------------------------------DIV128------------------------"<<std::endl;
    double d_div128[TEST_SIZE];
    float f_div128[TEST_SIZE];
	half h_div128[TEST_SIZE];
    ps_t p_div128[TEST_SIZE];
    double d_p_div128[TEST_SIZE];
    div2DoubleArrays(double_in, 128, d_div128, TEST_SIZE);
    div2FloatArrays(float_in, 128, f_div128, TEST_SIZE);
	div2HalfArrays(half_in, 128, h_div128, TEST_SIZE);
    div2PositArrays(posit_in, -7, p_div128, TEST_SIZE);
    convertPositArrayToDoubleArray(p_div128, d_p_div128, TEST_SIZE);
    double SNR_DIV128_P = computeSNR(d_div128,d_p_div128,TEST_SIZE);
    double SNR_DIV128_F = computeSNRFloat(d_div128, f_div128, TEST_SIZE);
	double SNR_DIV128_H = computeSNRHalf(d_div128, h_div128, TEST_SIZE);
    
    std::cout<<"SNR DIV128 POSIT VS DOUBLE :"<<SNR_DIV128_P<<std::endl;
    std::cout<<"SNR DIV128 FLOAT VS DOUBLE :"<<SNR_DIV128_F<<std::endl;
	std::cout<<"SNR DIV128 HALF VS DOUBLE :"<<SNR_DIV128_H<<std::endl;
		//////////////////////////TEST DIV1024//////////////////////////
	std::cout<<"--------------------------------DIV1024------------------------"<<std::endl;
    double d_div1024[TEST_SIZE];
    float f_div1024[TEST_SIZE];
	half h_div1024[TEST_SIZE];
    ps_t p_div1024[TEST_SIZE];
    double d_p_div1024[TEST_SIZE];
    div2DoubleArrays(double_in, 1024, d_div1024, TEST_SIZE);
    div2FloatArrays(float_in, 1024, f_div1024, TEST_SIZE);
	div2HalfArrays(half_in, 1024, h_div1024, TEST_SIZE);
    div2PositArrays(posit_in, -10, p_div1024, TEST_SIZE);
    convertPositArrayToDoubleArray(p_div1024, d_p_div1024, TEST_SIZE);
    double SNR_DIV1024_P = computeSNR(d_div1024,d_p_div1024,TEST_SIZE);
    double SNR_DIV1024_F = computeSNRFloat(d_div1024, f_div1024, TEST_SIZE);
	double SNR_DIV1024_H = computeSNRHalf(d_div1024, h_div1024, TEST_SIZE);
    std::cout<<"SNR DIV1024 POSIT VS DOUBLE :"<<SNR_DIV1024_P<<std::endl;
    std::cout<<"SNR DIV1024 FLOAT VS DOUBLE :"<<SNR_DIV1024_F<<std::endl;
	std::cout<<"SNR DIV1024 HALF VS DOUBLE :"<<SNR_DIV1024_H<<std::endl;
		//////////////////////////TEST DIV4096//////////////////////////
	std::cout<<"--------------------------------DIV4096------------------------"<<std::endl;
    double d_div4096[TEST_SIZE];
    float f_div4096[TEST_SIZE];
	half h_div4096[TEST_SIZE];
    ps_t p_div4096[TEST_SIZE];
    double d_p_div4096[TEST_SIZE];
    div2DoubleArrays(double_in, 4096, d_div4096, TEST_SIZE);
    div2FloatArrays(float_in, 4096, f_div4096, TEST_SIZE);
	div2HalfArrays(half_in, 4096, h_div4096, TEST_SIZE);
    div2PositArrays(posit_in, -12, p_div4096, TEST_SIZE);
    convertPositArrayToDoubleArray(p_div4096, d_p_div4096, TEST_SIZE);
    double SNR_DIV4096_P = computeSNR(d_div4096,d_p_div4096,TEST_SIZE);
    double SNR_DIV4096_F = computeSNRFloat(d_div4096, f_div4096, TEST_SIZE);
	double SNR_DIV4096_H = computeSNRHalf(d_div4096, h_div4096, TEST_SIZE);
    std::cout<<"SNR DIV4096 POSIT VS DOUBLE :"<<SNR_DIV4096_P<<std::endl;
    std::cout<<"SNR DIV4096 FLOAT VS DOUBLE :"<<SNR_DIV4096_F<<std::endl;
	std::cout<<"SNR DIV4096 HALF VS DOUBLE :"<<SNR_DIV4096_H<<std::endl;

	//////////////////////////TEST MOD//////////////////////////
    std::cout<<"--------------------------------MOD 2PI------------------------"<<std::endl;
    double d_mod[TEST_SIZE];
    float f_mod[TEST_SIZE];
	half h_mod[TEST_SIZE];
    ps_t p_mod[TEST_SIZE];
    double d_p_mod[TEST_SIZE];
    modDoubleArrays(double_in , d_mod, TEST_SIZE);
    modFloatArrays(float_in, f_mod, TEST_SIZE);
	modHalfArrays(half_in, h_mod, TEST_SIZE);
    modPositArrays(posit_in, p_mod, TEST_SIZE);
    convertPositArrayToDoubleArray(p_mod, d_p_mod, TEST_SIZE);
    double SNR_MOD_P = computeSNR(d_mod,d_p_mod,TEST_SIZE);
    double SNR_MOD_F = computeSNRFloat(d_mod, f_mod, TEST_SIZE);
	double SNR_MOD_H = computeSNRHalf(d_mod, h_mod, TEST_SIZE);
    std::cout<<"SNR MOD POSIT VS DOUBLE :"<<SNR_MOD_P<<std::endl;
    std::cout<<"SNR MOD FLOAT VS DOUBLE :"<<SNR_MOD_F<<std::endl;
	outFile<<"IN;ADD;SUB;MUL;DIV;DIV8;DIV32;DIV128;DIV1024;DIV4096;MOD2PI"<<std::endl;
	outFile<<SNR_IN_P<<";"<<SNR_ADD_P<<";"<<SNR_SUB_P<<";"<<SNR_MUL_P<<";"<<SNR_DIV_P<<";"<<SNR_DIV8_P<<";"<<SNR_DIV32_P<<";"<<SNR_DIV128_P<<";"<<SNR_DIV1024_P<<";"<<SNR_DIV4096_P<<";"<<SNR_MOD_P<<std::endl;
	outFile<<SNR_IN_F<<";"<<SNR_ADD_F<<";"<<SNR_SUB_F<<";"<<SNR_MUL_F<<";"<<SNR_DIV_F<<";"<<SNR_DIV8_F<<";"<<SNR_DIV32_F<<";"<<SNR_DIV128_F<<";"<<SNR_DIV1024_F<<";"<<SNR_DIV4096_F<<";"<<SNR_MOD_F<<std::endl;
	outFile<<SNR_IN_H<<";"<<SNR_ADD_H<<";"<<SNR_SUB_H<<";"<<SNR_MUL_H<<";"<<SNR_DIV_H<<";"<<SNR_DIV8_H<<";"<<SNR_DIV32_H<<";"<<SNR_DIV128_H<<";"<<SNR_DIV1024_H<<";"<<SNR_DIV4096_H<<";"<<SNR_MOD_H<<std::endl;
	outFile.close();

    return 0;

}