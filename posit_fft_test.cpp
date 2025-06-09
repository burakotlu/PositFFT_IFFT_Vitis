#include "posit.hpp"
#include <climits>
#include <etc/hls_sqrt_apfixed.h>
#include <ostream>
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <cstdint>
#include <filesystem>



void loadSignalToArray(const std::string& filename,
                       float floatArray[],
                       half halfArray[],
                       double doubleArray[],
                       ps_t positArray[],
                       int size) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    double sample;
    int i = 0;

    while (file >> sample && i < size) {
        floatArray[i] = static_cast<float>(sample);
        halfArray[i] =(half)sample;
        doubleArray[i] = sample;
        positArray[i] = double2posit(sample);
        ++i;
    }

    file.close();

    std::cout << "Loaded " << i << " values into arrays" << std::endl;
}

void writeToFileFloat(const std::string& filename, const float data_array[], size_t size) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }

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

    for (size_t i = 0; i < size; ++i) {
        outFile <<data_array[i] << "\n";
    }

    outFile.close();
    std::cout << "Saved " << size << " values to " << filename << std::endl;
}

void convertPositArrayToDoubleArray(ps_t positArray[], double doubleArray[], size_t size) {
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
    }
}


double computeSNR(const double original[], const double reconstructed[], size_t size) {
    double signalPower = 0.0, noisePower = 0.0;

    for (size_t i = 0; i < size; ++i) {
        signalPower += original[i] * original[i];
        double error = original[i] - reconstructed[i];
        noisePower += error * error;
    }

    if (noisePower == 0.0) return std::numeric_limits<double>::infinity(); 

    return 10.0 * log10(signalPower / noisePower);
}
double  computeRMSE(const double original[], const double reconstructed[], size_t size) {
    double sumSquaredError = 0.0;

    for (size_t i = 0; i < size; ++i) {
        double error = original[i] - reconstructed[i];
        sumSquaredError += error * error;
    }

    return sqrt(sumSquaredError / size)*1000000;
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
double computeRMSEFloat(const double original[], const float reconstructed[], size_t size) {
    float sumSquaredError = 0.0;

    for (size_t i = 0; i < size; ++i) {
        float error = original[i] - reconstructed[i];
        sumSquaredError += error * error;
    }

    return sqrt(sumSquaredError / size)*1000000;
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
double computeRMSEHalf(const double original[], const half reconstructed[], size_t size) {
    double sumSquaredError = 0.0;

    for (size_t i = 0; i < size; ++i) {
        double error = original[i] - (double)reconstructed[i];
        sumSquaredError += error * error;
    }

    return sqrt(sumSquaredError / size)*1000000;
}
// Converts a double array to float and writes to a file
std::vector<float> convertDoubleArrayToFloat(const std::vector<double>& doubleArray) {
    std::vector<float> floatArray;
    floatArray.reserve(doubleArray.size());

    for (double value : doubleArray) {
        float converted = static_cast<float>(value);
        floatArray.push_back(converted);
    }

    return floatArray;
}
void calculateSinArrays(
    const float floatSinCosIn[], float SinFloat[],
	const half halfSinCosIn[], half SinHalf[],
    const double doubleSinCosIn[], double SinDouble[],
    const ps_t positSinCosIn[], ps_t SinPosit[],
    int size)
{
    for (int i = 0; i < size; ++i) {
        // Calculate sine for float
        SinFloat[i] = fTailorSin(floatSinCosIn[i]);
		SinHalf[i] = hTailorSin(halfSinCosIn[i]);
        // Calculate sine for double
        SinDouble[i] = dTailorSin(doubleSinCosIn[i]);

        // Calculate sine for posit
        SinPosit[i] = positSin(positSinCosIn[i]);  
    }
}
void calculateCosArrays(
    const float floatSinCosIn[], float CosFloat[],
	const half halfSinCosIn[], half CosHalf[],
    const double doubleSinCosIn[], double CosDouble[],
    const ps_t positSinCosIn[], ps_t CosPosit[],
    int size)
{
    for (int i = 0; i < size; ++i) {
        // Calculate sine for float
        CosFloat[i] = fTailorCos(floatSinCosIn[i]);
		CosHalf[i] = hTailorCos(halfSinCosIn[i]);
        // Calculate sine for double
        CosDouble[i] = dTailorCos(doubleSinCosIn[i]);

        // Calculate sine for posit
        CosPosit[i] = positCos(positSinCosIn[i]);  
    }
}
void fillRandomDoublesArray(double* arr, int size, double minVal = -5.0, double maxVal = 5.0) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(minVal, maxVal);

    for (int i = 0; i < size; ++i) {
        arr[i] = dist(gen);
    }
}
int main() {
    std::cout<<"double PI: "<<PI<<std::endl;
   std::cout << "ONE: "  <<posit2double(ONE) << std::endl;
   std::cout << "mONE: "  <<posit2double(mONE) << std::endl;
    std::cout << "POSIT_PI_OVER2: " << posit2double(POSIT_PI_OVER2) << std::endl;
    std::cout << "PI: " << posit2double(POSIT_PI) << std::endl;
    std::cout << "2PI: " << posit2double(POSIT_2PI) << std::endl;
    std::cout << "POSIT_M_PI_OVER2: " << posit2double(POSIT_M_PI_OVER2) << std::endl;
    std::cout << "POSIT_M_PI: " << posit2double(POSIT_M_PI) << std::endl;
    std::cout << "POSIT_M_2PI: " << posit2double(POSIT_M_2PI) << std::endl;

    std::cout << "ZERO+mONE: "  <<posit2double(positAdd(ZERO,mONE)) << std::endl;
    std::cout << "(ZERO+mONE)+mONE: "  <<posit2double(positAdd(positAdd(ZERO,mONE),mONE)) << std::endl;
    
    std::cout<<"------------------------"<<std::endl;
    
    double ang = -278.0101318359375;
    std::cout << "positMod: "  <<posit2double(positMod(double2posit(ang),POSIT_2PI)) << std::endl;
    std::cout << "doubleMod: "  <<fmod(ang,PI*2) << std::endl;
     std::cout<<"------------------------"<<std::endl;
    std::cout << "positCos: "  <<posit2double(positCos(double2posit(ang)) )<< std::endl;
    std::cout << "doubleCos: "  <<dTailorCos(ang) << std::endl;
    std::cout << "floatCos: "  <<fTailorCos(ang) << std::endl;


    double doublekFac[IN_SIZE];
    ps_t positkFac[IN_SIZE];
    double positkFacDouble[IN_SIZE];
    for(int i=0;i<IN_SIZE;i++){
        positkFac[i] = calculateKFactor(i);
        doublekFac[i] = (double) i/IN_SIZE;

    }
    convertPositArrayToDoubleArray(positkFac, positkFacDouble, IN_SIZE);
    double snrPosit_F = computeSNR(doublekFac, positkFacDouble, IN_SIZE);
    double rmsePosit_F = computeRMSE(doublekFac, positkFacDouble, IN_SIZE);


    double dIn[1000];
    ps_t pIn[1000];
    double ptodIn[1000];
    double doubleMd[1000];
    ps_t positMd[1000];
    double positMdDouble[1000];
    fillRandomDoublesArray(dIn,1000,-250,250);
    convertDoubleArrayToPositArray(pIn,dIn,1000);
    convertPositArrayToDoubleArray(pIn, ptodIn, 1000);
    for(int i=0;i<1000;i++){
        positMd[i] = positMod(pIn[i],POSIT_2PI);
        doubleMd[i] = fmod(dIn[i],2*PI);

    }
    convertPositArrayToDoubleArray(positMd, positMdDouble, 1000);
    double snrPosit_M_IN = computeSNR(dIn, ptodIn, IN_SIZE);
    double snrPosit_M = computeSNR(doubleMd, positMdDouble, IN_SIZE);
    double rmsePosit_M_IN = computeRMSE(dIn, ptodIn, IN_SIZE);
    double rmsePosit_M = computeRMSE(doubleMd, positMdDouble, IN_SIZE);

    std::string appr_suffix = (APPR_TAILOR == 1) ? "_APP" : "_NAPP";
    
    // Construct the base path
    std::string basePath = "C:/Users/Burak/Desktop/PFFT/PFFT/TERMS_" 
                            + std::to_string(TERMS) + "_N_" + std::to_string(N) 
                            + "_ES_" + std::to_string(ES)
                            + "_SIZE_" + std::to_string(IN_SIZE) + appr_suffix;
    std::string inputFile = "C:/Users/Burak/Desktop/PFFT/PFFT/random_numbers_" + std::to_string(IN_SIZE) + ".txt";  // Input file with the time-domain signal
    std::string inputSinCosFile = "C:/Users/Burak/Desktop/PFFT/PFFT/angle_values_" + std::to_string(IN_SIZE) + ".txt";  // Input file with the time-domain signal
	std::string signalFile = basePath + "/signal.txt";  // Input file with the time-domain signal
    std::string dOutputFileReal = basePath + "/double_fft_real.txt";
    std::string dOutputFileImag = basePath + "/double_fft_imag.txt";
    std::string fOutputFileReal = basePath + "/float_fft_real.txt";
    std::string fOutputFileImag = basePath + "/float_fft_imag.txt";
    std::string hOutputFileReal = basePath + "/half_fft_real.txt";
    std::string hOutputFileImag = basePath + "/half_fft_imag.txt";
    std::string pOutputFileReal = basePath + "/posit_fft_real.txt";
    std::string pOutputFileImag = basePath + "/posit_fft_imag.txt";
    std::string snrFFTOut = basePath + "/snr_fft.txt";
	std::string rmseFFTOut = basePath + "/rmse_fft.txt";
    std::string difftOutputFile = basePath + "/double_ifft_output.txt";
    std::string fifftOutputFile = basePath + "/float_ifft_output.txt";
	std::string hifftOutputFile = basePath + "/half_ifft_output.txt";
    std::string pifftOutputFile = basePath + "/posit_ifft_output.txt";
    std::string snrIFFTOut = basePath + "/snr_ifft.txt";
	std::string rmseIFFTOut = basePath + "/rmse_ifft.txt";
    std::string inDouble = basePath + "/in_double.txt";
    std::string inFloat = basePath + "/in_float.txt";
    std::string inHalf = basePath + "/in_half.txt";
    std::string inPosit = basePath + "/in_posit.txt";
	std::string inSinCosDouble = basePath + "/inSinCos_double.txt";
    std::string inSinCosFloat = basePath + "/inSinCos_float.txt";
	std::string inSinCosHalf = basePath + "/inSinCos_half.txt";
    std::string inSinCosPosit = basePath + "/inSinCos_posit.txt";
	std::string sinDouble = basePath + "/sin_double.txt";
    std::string sinFloat = basePath + "/sin_float.txt";
	std::string sinHalf = basePath + "/sin_half.txt";
    std::string sinPosit = basePath + "/sin_posit.txt";
	std::string cosDouble = basePath + "/cos_double.txt";
    std::string cosFloat = basePath + "/cos_float.txt";
	std::string cosHalf = basePath + "/cos_half.txt";
    std::string cosPosit = basePath + "/cos_posit.txt";
	std::string snrSINOut = basePath + "/snr_sin.txt";
	std::string snrCOSOut = basePath + "/snr_cos.txt";
	std::string rmseSINOut = basePath + "/rmse_sin.txt";
	std::string rmseCOSOut = basePath + "/rmse_cos.txt";
    std::cout << "Before Load" << std::endl;

    // Step 1: Load the signal from the txt file into a HLS stream
    float floatSignal[IN_SIZE];
    half halfSignal[IN_SIZE];
    double doubleSignal[IN_SIZE];
    ps_t positSignal[IN_SIZE];
    double positSignalDouble[IN_SIZE];
    loadSignalToArray(inputFile, floatSignal,halfSignal, doubleSignal, positSignal, IN_SIZE);

    convertPositArrayToDoubleArray(positSignal, positSignalDouble, IN_SIZE);
    writeToFile(inDouble, doubleSignal, IN_SIZE);
    writeToFileFloat(inFloat, floatSignal, IN_SIZE);
	writeToFileHalf(inHalf,halfSignal,IN_SIZE);
    writeToFile(inPosit, positSignalDouble, IN_SIZE);
	
    double snrPosit_IN = computeSNR(doubleSignal, positSignalDouble, IN_SIZE);
    double snrFloat_IN = computeSNRFloat(doubleSignal,floatSignal, IN_SIZE);
	double snrHalf_IN = computeSNRHalf(doubleSignal,halfSignal, IN_SIZE);

    double rmsePosit_IN = computeRMSE(doubleSignal, positSignalDouble, IN_SIZE);
    double rmseFloat_IN = computeRMSEFloat(doubleSignal,floatSignal, IN_SIZE);
	double rmseHalf_IN = computeRMSEHalf(doubleSignal,halfSignal, IN_SIZE);
	
    std::cout << "After Load" << std::endl;
    float floatSinCosIn[IN_SIZE];
	half halfSinCosIn[IN_SIZE];
    double doubleSinCosIn[IN_SIZE];
    ps_t positSinCosIn[IN_SIZE];
    double positSinCosInDouble[IN_SIZE];
    loadSignalToArray(inputSinCosFile, floatSinCosIn,halfSinCosIn, doubleSinCosIn, positSinCosIn, IN_SIZE);
	convertPositArrayToDoubleArray(positSinCosIn, positSinCosInDouble, IN_SIZE);
	writeToFile(inSinCosDouble, doubleSinCosIn, IN_SIZE);
    writeToFileFloat(inSinCosFloat, floatSinCosIn, IN_SIZE);
	writeToFileHalf(inSinCosHalf, halfSinCosIn, IN_SIZE);
    writeToFile(inSinCosPosit, positSinCosInDouble, IN_SIZE);
	std::cout << "----------SIN--------------" << std::endl;
    float SinFloat[IN_SIZE];
	half SinHalf[IN_SIZE];
    double SinDouble[IN_SIZE];
    ps_t SinPosit[IN_SIZE];
	double SinPositDouble[IN_SIZE];
    calculateSinArrays(floatSinCosIn,SinFloat,halfSinCosIn, SinHalf,doubleSinCosIn,SinDouble,positSinCosIn,SinPosit,IN_SIZE);
	convertPositArrayToDoubleArray(SinPosit, SinPositDouble, IN_SIZE);
	writeToFile(sinDouble, SinDouble, IN_SIZE);        
    writeToFile(sinPosit, SinPositDouble, IN_SIZE);       
    writeToFileFloat(sinFloat, SinFloat, IN_SIZE);  
	writeToFileHalf(sinHalf, SinHalf, IN_SIZE); 
	
    double snrPosit_SIN = computeSNR(SinDouble, SinPositDouble, IN_SIZE);
    double snrFloat_SIN = computeSNRFloat(SinDouble,SinFloat, IN_SIZE); 
	double snrHalf_SIN = computeSNRHalf(SinDouble,SinHalf, IN_SIZE);   
    std::ofstream outFileSIN(snrSINOut);
    outFileSIN << "Posit: " << snrPosit_SIN << std::endl;
    outFileSIN << "Float: " << snrFloat_SIN << std::endl;
	outFileSIN << "Half: " << snrHalf_SIN << std::endl;
    outFileSIN.close();
	double rmsePosit_SIN = computeRMSE(SinDouble, SinPositDouble, IN_SIZE);
    double rmseFloat_SIN = computeRMSEFloat(SinDouble,SinFloat, IN_SIZE); 
	double rmseHalf_SIN = computeRMSEHalf(SinDouble,SinHalf, IN_SIZE);   
    std::ofstream outFileSinRMSE(rmseSINOut);
    outFileSinRMSE << "Posit: " << rmsePosit_SIN << std::endl;
    outFileSinRMSE << "Float: " << rmseFloat_SIN << std::endl;
	outFileSinRMSE << "Half: " << rmseHalf_SIN << std::endl;
    outFileSinRMSE.close();
	std::cout << "----------COS--------------" << std::endl;
    float CosFloat[IN_SIZE];
	half CosHalf[IN_SIZE];
    double CosDouble[IN_SIZE];
    ps_t CosPosit[IN_SIZE];
	double CosPositDouble[IN_SIZE];
    calculateCosArrays(floatSinCosIn,CosFloat,halfSinCosIn,CosHalf,doubleSinCosIn,CosDouble,positSinCosIn,CosPosit,IN_SIZE);
	convertPositArrayToDoubleArray(CosPosit, CosPositDouble, IN_SIZE);
	writeToFile(cosDouble, CosDouble, IN_SIZE);        
    writeToFile(cosPosit, CosPositDouble, IN_SIZE);       
    writeToFileFloat(cosFloat, CosFloat, IN_SIZE);  
	writeToFileHalf(cosHalf, CosHalf, IN_SIZE);
	
	double snrPosit_COS = computeSNR(CosDouble, CosPositDouble, IN_SIZE);
    double snrFloat_COS = computeSNRFloat(CosDouble,CosFloat, IN_SIZE);  
	double snrHalf_COS = computeSNRHalf(CosDouble,CosHalf, IN_SIZE); 	
    std::ofstream outFileCOS(snrCOSOut);
    outFileCOS << "Posit: " << snrPosit_COS << std::endl;
    outFileCOS << "Float: " << snrFloat_COS << std::endl;
	outFileCOS << "Half: " << snrHalf_COS << std::endl;
    outFileCOS.close();
	
	double rmsePosit_COS = computeRMSE(CosDouble, CosPositDouble, IN_SIZE);
    double rmseFloat_COS = computeRMSEFloat(CosDouble,CosFloat, IN_SIZE);  
	double rmseHalf_COS = computeRMSEHalf(CosDouble,CosHalf, IN_SIZE); 	
    std::ofstream outFileCosRMSE(rmseCOSOut);
    outFileCosRMSE << "Posit: " << rmsePosit_COS << std::endl;
    outFileCosRMSE << "Float: " << rmseFloat_COS << std::endl;
	outFileCosRMSE << "Half: " << rmseHalf_COS << std::endl;
    outFileCosRMSE.close();	
    std::cout << "----------FFT--------------" << std::endl;
    std::cout << "Calculating double results..." << std::endl;
    dFFTResult d_fftResult;
    dFFT(doubleSignal, d_fftResult);

    std::cout << "Calculating float results..." << std::endl;
    fFFTResult f_fftResult;
    fFFT(floatSignal, f_fftResult);
	
    std::cout << "Calculating half results..." << std::endl;
    hFFTResult h_fftResult;
    hFFT(halfSignal, h_fftResult);
    std::cout << "Calculating posit results..." << std::endl;
    pFFTResult pfftResult;
    pFFT(positSignal, pfftResult);

    double RC_real[IN_SIZE];  
    double RC_imag[IN_SIZE];

    convertPositArrayToDoubleArray(pfftResult.real, RC_real, IN_SIZE);
    convertPositArrayToDoubleArray(pfftResult.imag, RC_imag, IN_SIZE);

    std::cout << "Writing results in files..." << std::endl;
// Assuming data is stored in arrays like RC_real, RC_imag, d_fftResult.real, d_fftResult.imag
    writeToFile(signalFile, doubleSignal, IN_SIZE);        // For signalArray
    writeToFile(pOutputFileReal, RC_real, IN_SIZE);       // For real part of posit FFT result
    writeToFile(pOutputFileImag, RC_imag, IN_SIZE);      // For imaginary part of posit FFT result
    writeToFile(dOutputFileReal, d_fftResult.real, IN_SIZE);  // For real part of double FFT result
    writeToFile(dOutputFileImag, d_fftResult.imag, IN_SIZE);  // For imaginary part of double FFT result

    writeToFileFloat(fOutputFileReal, f_fftResult.real,IN_SIZE);
    writeToFileFloat(fOutputFileImag, f_fftResult.imag,IN_SIZE);
    writeToFileHalf(hOutputFileReal, h_fftResult.real,IN_SIZE);
    writeToFileHalf(hOutputFileImag, h_fftResult.imag,IN_SIZE);
    // Compute SNRs

    double snrPosit_FFT = computeSNR(d_fftResult.real, RC_real, IN_SIZE);
    double snrFloat_FFT = computeSNRFloat(d_fftResult.real,f_fftResult.real, IN_SIZE);
    double snrHalf_FFT = computeSNRHalf(d_fftResult.real,h_fftResult.real, IN_SIZE);
   
	std::ofstream outFileFFT(snrFFTOut);
    outFileFFT << "Posit: " << snrPosit_FFT << std::endl;
    outFileFFT << "Float: " << snrFloat_FFT << std::endl;
	outFileFFT << "Half: " << snrHalf_FFT << std::endl;
    outFileFFT.close();
    // Compute RMSE

    double rmsePosit_FFT = computeRMSE(d_fftResult.real, RC_real, IN_SIZE);
    double rmseFloat_FFT = computeRMSEFloat(d_fftResult.real,f_fftResult.real, IN_SIZE);
    double rmseHalf_FFT = computeRMSEHalf(d_fftResult.real,h_fftResult.real, IN_SIZE);
   
	std::ofstream outFileFFTRMSE(rmseFFTOut);
    outFileFFTRMSE << "Posit: " << rmsePosit_FFT << std::endl;
    outFileFFTRMSE << "Float: " << rmseFloat_FFT << std::endl;
	outFileFFTRMSE << "Half: " << rmseHalf_FFT << std::endl;
    outFileFFTRMSE.close();

    std::cout << "Computing double IFFT..." << std::endl;
    double d_reconstructedSignal[IN_SIZE];
    double d_real[IN_SIZE], d_imag[IN_SIZE];
    for(int i=0;i<IN_SIZE;i++){
        d_real[i] = d_fftResult.real[i];
        d_imag[i] = -d_fftResult.imag[i];
        d_fftResult.imag[i] = d_imag[i];
        //d_real[i] = 1.0f;
        //d_imag[i] = 1.0f;
    }
    dIFFT(d_fftResult, d_reconstructedSignal, IN_SIZE);
    std::cout << "Writing IFFT results to file..." << std::endl;
    writeToFile(difftOutputFile, d_reconstructedSignal, IN_SIZE);
    // Conversion buffers
    float f_real[IN_SIZE], f_imag[IN_SIZE];
	half h_real[IN_SIZE], h_imag[IN_SIZE];
    ps_t p_real[IN_SIZE], p_imag[IN_SIZE];
	fFFTResult f_fftResult2;
	hFFTResult h_fftResult2;
	pFFTResult p_fftResult2;
    // Convert double FFT result to float and posit
    for (int i = 0; i < IN_SIZE; i++) {

        f_fftResult2.real[i] = static_cast<float>(d_real[i]);
        f_fftResult2.imag[i] = static_cast<float>(d_imag[i]);

        h_fftResult2.real[i] = static_cast<half>(d_real[i]);
        h_fftResult2.imag[i] = static_cast<half>(d_imag[i]);
		
        p_fftResult2.real[i] = double2posit(d_real[i]);
        p_fftResult2.imag[i] = double2posit(d_imag[i]);
    }
    std::cout << "Computing float IFFT..." << std::endl;
    float f_reconstructedSignal[IN_SIZE];
    fIFFT(f_fftResult2, f_reconstructedSignal, IN_SIZE);
    writeToFileFloat(fifftOutputFile, f_reconstructedSignal, IN_SIZE);

    half h_reconstructedSignal[IN_SIZE];
    hIFFT(h_fftResult2, h_reconstructedSignal, IN_SIZE);
    writeToFileHalf(hifftOutputFile, h_reconstructedSignal, IN_SIZE);
	
    std::cout << "Computing posit IFFT..." << std::endl;
    ps_t p_reconstructedSignal[IN_SIZE];

    pIFFT(p_fftResult2, p_reconstructedSignal, IN_SIZE);

    double RC_result[IN_SIZE];
    convertPositArrayToDoubleArray(p_reconstructedSignal, RC_result, IN_SIZE);
    writeToFile(pifftOutputFile, RC_result, IN_SIZE);
     // Compute SNRs for IFFT
    double snrPosit_IFFT = computeSNR(d_reconstructedSignal, RC_result, IN_SIZE);
    double snrFloat_IFFT = computeSNRFloat(d_reconstructedSignal, f_reconstructedSignal, IN_SIZE);
    double snrHalf_IFFT = computeSNRHalf(d_reconstructedSignal, h_reconstructedSignal, IN_SIZE);
   
	std::ofstream outFileIFFT(snrIFFTOut);

    outFileIFFT << "-------------WRT ORIG-------------------------" << std::endl;
    outFileIFFT << "Posit: " << snrPosit_IFFT << std::endl;
    outFileIFFT << "Float: " << snrFloat_IFFT << std::endl;
	outFileIFFT << "Half: " << snrHalf_IFFT << std::endl;
    outFileIFFT.close();
	
     // Compute RMSE for IFFT
    double rmsePosit_IFFT = computeRMSE(d_reconstructedSignal, RC_result, IN_SIZE);
    double rmseFloat_IFFT = computeRMSEFloat(d_reconstructedSignal, f_reconstructedSignal, IN_SIZE);
    double rmseHalf_IFFT = computeRMSEHalf(d_reconstructedSignal, h_reconstructedSignal, IN_SIZE);
   
	std::ofstream outFileIFFT_RMSE(rmseIFFTOut);

    outFileIFFT_RMSE << "-------------WRT ORIG-------------------------" << std::endl;
    outFileIFFT_RMSE << "Posit: " << rmsePosit_IFFT << std::endl;
    outFileIFFT_RMSE << "Float: " << rmseFloat_IFFT << std::endl;
	outFileIFFT_RMSE << "Half: " << rmseHalf_IFFT << std::endl;
    outFileIFFT_RMSE.close();
	
    // Compute SNRs with respect to original signal
    double snrPosit_IFFT_ORIG = computeSNR(doubleSignal, RC_result,IN_SIZE);
    double snrFloat_IFFT_ORIG = computeSNRFloat(doubleSignal, f_reconstructedSignal,IN_SIZE);
    double snrHalf_IFFT_ORIG = computeSNRHalf(doubleSignal, h_reconstructedSignal,IN_SIZE);
	double snrDouble_IFFT_ORIG = computeSNR(doubleSignal, d_reconstructedSignal,IN_SIZE);
    // Compute RMSE with respect to original signal
    double rmsePosit_IFFT_ORIG = computeRMSE(doubleSignal, RC_result,IN_SIZE);
    double rmseFloat_IFFT_ORIG = computeRMSEFloat(doubleSignal, f_reconstructedSignal,IN_SIZE);
    double rmseHalf_IFFT_ORIG = computeRMSEHalf(doubleSignal, h_reconstructedSignal,IN_SIZE);
	double rmseDouble_IFFT_ORIG = computeRMSE(doubleSignal, d_reconstructedSignal,IN_SIZE);
    std::cout << "----------------SNR RESULTS K Factor-------------------------" << std::endl;
    std::cout << "SNR k factor Calculation: " << snrPosit_F << " dB" << std::endl;
    std::cout << "----------------SNR RESULTS MOD-------------------------" << std::endl;
    std::cout << "SNR Before Calculation: " << snrPosit_M_IN << " dB" << std::endl;
    std::cout << "SNR Mod Calculation: " << snrPosit_M << " dB" << std::endl;
    std::cout << "----------------SNR RESULTS INPUT-------------------------" << std::endl;
    std::cout << "SNR (Posit Reconstruction): " << snrPosit_IN << " dB" << std::endl;
    std::cout << "SNR (Float Conversion): " << snrFloat_IN << " dB" << std::endl;
    std::cout << "SNR (Half Conversion): " << snrHalf_IN << " dB" << std::endl;    
	std::cout << "----------------SNR RESULTS SIN-------------------------" << std::endl;
    std::cout << "SNR (Posit Reconstruction): " << snrPosit_SIN << " dB" << std::endl;
    std::cout << "SNR (Float Conversion): " << snrFloat_SIN << " dB" << std::endl;
    std::cout << "SNR (Half Conversion): " << snrHalf_SIN << " dB" << std::endl;   
    std::cout << "----------------SNR RESULTS COS-------------------------" << std::endl;
    std::cout << "SNR (Posit Reconstruction): " << snrPosit_COS << " dB" << std::endl;
    std::cout << "SNR (Float Conversion): " << snrFloat_COS << " dB" << std::endl;
    std::cout << "SNR (Half Conversion): " << snrHalf_COS << " dB" << std::endl;    
	std::cout << "----------------SNR RESULTS FFT-------------------------" << std::endl;
    std::cout << "SNR (Posit Reconstruction): " << snrPosit_FFT << " dB" << std::endl;
    std::cout << "SNR (Float Conversion): " << snrFloat_FFT << " dB" << std::endl;
    std::cout << "SNR (Half Conversion): " << snrHalf_FFT << " dB" << std::endl;    
	std::cout << "----------------SNR RESULTS IFFT-------------------------" << std::endl;
    std::cout << "SNR (Posit Reconstruction): " << snrPosit_IFFT << " dB" << std::endl;
    std::cout << "SNR (Float Conversion): " << snrFloat_IFFT << " dB" << std::endl;
    std::cout << "SNR (Half Conversion): " << snrHalf_IFFT << " dB" << std::endl;    
	std::cout << "----------------SNR RESULTS IFFT WRT ORIG-----------------" << std::endl;
    std::cout << "SNR (Posit Reconstruction): " << snrPosit_IFFT_ORIG << " dB" << std::endl;
    std::cout << "SNR (Float Conversion): " << snrFloat_IFFT_ORIG << " dB" << std::endl;
    std::cout << "SNR (Half Conversion): " << snrHalf_IFFT_ORIG << " dB" << std::endl;    
	std::cout << "SNR (Double Conversion): " << snrDouble_IFFT_ORIG << " dB" << std::endl;
    std::cout << "----------------SNR RESULTS MOD-------------------------" << std::endl;
    std::cout << "SNR Before Calculation: " << snrPosit_M_IN << " dB" << std::endl;
    std::cout << "SNR Mod Calculation: " << snrPosit_M << " dB" << std::endl;

    std::cout << "----------------RMSE RESULTS K Factor-------------------------" << std::endl;
    std::cout << "RMSE k factor Calculation: " << rmsePosit_F << std::endl;
    std::cout << "----------------RMSE RESULTS MOD-------------------------" << std::endl;
    std::cout << "RMSE Before Calculation: " << rmsePosit_M_IN << std::endl;
    std::cout << "RMSE Mod Calculation: " << rmsePosit_M << std::endl;
    std::cout << "----------------RMSE RESULTS INPUT-------------------------" << std::endl;
    std::cout << "RMSE (Posit Reconstruction): " << rmsePosit_IN << std::endl;
    std::cout << "RMSE (Float Conversion): " << rmseFloat_IN << std::endl;
    std::cout << "RMSE (Half Conversion): " << rmseHalf_IN << std::endl;    
	std::cout << "----------------RMSE RESULTS SIN-------------------------" << std::endl;
    std::cout << "RMSE (Posit Reconstruction): " << rmsePosit_SIN << std::endl;
    std::cout << "RMSE (Float Conversion): " << rmseFloat_SIN << std::endl;
    std::cout << "RMSE (Half Conversion): " << rmseHalf_SIN << std::endl;   
    std::cout << "----------------RMSE RESULTS COS-------------------------" << std::endl;
    std::cout << "RMSE (Posit Reconstruction): " << rmsePosit_COS << std::endl;
    std::cout << "RMSE (Float Conversion): " << rmseFloat_COS << std::endl;
    std::cout << "RMSE (Half Conversion): " << rmseHalf_COS << std::endl;    
	std::cout << "----------------RMSE RESULTS FFT-------------------------" << std::endl;
    std::cout << "RMSE (Posit Reconstruction): " << rmsePosit_FFT << std::endl;
    std::cout << "RMSE (Float Conversion): " << rmseFloat_FFT << std::endl;
    std::cout << "RMSE (Half Conversion): " << rmseHalf_FFT << std::endl;    
	std::cout << "----------------RMSE RESULTS IFFT-------------------------" << std::endl;
    std::cout << "RMSE (Posit Reconstruction): " << rmsePosit_IFFT << std::endl;
    std::cout << "RMSE (Float Conversion): " << rmseFloat_IFFT << std::endl;
    std::cout << "RMSE (Half Conversion): " << rmseHalf_IFFT << std::endl;    
	std::cout << "----------------RMSE RESULTS IFFT WRT ORIG-----------------" << std::endl;
    std::cout << "RMSE (Posit Reconstruction): " << rmsePosit_IFFT_ORIG << std::endl;
    std::cout << "RMSE (Float Conversion): " << rmseFloat_IFFT_ORIG << std::endl;
    std::cout << "RMSE (Half Conversion): " << rmseHalf_IFFT_ORIG << std::endl;    
	std::cout << "RMSE (Double Conversion): " << rmseDouble_IFFT_ORIG << std::endl;
    std::cout << "----------------RMSE RESULTS MOD-------------------------" << std::endl;
    std::cout << "RMSE Before Calculation: " << rmsePosit_M_IN << std::endl;
    std::cout << "RMSE Mod Calculation: " << rmsePosit_M << std::endl;


    std::cout << "Processing complete." << std::endl;

    return 0;
}