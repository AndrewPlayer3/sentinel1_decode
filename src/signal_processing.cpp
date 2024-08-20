#include "signal_processing.hpp"

using namespace std;


vector<vector<float>> norm_2d(const vector<vector<complex<float>>>& complex_values, const bool& log_scale = false)
{
    float max_value = 0.0;    

    int rows = complex_values.size();
    int cols = complex_values[0].size();

    vector<vector<float>> norm(rows, vector<float>(cols));

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            float real = complex_values[i][j].real();
            float imag = complex_values[i][j].imag();

            float mag = sqrt((real * real) + (imag * imag));

            norm[i][j] = mag;

            if (mag > max_value) max_value = mag;
        }
    }

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            norm[i][j] = log_scale ? 20 * log10(norm[i][j] / max_value) : norm[i][j] / max_value;
        }
    }

    return norm;
}


vector<float> norm_1d(const vector<complex<float>>& complex_values, const bool& log_scale = false)
{
    int num_samples = complex_values.size();

    float max_value = 0;

    vector<float> norm(num_samples);

    for (int i = 0; i < num_samples; i++)
    {
        float real = complex_values[i].real();
        float imag = complex_values[i].imag();

        float mag = sqrt(pow(real, 2) + pow(imag, 2));

        norm[i] = mag;

        if (mag > max_value) max_value = mag;
    }

    for (int i = 0; i < num_samples; i++)
    {
        norm[i] = log_scale ? 20 * log10(norm[i] / max_value) : norm[i] / max_value;
    }

    return norm;
}


vector<float> magnitude_1d(const vector<complex<float>>& complex_values)
{
    int num_samples = complex_values.size();

    vector<float> magnitude(num_samples);

    #pragma omp parallel for
    for (int i = 0; i < num_samples; i++)
    {
        float real = complex_values[i].real();
        float imag = complex_values[i].imag();

        float mag = sqrt(pow(real, 2) + pow(imag, 2));

        magnitude[i] = mag;    
    }

    return magnitude;
}


vector<vector<float>> magnitude_2d(const vector<vector<complex<float>>>& complex_values)
{
    float max_value = 0;    

    int rows = complex_values.size();
    int cols = complex_values[0].size();

    vector<vector<float>> magnitude(rows, vector<float>(cols));

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            float real = complex_values[i][j].real();
            float imag = complex_values[i][j].imag();

            float norm_value = sqrt(pow(real, 2) + pow(imag, 2));

            magnitude[i][j] = norm_value;
        }
    }

    return magnitude;
}


vector<complex<float>> compute_1d_dft(const vector<complex<float>>& signal)
{
    cout << "Initializing 1D Complex Vector for FFTW" << endl;
    
    int fft_threads = 8;
    int fft_size    = signal.size();

    vector<complex<float>> signal_fft = signal;

    cout << "Executing 1D DFT Plan" << endl;

    fftwf_plan plan = fftwf_plan_dft_1d(
        fft_size,
        reinterpret_cast<fftwf_complex*>(signal_fft.data()),
        reinterpret_cast<fftwf_complex*>(signal_fft.data()),
        FFTW_FORWARD,
        FFTW_ESTIMATE
    );

    fftwf_execute(plan);

    fftwf_destroy_plan(plan);

    return signal_fft;
}


vector<vector<complex<float>>> compute_2d_dft(
    const vector<vector<complex<float>>>& signal,
    int fft_rows = 0,
    int fft_cols = 0,
    const bool& inverse = false
) {
    cout << "Initializing 1D Complex Vector for FFTW" << endl;

    if (not fft_rows) fft_rows = signal.size();
    if (not fft_cols) fft_cols = signal[0].size();

    vector<complex<float>> signal_fftw_io(fft_rows * fft_cols);

    for (int i = 0; i < fft_rows; i++)
    {
        for (int j = 0; j < fft_cols; j++)
        {
            signal_fftw_io[i*fft_cols+j] = signal[i][j];
        }
    }

    cout << "Executing 2D DFT Plan" << endl;

    int fft_threads   = 8;
    int fft_direction = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    fftwf_plan plan = fftwf_plan_dft_2d(
        fft_rows, fft_cols, 
        reinterpret_cast<fftwf_complex*>(signal_fftw_io.data()),
        reinterpret_cast<fftwf_complex*>(signal_fftw_io.data()),
        FFTW_FORWARD, FFTW_ESTIMATE
    );

    fftwf_execute(plan);

    cout << "Destroying DFT Plan" << endl;

    fftwf_destroy_plan(plan);

    cout << "Copying Complex Data into 2D Vector" << endl;

    vector<vector<complex<float>>> signal_fft(fft_rows, vector<complex<float>>(fft_cols));

    for (int i = 0; i < fft_rows; i++)
    {
        for (int j = 0; j < fft_cols; j++)
        {
            signal_fft[i][j] = signal_fftw_io[i*fft_cols+j];
        }
    }

    return signal_fft;
}
