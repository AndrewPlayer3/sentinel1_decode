#include "signal_processing.hpp"

using namespace std;


vector<float> flatten(const vector<vector<float>>& values)
{
    int rows = values.size();
    int cols = values[0].size();

    vector<float> flat(rows * cols);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            flat[i * cols + j] = values[i][j];
        }
    }

    return flat;
}


vector<vector<float>> norm_2d(const vector<vector<complex<float>>>& complex_values)
{
    float max_value = 0;    

    int rows = complex_values.size();
    int cols = complex_values[0].size();

    vector<vector<float>> norm(rows, vector<float>(cols));

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            float real = complex_values[i][j].real();
            float imag = complex_values[i][j].imag();

            float norm_value = sqrt(pow(real, 2) + pow(imag, 2));

            norm[i][j] = norm_value;

            if (norm_value > max_value) max_value = norm_value;
        }
    }

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            norm[i][j] = norm[i][j] / max_value;
        }
    }

    return norm;
}


vector<float> norm_1d(const vector<complex<float>>& complex_values)
{
    int num_samples = complex_values.size();

    float max_value = 0;

    vector<float> norm(num_samples);

    for (int i = 0; i < num_samples; i++)
    {
        float real = complex_values[i].real();
        float imag = complex_values[i].imag();

        float norm_value = sqrt(pow(real, 2) + pow(imag, 2));

        norm[i] = norm_value;

        if (norm_value > max_value) max_value = norm_value;
    }

    for (int i = 0; i < num_samples; i++)
    {
        norm[i] = norm[i] / max_value;
    }

    return norm;
}


vector<complex<float>> compute_1d_dft(const vector<complex<float>>& signal)
{
    cout << "Initializing 1D Complex Vector for FFTW" << endl;
    
    int fft_size = signal.size();

    vector<complex<float>> signal_fft = signal;

    cout << "Executing 1D DFT Plan" << endl;

    fftwf_plan plan;

    plan = fftwf_plan_dft_1d(
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


vector<vector<complex<float>>> compute_2d_dft(const vector<vector<complex<float>>>& signal)
{
    cout << "Initializing 1D Complex Vector for FFTW" << endl;

    int fft_rows = signal.size();
    int fft_cols = signal[0].size();

    vector<complex<float>> signal_fftw_io(fft_rows * fft_cols);

    for (int i = 0; i < fft_rows; i++)
    {
        for (int j = 0; j < fft_cols; j++)
        {
            signal_fftw_io[i*fft_cols+j] = signal[i][j];
        }
    }

    cout << "Executing 2D DFT Plan" << endl;

    fftwf_plan plan = fftwf_plan_dft_2d(
        fft_rows, fft_cols, 
        reinterpret_cast<fftwf_complex*>(signal_fftw_io.data()),
        reinterpret_cast<fftwf_complex*>(signal_fftw_io.data()),
        FFTW_FORWARD, FFTW_ESTIMATE
    );

    fftwf_execute(plan);

    cout << "Copying Complex Data into 2D Vector" << endl;

    vector<vector<complex<float>>> signal_fft(fft_rows, vector<complex<float>>(fft_cols));

    for (int i = 0; i < fft_rows; i++)
    {
        for (int j = 0; j < fft_cols; j++)
        {
            signal_fft[i][j] = signal_fftw_io[i*fft_cols+j];
        }
    }

    cout << "Destroying DFT Plan" << endl;

    fftwf_destroy_plan(plan);

    return signal_fft;
}
