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


vector<complex<float>> compute_1d_dft(const vector<complex<float>>& signal, int fft_size = 0, const bool& inverse = false)
{
    cout << "Initializing 1D Complex Vector for FFTW" << endl;

    if (fft_size <= 0) fft_size = signal.size();

    vector<complex<float>> fft_vector = signal;

    int fft_direction = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    cout << "Executing 1D DFT Plan" << endl;

    fftwf_plan plan = fftwf_plan_dft_1d(
        fft_size,
        reinterpret_cast<fftwf_complex*>(fft_vector.data()),
        reinterpret_cast<fftwf_complex*>(fft_vector.data()),
        fft_direction,
        FFTW_ESTIMATE
    );

    fftwf_execute(plan);

    cout << "Destroying DFT Plan" << endl;

    fftwf_destroy_plan(plan);

    if (inverse)
    {
        for (int i = 0; i < fft_size; i++)
        {
            fft_vector[i] /= fft_size;
        }
    }

    return fft_vector;
}


// TODO: Find matrix math library to do this with efficiency.
vector<vector<complex<float>>> compute_1d_dft(const vector<vector<complex<float>>>& signals, int fft_size = 0, const int& axis = 0, const bool& inverse = false)
{
    if (axis != 0 and axis != 1) 
    {
        throw invalid_argument("FFT axis must be 0 (rows) or 1 (cols).");
    }

    int signal_rows   = signals.size();
    int signal_cols   = signals[0].size();
    int min_fft_size  = 8;
    int max_fft_size  = axis ? signal_cols : signal_rows;
    int fft_direction = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    if (fft_size == 0)
    {
        fft_size = max_fft_size;
    }
    else if (fft_size < min_fft_size or fft_size > max_fft_size)
    {
        throw invalid_argument("FFT size must be in [8, signal_length_for_axis].");
    }

    cout << "Initializing Complex Arrays and FFTW Plan" << endl;
    
    // axis = 0 => fft_array_in  is (signal_cols, signal_rows)
    //             fft_array_out is (signal_cols, fft_size)
    //
    // axis = 1 => fft_array_in  is (signal_rows, signal_cols)
    //             fft_array_out is (signal_rows, fft_size)

    int out_array_rows = axis ? signal_rows : signal_cols;

    vector<vector<complex<float>>> fft_array_in;
    vector<vector<complex<float>>> fft_array_out(out_array_rows, vector<complex<float>>(fft_size));
    vector<fftwf_plan> plans;

    if (axis)
    {   
        fft_array_in = signals;
        for (int i = 0; i < signal_rows; i++)
        {
            plans.push_back(fftwf_plan_dft_1d(
                fft_size,
                reinterpret_cast<fftwf_complex*>(fft_array_in[i].data()),
                reinterpret_cast<fftwf_complex*>(fft_array_out[i].data()),
                fft_direction,
                FFTW_ESTIMATE
            ));
        }
    }
    else
    {
        fft_array_in.resize(signal_cols, vector<complex<float>>(signal_rows));
        for (int i = 0; i < signal_cols; i++)
        {
            for (int j = 0; j < signal_rows; j++)
            {
                fft_array_in[i][j] = signals[j][i];  // transpose the signal
            }
            plans.push_back(fftwf_plan_dft_1d(
                fft_size,
                reinterpret_cast<fftwf_complex*>(fft_array_in[i].data()),
                reinterpret_cast<fftwf_complex*>(fft_array_out[i].data()),
                fft_direction,
                FFTW_ESTIMATE
            ));
        }
    }

    if (inverse) cout << "Executing IFFTs" << endl;
    else         cout << "Executing FFTs"  << endl;

    for (fftwf_plan& plan : plans) fftwf_execute(plan);

    cout << "Destroying Plans" << endl;

    for (fftwf_plan& plan : plans) fftwf_destroy_plan(plan);

    cout << "Setting FFT Output Data" << endl;

    float norm_factor = inverse ? 1.0f / (fft_size) : 1.0f;

    if (not axis)
    {
        vector<vector<complex<float>>> array_out(fft_size, vector<complex<float>>(out_array_rows));
        for (int i = 0; i < fft_size; i++)
        {
            for (int j = 0; j < out_array_rows; j++)
            {
                array_out[i][j] = fft_array_out[j][i] * norm_factor;
            }
        }
        return array_out;
    }
    else if (inverse)
    {
        for (int i = 0; i < out_array_rows; i++)
        {
            for (int j = 0; j < fft_size; j++)
            {
                fft_array_out[i][j] *= norm_factor;
            }
        }
    }

    return fft_array_out;
}


vector<vector<complex<float>>> compute_2d_dft(
    const vector<vector<complex<float>>>& signal,
    const bool& inverse = false,
    int fft_rows = 0,
    int fft_cols = 0
) {
    cout << "Initializing 1D Complex Vector for FFTW" << endl;

    if (fft_rows == 0) fft_rows = signal.size();
    if (fft_cols == 0) fft_cols = signal[0].size();

    bool rows_out_of_lims = fft_rows < 0 or fft_rows > signal.size();
    bool cols_out_of_lims = fft_cols < 0 or fft_cols > signal[0].size();

    if (rows_out_of_lims or cols_out_of_lims)
    {
        throw invalid_argument("Invalid FFT size for signal.");
    }

    vector<complex<float>> signal_fftw(fft_rows * fft_cols);

    for (int i = 0; i < fft_rows; i++)
    {
        for (int j = 0; j < fft_cols; j++)
        {
            signal_fftw[i*fft_cols+j] = signal[i][j];
        }
    }

    fftwf_plan plan;

    if (inverse) cout << "Executing 2D IFFT Plan" << endl;
    else cout << "Executing 2D FFT Plan"  << endl;

    float norm_factor = inverse ? 1.0f / (fft_rows * fft_cols) : 1.0f;
    int fft_direction = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    plan = fftwf_plan_dft_2d(
        fft_rows, fft_cols, 
        reinterpret_cast<fftwf_complex*>(signal_fftw.data()),
        reinterpret_cast<fftwf_complex*>(signal_fftw.data()),
        fft_direction, FFTW_ESTIMATE
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
            signal_fft[i][j] = signal_fftw[i*fft_cols+j] * norm_factor;
        }
    }

    return signal_fft;
}
