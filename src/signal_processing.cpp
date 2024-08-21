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


vector<vector<complex<float>>> transpose(const vector<vector<complex<float>>>& arr)
{
    int rows = arr.size();
    int cols = arr[0].size();

    vector<vector<complex<float>>> arr_t(cols, vector<complex<float>>(rows));

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < cols; i++)
    {
        for (int j = 0; j < rows; j++)
        {
            arr_t[i][j] = arr[j][i];
        }
    }
    return arr_t;
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
    if (fft_size == 0)
    {
        fft_size = axis ? signals[0].size() : signals.size();
    }
    if (not axis)
    {
        cout << "Tranposing Vector for Row-Axis 1D FFT" << endl;;

        return transpose(_compute_1d_dft(transpose(signals), fft_size, inverse));
    }
    return _compute_1d_dft(signals, fft_size, inverse);
}


vector<vector<complex<float>>> _compute_1d_dft(const vector<vector<complex<float>>>& signals, int fft_size, const bool& inverse)
{
    int signal_rows   = signals.size();
    int signal_cols   = signals[0].size();
    int min_fft_size  = 8;
    int fft_direction = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    cout << "Initializing Complex Arrays" << endl;

    vector<vector<complex<float>>> fft_array_in = signals;
    vector<vector<complex<float>>> fft_array_out(signal_rows, vector<complex<float>>(fft_size));
    vector<fftwf_plan> plans(signal_rows);

    cout << "Initializing Plans, Excecuting FFTs" << endl;

    for (int row = 0; row < signal_rows; row++)
    {
        fftwf_plan plan = fftwf_plan_dft_1d(
            fft_size,
            reinterpret_cast<fftwf_complex*>(fft_array_in[row].data()),
            reinterpret_cast<fftwf_complex*>(fft_array_out[row].data()),
            fft_direction,
            FFTW_ESTIMATE
        );
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
    }

    cout << "Setting FFT Output Data" << endl;

    float norm_factor = inverse ? 1.0f / (fft_size) : 1.0f;

    if (inverse)
    {
        for (int i = 0; i < signal_rows; i++)
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
