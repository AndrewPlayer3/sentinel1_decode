#include "signal_processing.h"


CF_VEC_1D conjugate(const CF_VEC_1D& complex_samples)
{
    int num_samples = complex_samples.size();

    CF_VEC_1D complex_conj(num_samples);

    for (int i = 0; i < num_samples; i++)
    {
        complex_conj[i] = complex_samples[i].real() - complex_samples[i].imag();
    }
    return complex_conj;
}   


F_VEC_1D hanning_window(const int& num_samples)
{
    F_VEC_1D window(num_samples);

    for (int i = 0; i < num_samples; i++)
    {
        window[i] = sin(PI * i / num_samples) * sin(PI * i / num_samples);
    }

    return window;
}


CF_VEC_1D apply_hanning_window(const CF_VEC_1D& complex_samples)
{
    int num_samples = complex_samples.size();

    F_VEC_1D  window = hanning_window(num_samples);
    CF_VEC_1D filtered_samples(num_samples);

    for (int i = 0; i < num_samples; i++)
    {
        filtered_samples[i] = complex_samples[i] * window[i];
    }
    return filtered_samples;
}


void apply_hanning_window_in_place(CF_VEC_1D& complex_samples)
{
    int num_samples = complex_samples.size();

    F_VEC_1D window = hanning_window(num_samples);

    for (int i = 0; i < num_samples; i++)
    {
        complex_samples[i] *= window[i];
    }
}


F_VEC_1D flatten(const F_VEC_2D& values)
{
    int rows = values.size();
    int cols = values[0].size();

    F_VEC_1D flat(rows * cols);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            flat[i * cols + j] = values[i][j];
        }
    }

    return flat;
}


CF_VEC_1D flatten(const CF_VEC_2D& values)
{
    int rows = values.size();
    int cols = values[0].size();

    CF_VEC_1D flat(rows * cols);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            flat[i * cols + j] = values[i][j];
        }
    }

    return flat;
}


F_VEC_2D norm_2d(
    const CF_VEC_2D& complex_values,
    const bool& log_scale = false
) {
    float max_value = 0.0;    

    int rows = complex_values.size();
    int cols = complex_values[0].size();

    F_VEC_2D norm(rows, F_VEC_1D(cols));

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


F_VEC_1D norm_1d(
    const CF_VEC_1D& complex_values,
    const bool& log_scale = false
) {
    int num_samples = complex_values.size();

    float max_value = 0;

    F_VEC_1D norm(num_samples);

    for (int i = 0; i < num_samples; i++)
    {
        float real = complex_values[i].real();
        float imag = complex_values[i].imag();

        float mag = sqrt(real*real + imag*imag);

        norm[i] = mag;

        if (mag > max_value) max_value = mag;
    }

    for (int i = 0; i < num_samples; i++)
    {
        norm[i] = log_scale ? 20 * log10(norm[i] / max_value) : norm[i] / max_value;
    }

    return norm;
}


F_VEC_1D magnitude_1d(const CF_VEC_1D& complex_values)
{
    int num_samples = complex_values.size();

    F_VEC_1D magnitude(num_samples);

    #pragma omp parallel for
    for (int i = 0; i < num_samples; i++)
    {
        float real = complex_values[i].real();
        float imag = complex_values[i].imag();

        float mag = sqrt(real*real + imag*imag);

        magnitude[i] = mag;    
    }

    return magnitude;
}


F_VEC_2D magnitude_2d(const CF_VEC_2D& complex_values)
{
    float max_value = 0;    

    int rows = complex_values.size();
    int cols = complex_values[0].size();

    F_VEC_2D magnitude(rows, F_VEC_1D(cols));

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            float real = complex_values[i][j].real();
            float imag = complex_values[i][j].imag();

            magnitude[i][j] = sqrt(real*real + imag*imag);
        }
    }

    return magnitude;
}


CF_VEC_2D transpose(const CF_VEC_2D& arr)
{
    int rows = arr.size();
    int cols = arr[0].size();

    CF_VEC_2D arr_t(cols, CF_VEC_1D(rows));

    for (int i = 0; i < cols; i++)
    {
        for (int j = 0; j < rows; j++)
        {
            arr_t[i][j] = arr[j][i];
        }
    }
    return arr_t;
}


CF_VEC_1D compute_1d_dft(
    const CF_VEC_1D& signal,
          int   fft_size = 0,
    const bool& inverse  = false
) {
    if (fft_size <= 0) fft_size = signal.size();

    CF_VEC_1D fft_vector = signal;

    int fft_direction = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    fftwf_plan plan = fftwf_plan_dft_1d(
        fft_size,
        reinterpret_cast<fftwf_complex*>(fft_vector.data()),
        reinterpret_cast<fftwf_complex*>(fft_vector.data()),
        fft_direction,
        FFTW_ESTIMATE
    );

    fftwf_execute(plan);

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


CF_VEC_2D compute_axis_dft(
    CF_VEC_2D&  signals,
          int   fft_size = 0,
    const int&  axis     = 0,
    const bool& inverse  = false
) {
    if (axis != 0 and axis != 1) 
    {
        throw std::invalid_argument("FFT axis must be 0 (rows) or 1 (cols).");
    }
    if (fft_size == 0)
    {
        fft_size = axis ? signals[0].size() : signals.size();
    }
    if (not axis)
    {
        std::cout << "Tranposing Vector for Row-Axis 1D FFT" << std::endl;;

        CF_VEC_2D signals_out = transpose(signals);
        signals_out = _compute_axis_dft(signals_out, fft_size, inverse);
        return transpose(signals_out);
    }
    return _compute_axis_dft(signals, fft_size, inverse);
}


CF_VEC_2D _compute_axis_dft(
    CF_VEC_2D&  signals,
          int&  fft_size,
    const bool& inverse
) {
    int signal_rows   = signals.size();
    int signal_cols   = signals[0].size();
    int min_fft_size  = 8;
    int fft_direction = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    std::cout << "Initializing Complex Arrays" << std::endl;

    CF_VEC_2D fft_array_out(signal_rows, CF_VEC_1D(fft_size));

    std::cout << "Initializing Plans, Excecuting FFTs" << std::endl;

    for (int row = 0; row < signal_rows; row++)
    {
        fftwf_plan plan = fftwf_plan_dft_1d(
            fft_size,
            reinterpret_cast<fftwf_complex*>(signals[row].data()),
            reinterpret_cast<fftwf_complex*>(fft_array_out[row].data()),
            fft_direction,
            FFTW_ESTIMATE
        );
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
    }

    std::cout << "Setting FFT Output Data" << std::endl;

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


CF_VEC_2D compute_2d_dft(
    const CF_VEC_2D& signal,
    const bool& inverse = false,
    int fft_rows = 0,
    int fft_cols = 0
) {
    std::cout << "Initializing 1D Complex Vector for FFTW" << std::endl;

    if (fft_rows == 0) fft_rows = signal.size();
    if (fft_cols == 0) fft_cols = signal[0].size();

    bool rows_out_of_lims = fft_rows < 0 or fft_rows > signal.size();
    bool cols_out_of_lims = fft_cols < 0 or fft_cols > signal[0].size();

    if (rows_out_of_lims or cols_out_of_lims)
    {
        throw std::invalid_argument("Invalid FFT size for signal.");
    }

    CF_VEC_1D signal_fftw(fft_rows * fft_cols);

    for (int i = 0; i < fft_rows; i++)
    {
        for (int j = 0; j < fft_cols; j++)
        {
            signal_fftw[i*fft_cols+j] = signal[i][j];
        }
    }

    fftwf_plan plan;

    if (inverse) std::cout << "Executing 2D IFFT Plan" << std::endl;
    else         std::cout << "Executing 2D FFT Plan"  << std::endl;

    float norm_factor = inverse ? 1.0f / (fft_rows * fft_cols) : 1.0f;
    int fft_direction = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    plan = fftwf_plan_dft_2d(
        fft_rows, fft_cols, 
        reinterpret_cast<fftwf_complex*>(signal_fftw.data()),
        reinterpret_cast<fftwf_complex*>(signal_fftw.data()),
        fft_direction, FFTW_ESTIMATE
    );

    fftwf_execute(plan);

    std::cout << "Destroying DFT Plan" << std::endl;

    fftwf_destroy_plan(plan);

    std::cout << "Copying Complex Data into 2D Vector" << std::endl;

    CF_VEC_2D signal_fft(fft_rows, CF_VEC_1D(fft_cols));

    for (int i = 0; i < fft_rows; i++)
    {
        for (int j = 0; j < fft_cols; j++)
        {
            signal_fft[i][j] = signal_fftw[i*fft_cols+j] * norm_factor;
        }
    }

    return signal_fft;
}


F_VEC_1D scale(const CF_VEC_1D& signal, const std::string& scaling_mode)
{
    F_VEC_1D samples(signal.size());

    if      (scaling_mode == "norm_log") samples = norm_1d(signal, true);
    else if (scaling_mode == "norm"    ) samples = norm_1d(signal, false);
    else if (scaling_mode == "mag"     ) samples = magnitude_1d(signal);      
    else if (scaling_mode == "real" or scaling_mode == "imag")
    {
        bool real = scaling_mode == "real";
        
        for (int i = 0; i < signal.size(); i++)
        {
            samples[i] = real ? signal[i].real() : signal[i].imag();
        }
    }
    else
    {
        throw std::invalid_argument(scaling_mode + " is not a valid scaling mode.");
    }
    return samples;
}


F_VEC_1D scale(const CF_VEC_2D& signal, const std::string& scaling_mode)
{
    int rows = signal.size();
    int cols = signal[0].size();
    
    F_VEC_1D samples(rows*cols);

    if      (scaling_mode == "norm_log") samples = flatten(norm_2d(signal, true));
    else if (scaling_mode == "norm"    ) samples = flatten(norm_2d(signal, false));
    else if (scaling_mode == "mag"     ) samples = flatten(magnitude_2d(signal));    
    else if (scaling_mode == "real" or scaling_mode == "imag")
    {
        bool real = scaling_mode == "real";
        int  size = signal.size() * signal[0].size();
        
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                samples[i*cols+j] = real ? signal[i][j].real() : signal[i][j].imag();
            }
        }
    }
    else
    {
        throw std::invalid_argument(scaling_mode + " is not a valid scaling mode.");
    }
    return samples;
}