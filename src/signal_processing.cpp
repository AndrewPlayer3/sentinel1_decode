#include "signal_processing.h"


CF_VEC_1D linspace(const std::complex<float>& start, const std::complex<float>& end, const int& size)
{
    CF_VEC_1D range(size);

    std::complex<float> step = (end - start) / (static_cast<float>(size) - 1);

    for (int i = 0; i < size; ++i) {
        range[i] = start + static_cast<float>(i) * step;
    }

    return range;
}


F_VEC_1D linspace(const float& start, const float& end, const int& size)
{
    F_VEC_1D range(size);

    float step = (end - start) / (size - 1);

    for (int i = 0; i < size; ++i) {
        range[i] = start + i * step;
    }

    return range;
}


CF_VEC_1D conjugate(const CF_VEC_1D& complex_samples)
{
    int num_samples = complex_samples.size();

    CF_VEC_1D complex_conj = complex_samples;

    conjugate_in_place(complex_conj);

    return complex_conj;
}   


void conjugate_in_place(CF_VEC_1D& complex_samples)
{
    std::for_each(
        complex_samples.begin(), complex_samples.end(),
            [] (std::complex<float>& n) { n = std::conj(n); }
    );
}  


F_VEC_1D hanning_window(const int& num_samples)
{
    F_VEC_1D window(num_samples);

    std::iota(window.begin(), window.end(), 0.0f);

    std::for_each(
        window.begin(), window.end(),
            [num_samples] (float& n) { n = sin(PI * n / num_samples) * sin(PI * n / num_samples); }
    );
    return window;
}


CF_VEC_1D apply_hanning_window(const CF_VEC_1D& complex_samples)
{
    int num_samples = complex_samples.size();

    F_VEC_1D  window = hanning_window(num_samples);
    CF_VEC_1D filtered_samples(num_samples);

    std::transform(
        complex_samples.begin(), complex_samples.end(), window.begin(), 
            filtered_samples.begin(),
                [] (const std::complex<float>& n, const float& m) { return n * m; }
    );
    return filtered_samples;
}


void apply_hanning_window_in_place(CF_VEC_1D& complex_samples)
{
    int num_samples = complex_samples.size();

    F_VEC_1D window = hanning_window(num_samples);

    std::transform(
        window.begin(), window.end(), 
            complex_samples.begin(), complex_samples.begin(),
                [] (float& w, std::complex<float>& n) { return n * w; }
    );
}


F_VEC_1D flatten(const F_VEC_2D& values)
{
    int rows = values.size();
    int cols = values[0].size();

    F_VEC_1D flat(rows * cols);

    for (int i = 0; i < rows; i++)
    {
        const F_VEC_1D& value_row = values[i];
        for (int j = 0; j < cols; j++)
        {
            flat[i * cols + j] = value_row[j];
        }
    }
    return flat;
}


CF_VEC_1D flatten(const CF_VEC_2D& values)
{
    int rows = values.size();
    int cols = values[0].size();

    CF_VEC_1D flat(rows * cols);

    for (int i = 0; i < rows; i++)
    {
        const CF_VEC_1D& value_row = values[i];
        for (int j = 0; j < cols; j++)
        {
            flat[i * cols + j] = value_row[j];
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
        const CF_VEC_1D& complex_value_row = complex_values[i];
        F_VEC_1D& norm_row = norm[i];

        for (int j = 0; j < cols; j++)
        {
            float mag = std::abs(complex_value_row[j]);
            norm_row[j] = mag;

            if (mag > max_value) max_value = mag;
        }
    }

    for (int i = 0; i < rows; i++)
    {
        F_VEC_1D& norm_row = norm[i];
        for (int j = 0; j < cols; j++)
        {
            norm_row[j] = log_scale ? 20 * log10(norm_row[j] / max_value) : norm_row[j] / max_value;
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
        float mag = std::abs(complex_values[i]);
        norm[i]   = mag;

        if (mag > max_value) max_value = mag;
    }

    std::for_each(
        norm.begin(), norm.end(), 
            [log_scale, max_value](float &n) 
            { 
                log_scale ? 20 * log10(n / max_value) : n / max_value;
            }
    );
    return norm;
}


F_VEC_1D magnitude_1d(const CF_VEC_1D& complex_values)
{
    F_VEC_1D magnitude(complex_values.size());

    std::transform(
        complex_values.begin(), complex_values.end(),
            magnitude.begin(),
                [] (const std::complex<float>& n) { return std::abs(n); }
    );
    return magnitude;
}


F_VEC_2D magnitude_2d(const CF_VEC_2D& complex_values)
{
    int rows = complex_values.size();
    int cols = complex_values[0].size();

    F_VEC_2D magnitude(rows, F_VEC_1D(cols));

    for (int i = 0; i < rows; i++)
    {
        const CF_VEC_1D& complex_values_row = complex_values[i];
               F_VEC_1D& magnitude_row = magnitude[i];

        std::transform(
            complex_values_row.begin(), complex_values_row.end(),
                magnitude_row.begin(),
                    [] (const std::complex<float>& n) { return std::abs(n); }
        );
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
        CF_VEC_1D& arr_t_row = arr_t[i];
        for (int j = 0; j < rows; j++)
        {
            arr_t_row[j] = arr[j][i];
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

    float norm_factor = 1.0f / static_cast<float>(fft_size);

    if (inverse)
    {
        std::for_each(
            fft_vector.begin(), fft_vector.end(), 
                [norm_factor](std::complex<float> &n) { n *= norm_factor; }
        );
    }
    return fft_vector;
}


void compute_1d_dft_in_place(
          CF_VEC_1D& signal,
          int   fft_size = 0,
    const bool& inverse  = false
) {
    if (fft_size <= 0) fft_size = signal.size();

    int fft_direction = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    fftwf_plan plan = fftwf_plan_dft_1d(
        fft_size,
        reinterpret_cast<fftwf_complex*>(signal.data()),
        reinterpret_cast<fftwf_complex*>(signal.data()),
        fft_direction,
        FFTW_ESTIMATE
    );
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    float norm_factor = 1.0f / static_cast<float>(fft_size);

    if (inverse)
    {
        std::for_each(
            signal.begin(), signal.end(), 
                [norm_factor](std::complex<float>& n) { n *= norm_factor; }
        );
    }
}


void compute_axis_dft_in_place(
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
        signals = transpose(signals);
        _compute_axis_dft(signals, fft_size, inverse);
        signals = transpose(signals);
    }
    else
    {
        _compute_axis_dft(signals, fft_size, inverse);
    }
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
        CF_VEC_2D signals_out = transpose(signals);
        _compute_axis_dft(signals_out, fft_size, inverse);
        return transpose(signals_out);
    }
    CF_VEC_2D signals_out = signals;
    _compute_axis_dft(signals_out, fft_size, inverse);
    return signals_out;
}


void _compute_axis_dft(
    CF_VEC_2D&  signals,
          int&  fft_size,
    const bool& inverse
) {
    int signal_rows   = signals.size();
    int signal_cols   = signals[0].size();
    int min_fft_size  = 8;
    int fft_direction = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    std::cout << "Initializing Plans and Excecuting Axis FFTs" << std::endl;
    for (int row = 0; row < signal_rows; row++)
    {
        CF_VEC_1D& signal_row = signals[row];
        fftwf_plan plan = fftwf_plan_dft_1d(
            fft_size,
            reinterpret_cast<fftwf_complex*>(signal_row.data()),
            reinterpret_cast<fftwf_complex*>(signal_row.data()),
            fft_direction,
            FFTW_ESTIMATE
        );
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
    }
    float norm_factor = inverse ? 1.0f / (fft_size) : 1.0f;

    if (inverse)
    {
        std::cout << "Normalizing Inverse Axis FFT Output" << std::endl;
        for (int i = 0; i < signal_rows; i++)
        {
            std::for_each(
                signals[i].begin(), signals[i].end(), 
                    [norm_factor](std::complex<float> &n) { n *= norm_factor; }
            );
        }
    }
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
