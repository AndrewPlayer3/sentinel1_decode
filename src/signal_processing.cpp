#include "signal_processing.h"


CF_VEC_1D linspace(const std::complex<double>& start, const std::complex<double>& end, const int& size)
{
    CF_VEC_1D range(size);

    std::complex<double> step = (end - start) / (static_cast<double>(size) - 1);

    for (int i = 0; i < size; ++i) {
        range[i] = start + static_cast<double>(i) * step;
    }

    return range;
}


F_VEC_1D linspace(const double& start, const double& end, const int& size)
{
    F_VEC_1D range(size);

    double step = (end - start) / (size - 1);

    for (int i = 0; i < size; ++i) {
        range[i] = start + i * step;
    }

    return range;
}


F_VEC_1D linear_resample(const F_VEC_1D& arr, const int& num_output_samples)
{
    int num_input_samples = arr.size();

    if (num_input_samples == num_output_samples)
    {
        return arr;
    }

    F_VEC_1D out(num_output_samples);

    double x = 0.0;
    double dx_in = 1.0 / (num_input_samples - 1);
    double dx_out = 1.0 / (num_output_samples - 1);

    int index = 0;

    for (int i = 0; i < num_input_samples; i++)
    {
        double x0 = i * dx_in;
        double x1 = (i + 1) * dx_in;
        double y0 = arr[i];
        double y1 = arr[i+1];

        while (x < x1 and index < num_output_samples)
        {
            out[index] = y0 + ((x - x0) * (y1 - y0) / (x1 - x0));
            x += dx_out;
            index += 1;
        }
    }

    out[-1] = out[-1];
 
    return out;
}


CF_VEC_1D linear_resample(const CF_VEC_1D& arr, const int& num_output_samples)
{
    int num_input_samples = arr.size();

    if (num_input_samples == num_output_samples)
    {
        return arr;
    }

    CF_VEC_1D out(num_output_samples);

    double x = 0.0;
    double dx_in = 1.0 / (num_input_samples - 1);
    double dx_out = 1.0 / (num_output_samples - 1);

    int index = 0;

    for (int i = 0; i < num_input_samples; i++)
    {
        double x0 = i * dx_in;
        double x1 = (i + 1) * dx_in;

        std::complex<double> y0 = arr[i];
        std::complex<double> y1 = arr[i+1];

        while (x < x1 and index < num_output_samples)
        {
            out[index] = y0 + ((x - x0) * (y1 - y0) / (x1 - x0));
            x += dx_out;
            index += 1;
        }
    }

    out[-1] = out[-1];
 
    return out;
}


CF_VEC_1D quadratic_resample(const CF_VEC_1D& arr, const int& num_output_samples)
{
    int num_input_samples = arr.size();

    if (num_input_samples == num_output_samples)
    {
        return arr;
    }

    if (num_input_samples < 3)
    {
        throw std::invalid_argument("Need at least 3 points for quadratic interpolation.");
    }

    CF_VEC_1D out(num_output_samples);

    double dx_in = 1.0 / (num_input_samples - 1);
    double dx_out = 1.0 / (num_output_samples - 1);

    for (int j = 0; j < num_output_samples; ++j)
    {
        double x = j * dx_out;

        int i = std::min(int(x / dx_in), num_input_samples - 2);
        
        int i0 = std::max(i - 1, 0);
        int i1 = i;
        int i2 = std::min(i + 1, num_input_samples - 1);

        double x0 = i0 * dx_in;
        double x1 = i1 * dx_in;
        double x2 = i2 * dx_in;

        std::complex<double> y0 = arr[i0];
        std::complex<double> y1 = arr[i1];
        std::complex<double> y2 = arr[i2];

        // Lagrange's form
        double L0 = (x - x1)*(x - x2)/((x0 - x1)*(x0 - x2));
        double L1 = (x - x0)*(x - x2)/((x1 - x0)*(x1 - x2));
        double L2 = (x - x0)*(x - x1)/((x2 - x0)*(x2 - x1));

        out[j] = L0*y0 + L1*y1 + L2*y2;
    }

    return out;
}


void conjugate_in_place(CF_VEC_1D& complex_samples)
{
    std::for_each(
        complex_samples.begin(), complex_samples.end(),
            [] (std::complex<double>& n) { n = std::conj(n); }
    );
}  



void apply_hanning_window_in_place(CF_VEC_1D& complex_samples)
{
    int num_samples = complex_samples.size();

    F_VEC_1D window(num_samples);

    std::iota(window.begin(), window.end(), 0.0);

    std::transform(
        window.begin(), window.end(), 
            complex_samples.begin(), complex_samples.begin(),
                [num_samples] (double& w, std::complex<double>& n) { 
                    return sin(PI * w / num_samples) * sin(PI * w / num_samples) * n;
                }
    );
}


void apply_hanning_window_in_place(CF_VEC_2D& complex_samples)
{
    int num_samples = complex_samples[0].size();

    F_VEC_1D window(num_samples);

    std::iota(window.begin(), window.end(), 0.0);

    for (CF_VEC_1D& row : complex_samples)
    {
        std::transform(
            window.begin(), window.end(), 
                row.begin(), row.begin(),
                    [num_samples] (double& w, std::complex<double>& n) { 
                        return sin(PI * w / num_samples) * sin(PI * w / num_samples) * n;
                    }
        );
    }
}


std::vector<float> flatten(const std::vector<std::vector<float>>& values)
{
    int rows = values.size();
    int cols = values[0].size();

    std::vector<float> flat(rows * cols);

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


std::vector<std::vector<float>> norm_2d(
    const CF_VEC_2D& complex_values,
    const bool& log_scale = false
) {
    double max_value = 0.0;    

    int rows = complex_values.size();
    int cols = complex_values[0].size();

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            max_value = std::max(max_value, std::abs(complex_values[i][j]));
        }
    }

    std::cout << "Normalization Value: " << max_value << std::endl;

    std::vector<std::vector<float>> norm(rows, std::vector<float>(cols));

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            double val = std::abs(complex_values[i][j]) + 1e-12;
            norm[i][j] = log_scale ? 20.0 * log10(val / max_value) : val / max_value;
        }
    }

    return norm;
}


F_VEC_1D norm_1d(
    const CF_VEC_1D& complex_values,
    const bool& log_scale = false
) {
    int num_samples = complex_values.size();
    double max_value = 0;

    F_VEC_1D norm(num_samples);

    for (int i = 0; i < num_samples; i++)
    {
        double mag = std::abs(complex_values[i]);
        norm[i]   = mag;

        if (mag > max_value) max_value = mag;
    }

    std::for_each(
        norm.begin(), norm.end(), 
            [log_scale, max_value](double &n) 
            { 
                return log_scale ? 20 * log10(n / max_value) : n / max_value;
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
                [] (const std::complex<double>& n) { return std::pow(std::abs(n), 2.0); }
    );
    return magnitude;
}


std::vector<std::vector<float>> magnitude_2d(const CF_VEC_2D& complex_values)
{
    int rows = complex_values.size();
    int cols = complex_values[0].size();

    std::vector<std::vector<float>> magnitude(rows, std::vector<float>(cols));

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            magnitude[i][j] = std::pow(std::abs(complex_values[i][j]), 2.0);
        }
    }
    return magnitude;
}


CF_VEC_2D transpose(const CF_VEC_2D& arr)
{
    int rows = arr.size();
    int cols = arr[0].size();

    CF_VEC_2D arr_t(cols, CF_VEC_1D(rows));

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


F_VEC_2D transpose(const F_VEC_2D& arr)
{
    int rows = arr.size();
    int cols = arr[0].size();

    F_VEC_2D arr_t(cols, F_VEC_1D(rows));

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


F_VEC_1D fftfreq(int n, double d = 1.0) {
    std::vector<double> freqs(n);
    double val = 1.0 / (n * d);
    int N = (n - 1) / 2 + 1;

    for (int i = 0; i < N; ++i)
        freqs[i] = i * val;

    for (int i = N; i < n; ++i)
        freqs[i] = (i - n) * val;

    return freqs;
}


CF_VEC_1D compute_1d_dft(
    const CF_VEC_1D& signal,
          int   fft_size = 0,
    const bool& inverse  = false
) {
    if (fft_size <= 0) fft_size = signal.size();

    CF_VEC_1D fft_vector = signal;

    int fft_direction = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    fftw_plan plan = fftw_plan_dft_1d(
        fft_size,
        reinterpret_cast<fftw_complex*>(fft_vector.data()),
        reinterpret_cast<fftw_complex*>(fft_vector.data()),
        fft_direction,
        FFTW_ESTIMATE
    );
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    double norm_factor = 1.0 / static_cast<double>(fft_size);

    if (inverse)
    {
        std::for_each(
            fft_vector.begin(), fft_vector.end(), 
                [norm_factor](std::complex<double> &n) { n *= norm_factor; }
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

    fftw_plan plan = fftw_plan_dft_1d(
        fft_size,
        reinterpret_cast<fftw_complex*>(signal.data()),
        reinterpret_cast<fftw_complex*>(signal.data()),
        fft_direction,
        FFTW_ESTIMATE
    );
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    double norm_factor = 1.0 / static_cast<double>(fft_size);

    if (inverse)
    {
        std::for_each(
            signal.begin(), signal.end(), 
                [norm_factor](std::complex<double>& n) { n *= norm_factor; }
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
    std::vector<fftw_plan> plans(signal_rows);

    std::transform(
        signals.begin(), signals.end(),
            plans.begin(),
                [fft_size, fft_direction] (CF_VEC_1D& row) {
                    return fftw_plan_dft_1d(
                        fft_size,
                        reinterpret_cast<fftw_complex*>(row.data()),
                        reinterpret_cast<fftw_complex*>(row.data()),
                        fft_direction,
                        FFTW_ESTIMATE
                    );
                }
    );

    #pragma omp parallel for
    for (int i = 0; i < signal_rows; i++)
    {
        fftw_execute(plans[i]);
    }

    std::for_each(
        plans.begin(), plans.end(),
            [] (fftw_plan& plan) { fftw_destroy_plan(plan); }
    );

    double norm_factor = inverse ? 1.0 / (fft_size) : 1.0;

    if (inverse)
    {
        std::cout << "Normalizing Inverse Axis FFT Output" << std::endl;
        #pragma omp parallel for num_threads(4)
        for (int i = 0; i < signal_rows; i++)
        {
            std::for_each(
                signals[i].begin(), signals[i].end(), 
                    [norm_factor](std::complex<double> &n) { n *= norm_factor; }
            );
        }
    }
}


std::vector<fftw_plan> get_fftw_plans(CF_VEC_2D& signals)
{
    std::vector<fftw_plan> plans(signals.size());
    std::transform(
        signals.begin(), signals.end(),
            plans.begin(),
                [] (CF_VEC_1D& row) {
                    return fftw_plan_dft_1d(
                        row.size(),
                        reinterpret_cast<fftw_complex*>(row.data()),
                        reinterpret_cast<fftw_complex*>(row.data()),
                        FFTW_FORWARD,
                        FFTW_ESTIMATE
                    );
                }
    );
    return plans;
}


void destroy_fftw_plans(std::vector<fftw_plan>& plans)
{
    std::for_each(
        plans.begin(), plans.end(),
            [] (fftw_plan& plan) { fftw_destroy_plan(plan); }
    );
}


void fftshift_in_place(CF_VEC_1D& signal)
{
    std::rotate(signal.begin(), signal.end()-(signal.size() / 2), signal.end());
}


void fftshift_in_place(CF_VEC_2D& signals)
{
    std::for_each(
        signals.begin(), signals.end(),
            [] (CF_VEC_1D& row) { 
                std::rotate(row.begin(), row.end()-(row.size() / 2), row.end());
            }
    );
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

    fftw_plan plan;

    if (inverse) std::cout << "Executing 2D IFFT Plan" << std::endl;
    else         std::cout << "Executing 2D FFT Plan"  << std::endl;

    double norm_factor = inverse ? 1.0 / (fft_rows * fft_cols) : 1.0;
    int fft_direction = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    plan = fftw_plan_dft_2d(
        fft_rows, fft_cols, 
        reinterpret_cast<fftw_complex*>(signal_fftw.data()),
        reinterpret_cast<fftw_complex*>(signal_fftw.data()),
        fft_direction, FFTW_ESTIMATE
    );

    fftw_execute(plan);

    std::cout << "Destroying DFT Plan" << std::endl;

    fftw_destroy_plan(plan);

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


std::vector<float> scale(const CF_VEC_2D& signal, const std::string& scaling_mode)
{
    int rows = signal.size();
    int cols = signal[0].size();

    std::vector<float> samples(rows*cols);

    std::cout << "First and Last 5 Before Scaling: " << std::endl;
    for (int i = 0; i < 5; i++) std::cout << signal[0][i] << " ";
    std::cout << std::endl;
    for (int i = 5; i > 0; i--) std::cout << signal[rows-1][cols-i] << " ";
    std::cout << std::endl;

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

    std::cout << "First and Last 5 After Scaling: " << std::endl;
    for (int i = 0; i < 5; i++) std::cout << samples[i] << " ";
    std::cout << std::endl;
    for (int i = 5; i > 0; i--) std::cout << samples[cols*rows-i] << " ";
    std::cout << std::endl;

    return samples;
}
