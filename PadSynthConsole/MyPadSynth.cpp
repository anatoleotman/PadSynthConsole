/*
  ==============================================================================

    MyPadSynth.cpp
    Created: 22 Mar 2021 9:12:13pm
    Author:  jonny

  ==============================================================================
*/
//#define DBG ( textToWrite ) JUCE_BLOCK_WITH_FORCED_SEMICOLON (juce::String tempDbgBuf; tempDbgBuf << textToWrite; juce::Logger::outputDebugString (tempDbgBuf);)
#include "MyPadSynth.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <armadillo>
#include <format>
// http://www.csounds.com/manual/html/MiscFormants.html

MyPadSynth::MyPadSynth( int N,
                        int samplerate,
                        int number_harmonics,
                        float fundamental_frequency,
                        std::pair<int, int> frequency_ratio,
                        float base_frequency ):
    samplerate_(samplerate),
    number_harmonics_(number_harmonics),
    fundamental_frequency_(fundamental_frequency),
    frequency_ratio_(frequency_ratio),
    base_frequency_(base_frequency),
    A_(std::vector<float>(number_harmonics, 0)),
    freq_amp_(std::vector<float>(N/2, 0)),
    prototype_signal_(std::vector<float>(samplerate, 0)),
    params_(PadSynthParameters()),
    timbre_type_("")
{
    this->set_parameters(N, samplerate, number_harmonics, fundamental_frequency, frequency_ratio, base_frequency);
}

MyPadSynth::~MyPadSynth() {
    //qDebug() << "MyPadSynth::~MyPadSynth()";
};

float MyPadSynth::relF(int N) {
    return N;
};

void MyPadSynth::setharmonic(int n, float value) {
    if ((n < 1) || (n >= this->A_.size())) return;
    this->A_[n] = value;
};

float MyPadSynth::getharmonic(int n) {
    if ((n < 1) || (n >= this->A_.size())) return 0.0;
    return this->A_[n];
};

float MyPadSynth::profile(float fi, float bwi) {
    float x = fi / bwi;
    x *= x;
    if (x > 14.71280603) return 0.0;//this avoids computing the e^(-x^2) where it's results are very close to zero
    return exp(-x) / bwi;
};

void MyPadSynth::synth(float f1, float bandwidth, float bandwidth_scale, std::vector<float>& samples) {
    /*  synth() generates the wavetable
    f		- the fundamental frequency (eg. 440 Hz)
    bw		- bandwidth in cents of the fundamental frequency (eg. 25 cents)
    bwscale	- how the bandwidth increase on the higher harmonics (recomanded value: 1.0)
    *smp	- a pointer to allocated memory that can hold N samples */
    for (int i = 0; i < this->wavetable_size_ / 2; i++) this->freq_amp_[i] = 0.0;//default, all the frequency amplitudes are zero

    for (int nh = 1; nh < this->A_.size(); nh++) {//for each harmonic
        float bw_Hz;//bandwidth of the current harmonic measured in Hz
        float bwi;
        float fi;
        float rF = f1 * this->relF(nh);

        bw_Hz = (pow(2.0, bandwidth / 1200.0) - 1.0) * f1 * pow(this->relF(nh), bandwidth_scale);

        bwi = bw_Hz / (2.0 * this->samplerate_);
        fi = rF / this->samplerate_;
        for (int i = 0; i < this->wavetable_size_ / 2; i++) {//here you can optimize, by avoiding to compute the profile for the full frequency (usually it's zero or very close to zero)
            float hprofile;
            hprofile = this->profile((i / (float)this->wavetable_size_) - fi, bwi);
            this->freq_amp_[i] += hprofile * this->A_[nh];
        };
    };

    std::vector<float> freq_real(this->wavetable_size_ / 2);
    std::vector<float> freq_imaginary(this->wavetable_size_ / 2);

    //Convert the freq_amp array to complex array (real/imaginary) by making the phases random
    for (int i = 0; i < this->wavetable_size_ / 2; i++) {
        float phase = RND() * 2.0 * 3.14159265358979;
        freq_real[i] = this->freq_amp_[i] * cos(phase);
        freq_imaginary[i] = this->freq_amp_[i] * sin(phase);
    };
    this->IFFT(freq_real, freq_imaginary, samples);

    //normalize the output
    this->normalize(samples);
};

/*
    Simple normalization function. It normalizes the sound to 1/sqrt(2)
*/
void MyPadSynth::normalize(std::vector<float>& samples) {
    float max = 0.0;
    for (int i = 0; i < samples.size(); i++) {
        if (fabs(samples[i]) > max) {
            max = fabs(samples[i]);
        }
    }
    if (max < 1e-5) {
        max = 1e-5;
    }
    for (int i = 0; i < samples.size(); i++) {
        samples[i] /= max * 1.4142;
    }
}

std::vector<float> MyPadSynth::fourier_series_coefficients(int numbers_of_coefficients, int sample_rate, const std::vector<float>& samples)
{
    //"""Calculate the Fourier series coefficients up to the Nth harmonic"""
    //    result = []
    //    T = len(period)
    //    t = np.arange(T)
    //    for n in range(N + 1) :
    //        an = 2 / T * (period * np.cos(2 * np.pi * n * t / T)).sum()
    //        bn = 2 / T * (period * np.sin(2 * np.pi * n * t / T)).sum()
    //        result.append((an, bn))
    //        return np.array(result)
    std::vector<float> result;
    //std::vector<double> result_double;
    result.reserve(numbers_of_coefficients);
    //result_double.reserve(2 * numbers_of_coefficients);
    int T = sample_rate; // prototype signal must be of period number of samples == sample_rate
    //result.push_back(1.f);
    //result_double.push_back(1.);
    for (int n = 0; n != numbers_of_coefficients; ++n) {
        double an = 0;
        double bn = 0;
        double coefficient = 0;
        double cos_member = 0;
        double sin_member = 0;
        double sample = 0;
        for (int t = 0; t != T; ++t) {
            cos_member = cos(2 * (M_PI * n * t) / (T));
            sin_member = sin(2 * (M_PI * n * t) / (T));
            sample = (double)samples.at(t);
            an += sample * cos_member;
            bn += sample * sin_member;
        }
        an *= 2. / T;
        bn *= 2. / T;
        //result.push_back((float)an);
        //result.push_back((float)bn);
        coefficient = sqrt(pow(an, 2) + pow(bn, 2));
        result.push_back((float)coefficient);
        //result_double.push_back(an);
        //result_double.push_back(bn);
        std::cout << "MyPadSynth::fourier_series_coefficients: rank " << n << " coef: " << coefficient;
    }
    //DBG("MyPadSynth::fourier_series_coefficients" << result);
    return result;
}

void MyPadSynth::build_harmonics(MyPadSynth::timbre_type harmonics_profile)
{
    this->A_.clear();
    this->A_.resize(this->number_harmonics_);
    for (int i = 0; i < this->A_.size(); i++) this->A_[i] = 0.0;
    this->A_[1] = 1.0; //default, the first harmonic has the amplitude 1.0

    switch (harmonics_profile) {
    case timbre_type::formants:
        for (int i = 1; i < this->number_harmonics_; i++) {
            this->A_[i] = 1.0 / i;
            float formants = exp(-pow(((float)i * this->fundamental_frequency_ - 600.0) / 150.0, 2.0)) +
                                exp(-pow(((float)i * this->fundamental_frequency_ - 900.0) / 250.0, 2.0)) +
                                exp(-pow(((float)i * this->fundamental_frequency_ - 2200.0) / 200.0, 2.0)) +
                                exp(-pow(((float)i * this->fundamental_frequency_ - 2600.0) / 250.0, 2.0)) +
                                exp(-pow(((float)i * this->fundamental_frequency_) / 3000.0, 2.0)) * 0.1;
            this->A_[i] *= formants;

            //if ((i % 2) == 0) this->A_[i] *= 2.0;
        }
        break;
    case timbre_type::sawtooth:
        this->timbre_type_ = "sawtooth";
        for (int i = 1; i < this->number_harmonics_; i++) {
            this->A_[i] = 1.0 / i;
            if ((i % 2) == 0) this->A_[i] *= 2.0;
        }
        break;
    case timbre_type::square:
        int T = this->samplerate_;
        //std::vector<float> x;
        //x.reserve(T);
        //this->prototype_signal_.reserve(T);
        float f = 1.f;
        float x = 0;
        for (int t = 0; t != T; ++t) {
            x = 2 * tanh(sin(2 * M_PI * f * t / T));
            if (x >= 1) x = 1;
            if (x <= -1) x = -1;
            //x = sin(2 * juce::MathConstants<float>::pi * f * t / T);
            //if (x >= 0) x = 1;
            //else x = -1;
            this->prototype_signal_[t] = x;
        }
        this->A_.clear();
        this->A_ = this->fourier_series_coefficients(this->number_harmonics_, T, this->prototype_signal_);
        break;
    }
}

void MyPadSynth::resample_harmonics(float resampling_factor)
{
    int resampled_harmonics_amplitude_number = (int)round(this->A_.size() * resampling_factor);
    if (resampled_harmonics_amplitude_number != this->A_.size()) {
        // resample this->harmonics_amplitude_
        //qDebug() << "PadSynthAlgorithm::PadSynthAlgorithm resample harmonics_amplitude size" << resampled_harmonics_amplitude_number;
        std::vector<float> harmonics_amplitude = std::vector<float>(this->A_.begin(), this->A_.end());
        //    qDebug() << "PadSynthAlgorithm::PadSynthAlgorithm harmonics_amplitude" << harmonics_amplitude.size() << harmonics_amplitude;
        arma::frowvec harmonics_amplitude_arma = arma::frowvec(harmonics_amplitude);
        //    harmonics_amplitude_arma.print("harmonics_amplitude_arma");
        arma::frowvec harmonics_amplitude_locations = arma::linspace<arma::frowvec>(1, harmonics_amplitude_arma.size(), harmonics_amplitude_arma.size());
        //        harmonics_amplitude_locations.print("harmonics_amplitude_locations");
        arma::frowvec resampled_harmonics_amplitude_locations = arma::linspace<arma::frowvec>(1, harmonics_amplitude_arma.size(), resampled_harmonics_amplitude_number);
        arma::frowvec resampled_harmonics_amplitude_arma;
        arma::interp1(harmonics_amplitude_locations, harmonics_amplitude_arma, resampled_harmonics_amplitude_locations, resampled_harmonics_amplitude_arma, "*linear");
        std::vector<float> resampled_harmonics_amplitude_stdvec(arma::conv_to<std::vector<float>>::from(resampled_harmonics_amplitude_arma));
        std::vector<float>resampled_harmonics_amplitude(resampled_harmonics_amplitude_stdvec.begin(), resampled_harmonics_amplitude_stdvec.end());
        this->A_.clear();
        this->A_.reserve(resampled_harmonics_amplitude_number);
        this->A_.insert(this->A_.end(), resampled_harmonics_amplitude.begin(), resampled_harmonics_amplitude.end());
        this->number_harmonics_ = resampled_harmonics_amplitude_number;
    }
}

const std::vector<float>& MyPadSynth::getWavetable() const
{
    return wavetable_;
}

void MyPadSynth::set_parameters(int N, int samplerate, int number_harmonics, float reference_frequency, std::pair<int, int> frequency_ratio, float base_frequency)
{
    if (N) {
        auto start = std::chrono::high_resolution_clock::now();
        this->wavetable_size_ = N;
        this->samplerate_ = samplerate;
        this->number_harmonics_ = number_harmonics;
        this->reference_frequency_ = reference_frequency;
        this->fundamental_frequency_ = reference_frequency * (float)frequency_ratio.first / (float)frequency_ratio.second;
        this->frequency_ratio_ = frequency_ratio;
        this->base_frequency_ = base_frequency;
        //    float base_frequency = 130.81;
        //    float f1 = base_frequency+4;
        this->current_index_ = 0;

        //this->A_.resize(number_harmonics);
        //for (int i = 0; i < this->A_.size(); i++) this->A_[i] = 0.0;
        //this->A_[1] = 1.0; //default, the first harmonic has the amplitude 1.0
        //this->build_harmonics(timbre_type::square);
        //this->build_harmonics(timbre_type::formants);
        this->build_harmonics(timbre_type::sawtooth);
        this->freq_amp_.resize(N / 2);
        this->wavetable_.clear();
        this->wavetable_.resize(N);
        //    float f1 = 260.737;
        //    float f1 = 440;
        //    float f1 = 880;
        float resampling_factor = base_frequency / this->fundamental_frequency_;
        if (resampling_factor != 1.) {
            this->resample_harmonics(resampling_factor);
        }
        float bandwidth = 100; // bandwidth unit is cents ?
        float bandwidth_scale = 1.;
        this->harmonic_profile_bandwidth_ = bandwidth;

        this->synth(this->fundamental_frequency_,
                    bandwidth, 
                    bandwidth_scale, 
                    this->wavetable_);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout <<"Time taken by MyPadSynth::set_parameters: " << duration.count() * .000001 << " sec" << std::endl;
    }
}

void MyPadSynth::set_parameters(PadSynthParameters params)
{
    //if (params.N) {
    //    auto start = std::chrono::high_resolution_clock::now();
    //    this->wavetable_size_ = params.N;
    //    this->samplerate_ = params.samplerate;
    //    this->number_harmonics_ = params.number_harmonics;
    //    this->fundamental_frequency_ = params.fundamental_frequency;
    //    this->base_frequency_ = params.base_frequency;
    //    //    float base_frequency = 130.81;
    //    //    float f1 = base_frequency+4;
    //    this->current_index_ = 0;

    //    this->A_.resize(params.number_harmonics);
    //    for (int i = 0; i < this->A_.size(); i++) this->A_[i] = 0.0;
    //    this->A_[1] = 1.0; //default, the first harmonic has the amplitude 1.0
    //    //this->build_harmonics(timbre_type::square);
    //    this->build_harmonics(timbre_type::formants);
    //    this->freq_amp_.resize(params.N / 2);
    //    this->wavetable_.clear();
    //    this->wavetable_.resize(params.N);
    //    //    float f1 = 260.737;
    //    //    float f1 = 440;
    //    //    float f1 = 880;
    //    float resampling_factor = params.base_frequency / params.fundamental_frequency;
    //    if (resampling_factor != 1.) {
    //        this->resample_harmonics(resampling_factor);
    //    }
    //    float bandwidth = 40; // bandwidth unit is cents ?
    //    float bandwidth_scale = 1.;

    //    this->synth(params.fundamental_frequency,
    //                bandwidth,
    //                bandwidth_scale,
    //                this->wavetable_);

    //    auto stop = std::chrono::high_resolution_clock::now();
    //    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    //    std::cout <<"Time taken by MyPadSynth::set_parameters: " << duration.count() * .000001 << " sec" << std::endl;
    //}
    this->set_parameters(   params.N, 
                            params.samplerate, 
                            params.number_harmonics, 
                            params.reference_frequency,
                            params.frequency_ratio, 
                            params.base_frequency);
    this->params_ = params;
}

MyPadSynth::PadSynthParameters MyPadSynth::get_parameters() {
    MyPadSynth::PadSynthParameters params(this->wavetable_size_, this->samplerate_, this->number_harmonics_, this->reference_frequency_, this->frequency_ratio_, this->base_frequency_);
    return params;
}

void MyPadSynth::save_wavetable()
{
    //---------------------------------------------------------------
    // 1. Let's setup our AudioFile instance

    AudioFile<float> a;
    a.setNumChannels(1);
    a.setNumSamplesPerChannel(this->wavetable_.size());

    //---------------------------------------------------------------
    // 2. Create some variables to help us generate a sine wave

    //const float sampleRate = 44100.f;
    //const float frequencyInHz = 440.f;

    //---------------------------------------------------------------
    // 3. Write the samples to the AudioFile sample buffer

    //for (int i = 0; i < a.getNumSamplesPerChannel(); i++)
    //{
    //    for (int channel = 0; channel < a.getNumChannels(); channel++)
    //    {
    //        a.samples[channel][i] = sin((static_cast<float> (i) / sampleRate) * frequencyInHz * 2.f * M_PI);
    //    }
    //}

    a.samples[0] = this->wavetable_;
    //---------------------------------------------------------------
    // 4. Save the AudioFile
    std::string filePath = std::format("dossier/{}_{}Hz_r{}-{}_rf{}_nh{}_bf{}_bw{}_{}.wav",
                                        this->params_.note_name,
                                        this->fundamental_frequency_,
                                        this->frequency_ratio_.first,
                                        this->frequency_ratio_.second,
                                        this->reference_frequency_,
                                        this->number_harmonics_,
                                        this->base_frequency_,
                                        this->harmonic_profile_bandwidth_,
                                        this->timbre_type_);
    //std::string filePath = "dossier/formants_ff" +
    //                        std::to_string(this->reference_frequency_) +
    //                        "_r" + std::to_string(this->frequency_ratio_.first) +
    //                        "-" + std::to_string(this->frequency_ratio_.second) +
    //                        "_bf" + std::to_string(this->base_frequency_) +
    //                        "_nh" + std::to_string(this->number_harmonics_) + ".wav";
    //set_parameters(int N, int samplerate, int number_harmonics, float fundamental_frequency, float base_frequency)
    bool file_saved = a.save(filePath, AudioFileFormat::Wave);
    //std::cout << "**********************" << std::endl;
    if (file_saved) {
        std::cout << "Write Wavetable To Audio File " << filePath << std::endl;
    }
    else {
        std::cout << "Failed to Write Wavetable To Audio File" << std::endl;
    }
    //std::cout << "**********************" << std::endl << std::endl;
}

float MyPadSynth::RND() {
    return (rand() / (RAND_MAX + 1.0));
};

/*
    Inverse Fast Fourier Transform
    You may replace it with any IFFT routine
*/
void MyPadSynth::IFFT(std::vector<float>& freq_amp, std::vector<float>& freq_phase, std::vector<float>& samples) {
    FFTwrapper fft(this->wavetable_size_);
    FFTFREQS fftfreqs;
    newFFTFREQS(&fftfreqs, this->wavetable_size_ / 2);

    for (int i = 0; i < this->wavetable_size_ / 2; i++) {
        fftfreqs.c[i] = freq_amp[i] * cos(freq_phase[i]);
        fftfreqs.s[i] = freq_amp[i] * sin(freq_phase[i]);
    }

    fft.freqs2smps(fftfreqs, samples.data());
    deleteFFTFREQS(&fftfreqs);
}
