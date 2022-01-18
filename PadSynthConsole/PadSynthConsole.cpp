// PadSynthConsole.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <sstream>

#include "MyPadSynth.h"
#include <format>
#include <termcolor.hpp>

void init_parameters(std::vector<MyPadSynth::PadSynthParameters>& all_params) {
    float reference_frequency = 440.;
    const float base_frequency = 440. / 5;
    const int sample_rate = 44100;
    const int duration = sample_rate * 6;
    const int number_harmonics = 64;
    all_params.push_back(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, reference_frequency, std::pair<int, int>(1, 1), base_frequency, "A4"));

    all_params.push_back(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, reference_frequency, std::pair<int, int>(1, 5), base_frequency, "F2"));
    all_params.push_back(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, reference_frequency, std::pair<int, int>(9, 8 * 5), base_frequency, "G2"));
    all_params.push_back(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, reference_frequency, std::pair<int, int>(1, 4), base_frequency, "A2"));
    all_params.push_back(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, reference_frequency, std::pair<int, int>(3, 2 * 5), base_frequency, "C3"));
    all_params.push_back(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, reference_frequency, std::pair<int, int>(3 * 9, 2 * 8 * 5), base_frequency, "D3_pythagorean_sixth"));
    all_params.push_back(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, reference_frequency, std::pair<int, int>(2, 5), base_frequency, "F3"));
    all_params.push_back(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, reference_frequency, std::pair<int, int>(9 * 2, 8 * 5), base_frequency, "G3"));
    all_params.push_back(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, reference_frequency, std::pair<int, int>(1, 2), base_frequency, "A3"));
    all_params.push_back(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, reference_frequency, std::pair<int, int>(3 * 2, 2 * 5), base_frequency, "C4"));
    all_params.push_back(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, reference_frequency, std::pair<int, int>(3 * 9 * 2, 2 * 8 * 5), base_frequency, "D4_pythagorean_sixth"));
    all_params.push_back(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, reference_frequency, std::pair<int, int>(4, 5), base_frequency, "F4"));
    all_params.push_back(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, reference_frequency, std::pair<int, int>(4 * 9, 5 * 8), base_frequency, "G4"));
}

void text_file_init_parameters(std::vector<MyPadSynth::PadSynthParameters>& all_params) {
    std::string line;
    std::ifstream myfile("example.txt");
    const int sample_rate = 44100;
    const int duration = sample_rate * 6;
    const int number_harmonics = 64;
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            std::istringstream iline(line);
            std::string field;
            std::vector<float> params_raw;
            if (line.rfind("//", 0) == 0) continue;
            while (getline(iline, field, ' ')) {
                float fieldf = std::stof(field);
                params_raw.push_back(fieldf);
                //std::cout << termcolor::blue << field << termcolor::reset << std::endl;
                //if (fieldf) std::cout << termcolor::bright_blue << fieldf << termcolor::reset << std::endl;
            }
            try {
                float reference_frequency = params_raw.at(0);
                float base_frequency = params_raw.at(3);
                std::pair<int, int> ratio((int)params_raw.at(1), (int)params_raw.at(2));
                all_params.push_back(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, reference_frequency, ratio, base_frequency));
                std::cout << termcolor::green << std::format("reference freq {} Hz \n Ratio {}/{}", reference_frequency, ratio.first, ratio.second) << termcolor::reset << std::endl;
            }
            catch (const std::out_of_range& e) {
                std::cout << termcolor::bright_red << "Out of Range error." << termcolor::reset << std::endl;
                //std::cout << "Out of Range error." << std::endl;
            }
        }
        myfile.close();
    }
    else std::cout << termcolor::red << "Unable to open file" << termcolor::reset << std::endl;
    //else std::cout << "Unable to open file" << std::endl;
}

int main()
{
    //std::cout << "Hello World!\n";
    const int sample_rate = 44100;
    const int duration = sample_rate * 6;
    //float freq = 220.;
    float freq = 260.137f;
    const int number_harmonics = 64;
    MyPadSynth pad_synth(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, freq, std::pair<int, int>(1, 1), freq, "H7"));
    pad_synth.save_wavetable();
    //MyPadSynth pad_synth2(MyPadSynth::PadSynthParameters(duration, sample_rate, number_harmonics, 260.137f * 2, std::pair<int, int>(1, 1), 260.137f));
    //pad_synth2.save_wavetable();
    //MyPadSynth pad_synth3(duration, sample_rate, number_harmonics, 260.137f, std::pair<int, int>(2, 1), 260.137f);
    //pad_synth3.save_wavetable();

    std::vector<MyPadSynth::PadSynthParameters> all_params;
    init_parameters(all_params);
    

    //MyPadSynth::PadSynthParameters A110(duration, sample_rate, number_harmonics, 110, std::pair<int, int>(1, 1), 50);
    //all_params.push_back(A110);
    
    for (const MyPadSynth::PadSynthParameters& params : all_params) {
        //MyPadSynth pad_synth(params);
        pad_synth.set_parameters(params);
        pad_synth.save_wavetable();
    }

    
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
