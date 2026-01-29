// e312_1809.cpp

#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/thread.hpp>
#include <zmq.hpp>
#include <boost/program_options.hpp>
#include <fftw3.h>
#include <chrono>
#include <thread>
#include <atomic>
#include <complex>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <csignal>
#include <vector>

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

struct BurstSteppedParams {
    // Frequency stepping parameters
    double center_freq = 5.4e9;         
    double step_bandwidth = 2e6;        // 2 MHz per step (reliable - max is ~3-4 MHz)
    double synthetic_bandwidth = 10e6;  // 10 MHz total synthetic bandwidth
    size_t num_freq_steps = 5;          // 5 × 2 MHz = 10 MHz
    
    // Burst parameters
    double burst_duration = 0.1;        // 100ms burst
    double burst_interval = 1.0;        // 1 second between burst starts
    double chirp_duration = 1e-3;       // 1ms chirps
    
    // Radio parameters
    double sample_rate = 2e6;           // 2 MS/s (matches step_bandwidth)
    double rx_gain = 30.0;				// Change this to 65-76 for long range applications (max is 76)
    double tx_gain = 20.0;				// Change this to 40-50 for long-range applications (max is 89.5)
    
    // Derived parameters
    size_t samples_per_chirp;
    size_t chirps_per_step;
    size_t total_chirps_per_burst;
    size_t range_bins_to_send = 256;
    size_t range_start_bin = 2;			// The step_bandwidth is 2e6, which has a range resolution of 75m. Therefore, range bin start of 2 excludes the first 150m, typically nearfield effects. Adjust this based on measurement distance to target.
    
    void update() {
        samples_per_chirp = static_cast<size_t>(sample_rate * chirp_duration);
        num_freq_steps = static_cast<size_t>(synthetic_bandwidth / step_bandwidth);
        
        // Divide burst time among frequency steps
        double time_per_step = burst_duration / num_freq_steps;
        chirps_per_step = static_cast<size_t>(time_per_step / chirp_duration);
        total_chirps_per_burst = chirps_per_step * num_freq_steps;
    }
    
    double get_freq_for_step(size_t step) {
        return center_freq - (synthetic_bandwidth/2) + (step * step_bandwidth);
    }
};

struct BurstHeader {
    uint32_t burst_id;
    uint32_t num_freq_steps;
    uint32_t chirps_per_step;
    uint32_t range_bins_per_chirp;
    float synthetic_bandwidth_mhz;
    uint64_t timestamp_us;
};

struct ChirpData {
    uint32_t freq_step;
    float center_freq_ghz;
};

class BurstSteppedCollector {
public:
    BurstSteppedCollector(const std::string& args, const std::string& pi_address)
        : context(1),
          socket(context, zmq::socket_type::push),
          pi_addr(pi_address) {
        
        // ZMQ setup
        int hwm = 100;
        socket.set(zmq::sockopt::sndhwm, hwm);
        int sndbuf = 8*1024*1024;
        socket.set(zmq::sockopt::sndbuf, sndbuf);
        
        std::string connect_addr = "tcp://" + pi_addr + ":5555";
        socket.connect(connect_addr);
        std::cout << "Connected to Pi at " << connect_addr << std::endl;
        
        // USRP setup
        std::string usrp_args = args.empty() ? "type=e3xx" : args;
        if (args.find("master_clock_rate") == std::string::npos) {
            usrp_args += ",master_clock_rate=16e6";  // Conservative MCR
        }
        
        usrp = uhd::usrp::multi_usrp::make(usrp_args);
        uhd::set_thread_priority_safe(1.0);
        
        std::cout << "USRP initialised for burst mode with frequency stepping" << std::endl;
    }
    
    ~BurstSteppedCollector() {
        if (fft_plan) fftwf_destroy_plan(fft_plan);
        if (fftw_in) fftwf_free(fftw_in);
        if (fftw_out) fftwf_free(fftw_out);
    }
    
    void configure(const BurstSteppedParams& params) {
        burst_params = params;
        burst_params.update();
        
        std::cout << "\n=== Burst + Frequency Stepping Configuration ===" << std::endl;
        std::cout << "Center frequency: " << params.center_freq/1e9 << " GHz" << std::endl;
        std::cout << "Burst duration: " << params.burst_duration*1000 << " ms" << std::endl;
        std::cout << "Burst interval: " << params.burst_interval << " s" << std::endl;
        std::cout << "Sample rate: " << params.sample_rate/1e6 << " MS/s" << std::endl;
        std::cout << "Frequency steps: " << params.num_freq_steps << std::endl;
        std::cout << "Step bandwidth: " << params.step_bandwidth/1e6 << " MHz" << std::endl;
        std::cout << "Synthetic bandwidth: " << params.synthetic_bandwidth/1e6 << " MHz" << std::endl;
        std::cout << "Chirps per step: " << burst_params.chirps_per_step << std::endl;
        std::cout << "Total chirps per burst: " << burst_params.total_chirps_per_burst << std::endl;
        
        // Configure radio
        usrp->set_clock_source("internal");
        usrp->set_time_source("internal");
        usrp->set_time_now(uhd::time_spec_t(0.0));
        
        // RX Configuration - start at center frequency
        usrp->set_rx_rate(burst_params.sample_rate, 0);
        usrp->set_rx_freq(uhd::tune_request_t(burst_params.center_freq), 0);
        usrp->set_rx_bandwidth(burst_params.step_bandwidth, 0);
        usrp->set_rx_antenna("RX2", 0);
        usrp->set_rx_gain(burst_params.rx_gain, 0);
        
        // TX Configuration - start at center frequency
        usrp->set_tx_rate(burst_params.sample_rate, 0);
        usrp->set_tx_freq(uhd::tune_request_t(burst_params.center_freq), 0);
        usrp->set_tx_bandwidth(burst_params.step_bandwidth, 0);
        usrp->set_tx_antenna("TX/RX", 0);
        usrp->set_tx_gain(burst_params.tx_gain, 0);
        
        // Stream args
        uhd::stream_args_t rx_args("fc32", "sc16");
        rx_args.channels = {0};
        rx_args.args["spp"] = "2000";
        rx_stream = usrp->get_rx_stream(rx_args);
        
        uhd::stream_args_t tx_args("fc32", "sc16");
        tx_args.channels = {0};
        tx_stream = usrp->get_tx_stream(tx_args);
        
        std::cout << "Actual sample rate: " << usrp->get_rx_rate(0)/1e6 << " MS/s" << std::endl;
        
        // Setup FFT
        setup_fft();
        
        // Generate reference chirp (single chirp, will be used for all frequencies)
        generate_reference_chirp();
        
        // Allocate buffers
        size_t max_samples = burst_params.total_chirps_per_burst * burst_params.samples_per_chirp;
        burst_buffer.resize(max_samples);
        compressed_burst.resize(burst_params.total_chirps_per_burst * burst_params.range_bins_to_send);
    }
    
    void start_collection(double total_duration) {
        std::cout << "\n=== Starting Burst Collection with Frequency Stepping ===" << std::endl;
        
        // Send metadata to Pi
        send_metadata(total_duration);
        std::cout << "Metadata sent to Pi" << std::endl;
        
        // Give Pi time to process metadata
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        
        auto collection_start = std::chrono::high_resolution_clock::now();
        size_t burst_count = 0;
        
        while (!stop_signal_called) {
            auto now = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = now - collection_start;
            if (elapsed.count() >= total_duration) break;
            
            auto burst_start = std::chrono::high_resolution_clock::now();
            
            std::cout << "\nBurst " << burst_count << ":" << std::endl;
            
            // Collect burst with TX/RX synchronisation
            if (collect_stepped_burst_synchronized()) {
                std::cout << "  Processing..." << std::endl;
                
                // Process burst (range compression)
                process_burst();
                
                std::cout << "  Sending..." << std::endl;
                
                // Send compressed data
                send_compressed_burst(burst_count);
                
                std::cout << "  Complete!" << std::endl;
                burst_count++;
            } else {
                std::cout << "  Failed!" << std::endl;
            }
            
            // Wait for next burst
            auto burst_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> burst_elapsed = burst_end - burst_start;
            
            double wait_time = burst_params.burst_interval - burst_elapsed.count();
            if (wait_time > 0) {
                std::cout << "  Waiting " << std::fixed << std::setprecision(1) 
                          << wait_time << "s for next burst..." << std::endl;
                std::this_thread::sleep_for(
                    std::chrono::milliseconds(static_cast<int>(wait_time * 1000))
                );
            }
        }
        
        send_end_marker();
        
        std::cout << "\nCollection complete: " << burst_count << " bursts" << std::endl;
    }
    
private:
    uhd::usrp::multi_usrp::sptr usrp;
    uhd::rx_streamer::sptr rx_stream;
    uhd::tx_streamer::sptr tx_stream;
    
    zmq::context_t context;
    zmq::socket_t socket;
    std::string pi_addr;
    
    fftwf_plan fft_plan = nullptr;
    fftwf_complex* fftw_in = nullptr;
    fftwf_complex* fftw_out = nullptr;
    
    BurstSteppedParams burst_params;
    
    std::vector<std::complex<float>> burst_buffer;
    std::vector<std::complex<float>> reference_chirp;
    std::vector<std::complex<float>> tx_buffer;
    std::vector<std::complex<float>> compressed_burst;
    std::vector<void*> tx_buffers;
    
    std::vector<size_t> freq_step_indices;
    
    void setup_fft() {
        size_t fft_size = next_power_of_2(burst_params.samples_per_chirp);
        
        fftw_in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fft_size);
        fftw_out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fft_size);
        
        fft_plan = fftwf_plan_dft_1d(fft_size, fftw_in, fftw_out,
                                     FFTW_FORWARD, FFTW_ESTIMATE);
        
        std::cout << "FFT ready: " << fft_size << " point FFT" << std::endl;
    }
    
    void generate_reference_chirp() {
        size_t n = burst_params.samples_per_chirp;
        reference_chirp.resize(n);
        tx_buffer.resize(n);
        
        double bw = burst_params.step_bandwidth;
        double t_chirp = burst_params.chirp_duration;
        
        for (size_t i = 0; i < n; i++) {
            double t = i / burst_params.sample_rate;
            double phase = M_PI * bw * t * t / t_chirp;
            reference_chirp[i] = std::complex<float>(cos(phase), sin(phase));
            tx_buffer[i] = 0.7f * reference_chirp[i];
        }
        
        tx_buffers = {tx_buffer.data()};
        std::cout << "Generated reference chirp" << std::endl;
    }
    
    bool collect_stepped_burst_synchronized() {
        // Synchronise TX and RX for each frequency step
        freq_step_indices.clear();
        size_t total_samples_collected = 0;
        
        // Collect data for each frequency step
        for (size_t step = 0; step < burst_params.num_freq_steps; step++) {
            // Tune both TX and RX to frequency for this step
            double freq = burst_params.get_freq_for_step(step);
            
            usrp->set_rx_freq(uhd::tune_request_t(freq), 0);
            usrp->set_tx_freq(uhd::tune_request_t(freq), 0);
            
            // Small delay for frequency settling
            std::this_thread::sleep_for(std::chrono::microseconds(100));
            
            std::cout << "  Step " << step << ": " << std::fixed << std::setprecision(3) 
                      << freq/1e9 << " GHz - ";
            
            // Transmit and receive chirps for this frequency step
            size_t samples_this_step = burst_params.chirps_per_step * burst_params.samples_per_chirp;
            
            // Start RX streaming for this step
            uhd::stream_cmd_t cmd(uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE);
            cmd.num_samps = samples_this_step;
            cmd.stream_now = true;
            rx_stream->issue_stream_cmd(cmd);
            
            // Transmit chirps while receiving
            uhd::tx_metadata_t tx_md;
            tx_md.start_of_burst = (step == 0);
            tx_md.end_of_burst = false;
            tx_md.has_time_spec = false;
            
            // TX and RX simultaneously for this step
            size_t samples_received = 0;
            size_t chirps_transmitted = 0;
            uhd::rx_metadata_t rx_md;
            
            while (samples_received < samples_this_step) {
                // Transmit chirp if needed
                if (chirps_transmitted < burst_params.chirps_per_step) {
                    tx_stream->send(tx_buffers, tx_buffer.size(), tx_md, 0);
                    tx_md.start_of_burst = false;
                    chirps_transmitted++;
                    
                    // Small delay between chirps
                    std::this_thread::sleep_for(std::chrono::microseconds(50));
                }
                
                // Receive samples
                size_t num_to_recv = std::min(burst_params.samples_per_chirp, 
                                             samples_this_step - samples_received);
                size_t num_rx = rx_stream->recv(
                    &burst_buffer[total_samples_collected + samples_received],
                    num_to_recv,
                    rx_md, 0.1
                );
                
                if (num_rx == 0 && rx_md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE) {
                    std::cout << "RX error: " << rx_md.strerror() << std::endl;
                    return false;
                }
                
                samples_received += num_rx;
            }
            
            // Track which samples belong to which frequency step
            for (size_t i = 0; i < burst_params.chirps_per_step; i++) {
                freq_step_indices.push_back(step);
            }
            
            total_samples_collected += samples_received;
            std::cout << chirps_transmitted << " chirps TX/RX" << std::endl;
        }
        
        // Send TX end of burst
        uhd::tx_metadata_t tx_md_end;
        tx_md_end.end_of_burst = true;
        tx_stream->send("", 0, tx_md_end);
        
        return true;
    }
    
    void process_burst() {
        size_t fft_size = next_power_of_2(burst_params.samples_per_chirp);
        size_t compressed_idx = 0;
        
        // Process each chirp in the burst
        for (size_t chirp = 0; chirp < burst_params.total_chirps_per_burst; chirp++) {
            size_t offset = chirp * burst_params.samples_per_chirp;
            size_t freq_step = freq_step_indices[chirp];
            
            // Dechirp with reference
            for (size_t i = 0; i < burst_params.samples_per_chirp; i++) {
                auto dechirped = burst_buffer[offset + i] * std::conj(reference_chirp[i]);
                fftw_in[i][0] = dechirped.real();
                fftw_in[i][1] = dechirped.imag();
            }
            
            // Zero pad
            for (size_t i = burst_params.samples_per_chirp; i < fft_size; i++) {
                fftw_in[i][0] = 0;
                fftw_in[i][1] = 0;
            }
            
            // FFT
            fftwf_execute(fft_plan);
            
            // Extract range bins
            for (size_t i = 0; i < burst_params.range_bins_to_send; i++) {
                size_t bin = burst_params.range_start_bin + i;
                compressed_burst[compressed_idx++] = 
                    std::complex<float>(fftw_out[bin][0], fftw_out[bin][1]);
            }
        }
    }
    
    void send_compressed_burst(size_t burst_id) {
        // Send burst header
        BurstHeader header;
        header.burst_id = burst_id;
        header.num_freq_steps = burst_params.num_freq_steps;
        header.chirps_per_step = burst_params.chirps_per_step;
        header.range_bins_per_chirp = burst_params.range_bins_to_send;
        header.synthetic_bandwidth_mhz = burst_params.synthetic_bandwidth / 1e6;
        header.timestamp_us = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now().time_since_epoch()).count();
        
        zmq::message_t header_msg(&header, sizeof(header));
        socket.send(header_msg, zmq::send_flags::sndmore);
        
        // Prepare burst data with chirp metadata
        size_t total_size = burst_params.total_chirps_per_burst * 
                           (sizeof(ChirpData) + burst_params.range_bins_to_send * sizeof(std::complex<float>));
        
        std::vector<uint8_t> burst_data(total_size);
        size_t offset = 0;
        
        for (size_t chirp = 0; chirp < burst_params.total_chirps_per_burst; chirp++) {
            // Add chirp metadata
            ChirpData chirp_meta;
            chirp_meta.freq_step = freq_step_indices[chirp];
            chirp_meta.center_freq_ghz = burst_params.get_freq_for_step(chirp_meta.freq_step) / 1e9;
            
            memcpy(&burst_data[offset], &chirp_meta, sizeof(ChirpData));
            offset += sizeof(ChirpData);
            
            // Add range-compressed data
            size_t data_start = chirp * burst_params.range_bins_to_send;
            size_t data_size = burst_params.range_bins_to_send * sizeof(std::complex<float>);
            memcpy(&burst_data[offset], &compressed_burst[data_start], data_size);
            offset += data_size;
        }
        
        zmq::message_t data_msg(burst_data.data(), burst_data.size());
        socket.send(data_msg, zmq::send_flags::none);
    }
    
    void send_metadata(double duration) {
        std::string meta = "METADATA\n";
        meta += "mode=burst_stepped\n";
        meta += "version=1.0\n";
        meta += "center_freq=" + std::to_string(burst_params.center_freq) + "\n";
        meta += "step_bandwidth=" + std::to_string(burst_params.step_bandwidth) + "\n";
        meta += "synthetic_bandwidth=" + std::to_string(burst_params.synthetic_bandwidth) + "\n";
        meta += "num_freq_steps=" + std::to_string(burst_params.num_freq_steps) + "\n";
        meta += "sample_rate=" + std::to_string(burst_params.sample_rate) + "\n";
        meta += "burst_duration=" + std::to_string(burst_params.burst_duration) + "\n";
        meta += "burst_interval=" + std::to_string(burst_params.burst_interval) + "\n";
        meta += "chirps_per_step=" + std::to_string(burst_params.chirps_per_step) + "\n";
        meta += "total_chirps_per_burst=" + std::to_string(burst_params.total_chirps_per_burst) + "\n";
        meta += "range_bins=" + std::to_string(burst_params.range_bins_to_send) + "\n";
        meta += "range_start_bin=" + std::to_string(burst_params.range_start_bin) + "\n";
        meta += "total_duration=" + std::to_string(duration) + "\n";
        
        zmq::message_t msg(meta.data(), meta.size());
        socket.send(msg, zmq::send_flags::none);
    }
    
    void send_end_marker() {
        BurstHeader end_marker = {0xFFFFFFFF, 0, 0, 0, 0, 0};
        zmq::message_t msg(&end_marker, sizeof(end_marker));
        socket.send(msg, zmq::send_flags::none);
        std::cout << "End marker sent" << std::endl;
    }
    
    size_t next_power_of_2(size_t n) {
        size_t power = 1;
        while (power < n) power *= 2;
        return power;
    }
};

int UHD_SAFE_MAIN(int argc, char *argv[]) {
    signal(SIGINT, &sig_int_handler);
    
    std::string device_args, pi_address;
    double duration, freq, synthetic_bw, step_bw;
    double burst_dur, burst_int, rx_gain, tx_gain;
    
    po::options_description desc("Options");
    desc.add_options()
        ("help", "Show help")
        ("args", po::value<std::string>(&device_args)->default_value(""),
         "USRP args")
        ("pi-addr", po::value<std::string>(&pi_address)->default_value("192.168.1.100"),
         "Pi IP address")
        ("duration", po::value<double>(&duration)->default_value(30.0),
         "Total duration (s)")
        ("freq", po::value<double>(&freq)->default_value(5.4e9),
         "Center frequency (Hz) - default 5.4 GHz")
        ("synthetic-bw", po::value<double>(&synthetic_bw)->default_value(10e6),
         "Total synthetic bandwidth (Hz)")
        ("step-bw", po::value<double>(&step_bw)->default_value(2e6),
         "Bandwidth per frequency step (Hz)")
        ("burst-duration", po::value<double>(&burst_dur)->default_value(0.1),
         "Duration of each burst (s)")
        ("burst-interval", po::value<double>(&burst_int)->default_value(1.0),
         "Time between burst starts (s)")
        ("rx-gain", po::value<double>(&rx_gain)->default_value(30.0),
         "RX gain (dB)")
        ("tx-gain", po::value<double>(&tx_gain)->default_value(20.0),
         "TX gain (dB)");
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        std::cout << "E312 Burst Mode with Frequency Stepping" << std::endl;
        std::cout << desc << std::endl;
        std::cout << "\nExample (2 MHz steps, 10 MHz synthetic):" << std::endl;
        std::cout << "  ./e312_burst_stepped --freq 5.4e9 --synthetic-bw 10e6 --step-bw 2e6" << std::endl;
        return 0;
    }
    
    try {
        BurstSteppedCollector collector(device_args, pi_address);
        
        BurstSteppedParams params;
        params.center_freq = freq;
        params.synthetic_bandwidth = synthetic_bw;
        params.step_bandwidth = step_bw;
        params.sample_rate = step_bw;  // Must match to satisfy Shannon-Nyquist
        params.burst_duration = burst_dur;
        params.burst_interval = burst_int;
        params.rx_gain = rx_gain;
        params.tx_gain = tx_gain;
        
        collector.configure(params);
        collector.start_collection(duration);
        
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
