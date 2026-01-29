// rpi_1809.cpp
// Receives data from the E312 via ethernet and saves to SSD

#include <zmq.hpp>
#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <thread>
#include <iomanip>
#include <csignal>
#include <cstring>
#include <ctime>
#include <complex>

static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

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

int main(int argc, char* argv[]) {
    signal(SIGINT, &sig_int_handler);
    
    // Parse arguments
    std::string save_path = "/mnt/ssd";
    int port = 5555;
    
    if (argc > 1) save_path = argv[1];
    if (argc > 2) port = std::atoi(argv[2]);
    
    std::cout << "\n=== RPi Burst+Stepped Receiver ===" << std::endl;
    std::cout << "Save path: " << save_path << std::endl;
    std::cout << "Port: " << port << std::endl;
    
    // Create ZMQ socket
    zmq::context_t context(1);
    zmq::socket_t socket(context, zmq::socket_type::pull);
    
    // Large receive buffer for burst data
    int rcvbuf = 16*1024*1024;  // 16MB
    socket.set(zmq::sockopt::rcvbuf, rcvbuf);
    
    // Bind
    std::string bind_addr = "tcp://*:" + std::to_string(port);
    socket.bind(bind_addr);
    std::cout << "Listening on " << bind_addr << std::endl;
    
    // Generate timestamp
    auto now = std::chrono::system_clock::now();
    auto t = std::chrono::system_clock::to_time_t(now);
    char timestamp[100];
    std::strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", std::localtime(&t));
    
    // Open output files
    std::string data_file = save_path + "/burst_stepped_" + std::string(timestamp) + ".dat";
    std::string meta_file = save_path + "/metadata_" + std::string(timestamp) + ".txt";
    std::string header_file = save_path + "/headers_" + std::string(timestamp) + ".dat";
    
    std::ofstream data_out(data_file, std::ios::binary);
    std::ofstream meta_out(meta_file);
    std::ofstream header_out(header_file, std::ios::binary);
    
    if (!data_out.is_open() || !meta_out.is_open() || !header_out.is_open()) {
        std::cerr << "ERROR: Cannot open output files!" << std::endl;
        return 1;
    }
    
    std::cout << "Output files:" << std::endl;
    std::cout << "  Data:    " << data_file << std::endl;
    std::cout << "  Headers: " << header_file << std::endl;
    std::cout << "  Meta:    " << meta_file << std::endl;
    
    // Wait for metadata
    std::cout << "\nWaiting for metadata (5 minute timeout)..." << std::endl;
    
    zmq::message_t meta_msg;
    socket.set(zmq::sockopt::rcvtimeo, 300000); // Time in milliseconds. Adjust this for longer timeout if needed.
    
    auto result = socket.recv(meta_msg);
    if (!result.has_value() || result.value() == 0) {
        std::cerr << "ERROR: Timeout waiting for metadata!" << std::endl;
        return 1;
    }
    
    // Save metadata
    std::string metadata(static_cast<char*>(meta_msg.data()), meta_msg.size());
    meta_out << "Recording started: " << timestamp << "\n";
    meta_out << metadata << std::endl;
    meta_out.flush();
    
    std::cout << "Metadata received!" << std::endl;
    
    // Parse key parameters from metadata
    size_t range_bins = 256;
    size_t num_freq_steps = 5;
    size_t chirps_per_step = 20;
    
    std::istringstream meta_stream(metadata);
    std::string line;
    while (std::getline(meta_stream, line)) {
        if (line.find("range_bins=") != std::string::npos) {
            range_bins = std::stoi(line.substr(line.find("=") + 1));
        } else if (line.find("num_freq_steps=") != std::string::npos) {
            num_freq_steps = std::stoi(line.substr(line.find("=") + 1));
        } else if (line.find("chirps_per_step=") != std::string::npos) {
            chirps_per_step = std::stoi(line.substr(line.find("=") + 1));
        }
    }
    
    size_t chirps_per_burst = num_freq_steps * chirps_per_step;
    size_t bytes_per_chirp = sizeof(ChirpData) + range_bins * sizeof(std::complex<float>);
    size_t expected_burst_size = chirps_per_burst * bytes_per_chirp;
    
    std::cout << "Expected configuration:" << std::endl;
    std::cout << "  Frequency steps: " << num_freq_steps << std::endl;
    std::cout << "  Chirps per step: " << chirps_per_step << std::endl;
    std::cout << "  Total chirps per burst: " << chirps_per_burst << std::endl;
    std::cout << "  Range bins: " << range_bins << std::endl;
    std::cout << "  Expected burst size: " << expected_burst_size/1024.0 << " KB" << std::endl;
    
    // Increase timeout for burst mode (bursts come every second)
    socket.set(zmq::sockopt::rcvtimeo, 5000);  // 5 second timeout
    
    // Receive bursts
    std::cout << "\nReceiving burst data..." << std::endl;
    
    size_t total_bursts = 0;
    size_t total_bytes = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    auto last_burst_time = start_time;
    
    while (!stop_signal_called) {
        // Receive burst header
        zmq::message_t header_msg;
        auto header_result = socket.recv(header_msg);
        
        if (!header_result.has_value() || header_result.value() == 0) {
            auto now = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> since_last = now - last_burst_time;
            
            if (total_bursts > 0 && since_last.count() > 10) {
                std::cout << "\nNo bursts for 10 seconds, ending reception" << std::endl;
                break;
            }
            continue;
        }
        
        // Check for end marker
        if (header_msg.size() == sizeof(BurstHeader)) {
            BurstHeader* header = static_cast<BurstHeader*>(header_msg.data());
            
            if (header->burst_id == 0xFFFFFFFF) {
                std::cout << "\nEnd marker received" << std::endl;
                break;
            }
            
            // Receive burst data
            zmq::message_t data_msg;
            auto data_result = socket.recv(data_msg);
            
            if (!data_result.has_value() || data_result.value() == 0) {
                std::cerr << "Failed to receive burst data!" << std::endl;
                continue;
            }
            
            // Write burst header
            header_out.write(reinterpret_cast<char*>(header), sizeof(BurstHeader));
            
            // Write burst data
            data_out.write(static_cast<char*>(data_msg.data()), data_msg.size());
            
            total_bursts++;
            total_bytes += sizeof(BurstHeader) + data_msg.size();
            last_burst_time = std::chrono::high_resolution_clock::now();
            
            // Progress update
            auto now = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = now - start_time;
            
            std::cout << "\rBurst " << std::setw(3) << header->burst_id 
                      << " | Total: " << total_bursts
                      << " | Data: " << std::fixed << std::setprecision(1) 
                      << total_bytes/1e6 << " MB"
                      << " | Time: " << std::setprecision(0)
                      << elapsed.count() << "s"
                      << " | Synthetic BW: " << header->synthetic_bandwidth_mhz << " MHz";
            std::cout.flush();
        }
    }
    
    // Close files
    data_out.close();
    header_out.close();
    
    // Final statistics
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    
    meta_out << "\n=== Recording Statistics ===" << std::endl;
    meta_out << "Total bursts: " << total_bursts << std::endl;
    meta_out << "Total chirps: " << total_bursts * chirps_per_burst << std::endl;
    meta_out << "Total data: " << total_bytes/1e6 << " MB" << std::endl;
    meta_out << "Duration: " << duration.count() << " seconds" << std::endl;
    
    if (duration.count() > 0) {
        double data_rate = (total_bytes * 8.0) / (duration.count() * 1e6);
        meta_out << "Average data rate: " << data_rate << " Mbps" << std::endl;
        
        double duty_cycle = (total_bursts * 0.1) / duration.count() * 100;
        meta_out << "Collection duty cycle: " << duty_cycle << "%" << std::endl;
    }
    
    meta_out.close();
    
    std::cout << "\n\n=== Reception Complete ===" << std::endl;
    std::cout << "Total bursts: " << total_bursts << std::endl;
    std::cout << "Total chirps: " << total_bursts * chirps_per_burst << std::endl;
    std::cout << "Total data: " << std::fixed << std::setprecision(2) 
              << total_bytes/1e6 << " MB" << std::endl;
    std::cout << "Duration: " << duration.count() << " seconds" << std::endl;
    std::cout << "\nFiles saved to: " << save_path << std::endl;
    std::cout << "Transfer to MATLAB for processing" << std::endl;
    
    return 0;
}
