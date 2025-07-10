# StARRS-Ground
Repository for the StARRS-Ground project software, including the drone radar system and Ground-Based Synthetic Aperture Radar (GB-SAR) system.

Each component has their own section, including:

# Drone Radar System

This section provides the Software for:

- **Raw Frequency Modulated Continuous Wave (FMCW) radar data collection** - Ettus Universal Software Radio Peripheral (USRP) Hardware Driver (UHD) E312 C++ scripts.
- **SAR post-processing** - Range Migration Algorithm (RMA) MATLAB script.
- **Radio Frequency (RF) triggering** - RF telemetry scripts for remote data triggering.
- **Global Navigation Satellite System (GNSS) corrections** - Overview of the EMLID software suite for precise platform positioning using Post Process Kinematics (PPK).
- **Interferometric Synthetic Aperture Radar (InSAR) post-processing** - Scripts for InSAR processing of drone radar data.

# GB-SAR System

This section provides the Software for:

- **System control** - Arduino scripts for linear rail operation (movement and speed control).
- **Raw radar data collection** - Ettus USRP UHD N200+CBX C++ and LibreVNA Python scripts.
- **SAR post-processing** - Scripts for SAR image formation using the linear rail.

# Laboratory Testing

This section provides the Software for:

- **Materials testing** - Overview of the LibreVNA software for laboratory radar materials testing.
