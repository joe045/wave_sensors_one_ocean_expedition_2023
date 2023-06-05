# wave_sensors_one_ocean_expedition_2023
Processing of the wave data obtained during the One Ocean Expedition

Expedition plan https://oneoceanexpedition.com

Expedition data online http://metadata.nmdc.no/metadata-api/landingpage/16316eb3cea666a1871679a8b78568e1

# Sensor Composition
- 3 Altimeters
    * 2 ultrasonic altimeters: Automation Products Group Inc, IRU-3433-C60
    * 1 radar altimeter: Automation Products Group Inc, True EchoT M Pulse Radar Level Transmit- ter PRL-050-V024-C4-SS-S6-BF-F
- 3 IMUs:
    * 1 main IMU: VectorNav Technologies, VN-100
    * 2 extra IMUs: Adafruit, ISM330DHCX + LIS3MDL FeatherWing
- Processing unit: 
    * Raspberry Pi 4 Model 8GB
    * Arduino Due board
    * Redboard

# Timeline
- 16.06.2021: Wave system mounted on Statsraad Lehmkuhl
- 20.08.2021: Expedition start in Arendal
- 11.10.2021: extra IMUs mounted
- 14.01.2022: extra IMU1 destroyed 
- 05.05.2022: extra IMU0 destroyed and system breakdown
- 15.09.2022: Wave system fixed by replacing processing unit and main IMU
- 15.04.2023: Expedition end in Bergen

# Organisation
- Code includes processing of the data
   * wave_processing includes all functions needed to calculate the wave parameters
   * run_all_data goes trough each 30minute file, calculates the wave parameters and returns everything in a dictionary
- Additional data includes wind and the data needed for calculation of the doppler effect
   * all_pos.csv includes ship position, ship speed and heading, wind speed and direction with minutely resolution
   * wave_direction.csv includes wave direction with hourly resolution    
