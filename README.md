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
- 20.08.2021: Expedition start
- 11.10.2021: Extra IMUs mounted
- 14.01.2022: Extra IMU1 destroyed 
- 05.05.2022: Extra IMU0 destroyed and system breakdown
- 15.09.2022: Wave system fixed by replacing processing unit and Main IMU
- 15.04.2022: Expedition end

# Organisation
- Code includes processing of the data
   * pre_processing calculates the default values needed before running all data
   * wave_processing includes all functions needed to calculate the wave parameters
   * run_all_data goes trough each 30minute file, calculates the wave parameters and returns everything in a dictionary
- Data includes files from June 2021 - December 2022. The files are sorted into three periods and seperated at harbour and sailing. 
   * data includes all sailing days in the period October 2021 - December 2022
   * old_data includes all sailing and harbor days in the period June 2021 - August 2021
   * old_harbour includes all harbour days in the period September 2021 - Oktober 2021
   * harbour includes all harbour days in the period November 2021 - April 2022
   * new_harbour includes all harbour days in the period September 2022 - December 2022
   * bad_data includes some files that are not used
