# wave_sensors_one_ocean_expedition_2023
Processing of the wave data obtained during the One Ocean Expedition

Expedition plan https://oneoceanexpedition.com

Expedition data online http://metadata.nmdc.no/metadata-api/landingpage/16316eb3cea666a1871679a8b78568e1

# Sensor Composition
- 3 Altimeters
    * 2 ultrasonic altimeters: Automation Products Group Inc, IRU-3433-C60 https://www.apgsensors.com/sites/default/files/manuals/IRU-
    * 1 radar altimeter: Automation Products Group Inc, True EchoT M Pulse Radar Level Transmitter PRL-050-V024-C4-SS-S6-BF-F www.apgsensors.com/sites/default/files/manuals/PRL-PRS-manual.pdf
 
[system.pdf](https://github.com/joe045/wave_sensors_one_ocean_expedition_2023/files/12301219/system.pdf)

![ug1_bowsprit](https://github.com/joe045/wave_sensors_one_ocean_expedition_2023/assets/71880331/edf7663e-067c-4064-9438-65401e83e127)
![radar_starboard](https://github.com/joe045/wave_sensors_one_ocean_expedition_2023/assets/71880331/d49d8cfa-7ff1-4a91-8ba3-0fbcccb2bc25)
![radar_and_ug2](https://github.com/joe045/wave_sensors_one_ocean_expedition_2023/assets/71880331/a33bd630-95cc-4270-937e-a2b61dcb6a10)



- 3 IMUs:
    * 1 main IMU: VectorNav Technologies, VN-100 https://www.vectornav.com/products/detail/vn-
    * 2 extra IMUs: Adafruit, ISM330DHCX + LIS3MDL FeatherWing https://www.adafruit.com/product/4569
- Processing unit: 
    * Raspberry Pi 4 Model B 8GB datasheets.raspberrypi.com/rpi4/raspberry-pi-4-product-brief.pdf
    * Arduino Due board https://store.arduino.cc/products/arduino-due
    * Redboard https://www.sparkfun.com/products/15444

![new_processing_unit](https://github.com/joe045/wave_sensors_one_ocean_expedition_2023/assets/71880331/ec16c3ae-4ba4-4eb7-bdd8-cacc0600be8e)
[data_flow.pdf](https://github.com/joe045/wave_sensors_one_ocean_expedition_2023/files/12301222/data_flow.pdf)


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
   * run_all_data goes through each 30-minute file, calculates the wave parameters, and returns everything in a dictionary
   * arduino_logging.py contains a set of simple helper functions that allow to fully log the altimeters arduino output. This allows in       particular to easily: 1) find the arduino due logger, 2) grab packets for the gauges or the IMU from the arduino due logger, as          well as semaphores if gets out of alignment, 3) convert the ADC readings into actual distance measurements 4) format data for            sharing on UDP broadcasting network.
- Additional data includes additional data needed for the calculation of the Doppler effect
   * all_pos.csv includes ship position, ship speed and heading, wind speed and direction recorded at the ship with minute resolution
   * wave_direction.csv includes wave direction from the wave model ECWAM with hourly resolution    
