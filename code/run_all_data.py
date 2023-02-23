#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 12:06:18 2022

@author: fabianmk
@edited by: judiththuolberg
"""
import numpy as np
import pandas as pd
from datetime import datetime
import os

from wave_processing import *
#script that loads the raw data (needs arduino_logging.py) and does wave processing
 
def give_def_val(IMU_name,pitch_or_roll_or_yaw,IMU_pos):
    """Takes in IMU_name (IMU,extra_IMU_0,estraIMU_1), value (roll, pitch, yaw), IMU_pos (1-before 11.10.21, 2-before 5.5.22, 3-after 9.15.22 and returns the orientation parameters of the imus"""
    
    if IMU_pos==1:
        IMU_orientation="old"
    elif IMU_pos==2:
        IMU_orientation="new"
    else:
        IMU_orientation="newnew"
    
    
    def_dict={}
    def_dict["IMU"]={}
    def_dict["IMU"]["old"]={}
    def_dict["IMU"]["new"]={}
    def_dict["IMU"]["newnew"]={}
#    def_dict["IMU"]["old"]["roll"]=0.508
#    def_dict["IMU"]["old"]["pitch"]=-2.465
#    def_dict["IMU"]["old"]["yaw"]=3.7
    def_dict["IMU"]["old"]["roll"]=0.06792480580007744
    def_dict["IMU"]["old"]["pitch"]=-2.404658955566439
    def_dict["IMU"]["old"]["yaw"]=3.399999999995373
    def_dict["IMU"]["new"]["roll"]=-0.15671931142428283
    def_dict["IMU"]["new"]["pitch"]=-1.783959144381644
    def_dict["IMU"]["new"]["yaw"]=-3.0000000000042633
#    def_dict["IMU"]["new"]["roll"]=-0.24498794118246076
#    def_dict["IMU"]["new"]["pitch"]=-1.7619275213210193
#    def_dict["IMU"]["new"]["yaw"]=-2.7
    def_dict["IMU"]["newnew"]["roll"]=-0.1413391770025631
    def_dict["IMU"]["newnew"]["pitch"]=-2.4231007024933064
    def_dict["IMU"]["newnew"]["yaw"]=3.2999999999953786

    def_dict["extraIMU_0"]={}
    #def_dict["extraIMU_0"]["old"]={}
    def_dict["extraIMU_0"]["new"]={}
    #def_dict["extraIMU_0"]["old"]["pitch"]=-7.96277512364255
    #def_dict["extraIMU_0"]["old"]["roll"]=11.417037106673876
    #def_dict["extraIMU_0"]["old"]["yaw"]=30.11000000005958
    def_dict["extraIMU_0"]["new"]["roll"]=11.437288542801129
    def_dict["extraIMU_0"]["new"]["pitch"]=-7.981126028529418
    def_dict["extraIMU_0"]["new"]["yaw"]=31.899999999993753
    # def_dict["extraIMU_0"]["new"]["roll"]=11.419899780705149
    # def_dict["extraIMU_0"]["new"]["pitch"]=-7.948132384501466
    # def_dict["extraIMU_0"]["new"]["yaw"]=30.11000000005958
    
    def_dict["extraIMU_1"]={}
    #def_dict["extraIMU_1"]["old"]={}
    def_dict["extraIMU_1"]["new"]={}
    #def_dict["extraIMU_1"]["old"]["pitch"]=-30.961389387885745
    #def_dict["extraIMU_1"]["old"]["roll"]=11.053247661545374
    #def_dict["extraIMU_1"]["old"]["yaw"]=10.6 #given as 0?
    def_dict["extraIMU_1"]["new"]["roll"]=10.993627877532974
    def_dict["extraIMU_1"]["new"]["pitch"]=-30.970541641796522
    def_dict["extraIMU_1"]["new"]["yaw"]=-2.5000000000042917 #given as 0?
    # def_dict["extraIMU_1"]["new"]["roll"]=10.626786574787479
    # def_dict["extraIMU_1"]["new"]["pitch"]=-31.047273238809188
    # def_dict["extraIMU_1"]["new"]["yaw"]=0 
   
    return def_dict[IMU_name][IMU_orientation][pitch_or_roll_or_yaw]

def corr_pitch_n_roll_array(IMU_name,pitch_array,roll_array,time_array):
    """Takes IMU_name (IMU,extra_IMU_0,extraIMU_1), pitch, roll and time values as an array and returns values corrected for the imus orientation"""
        
    data1_pitch=pitch_array[np.where(time_array<datetime(2021,10,11,14,0))]
    data1_roll=roll_array[np.where(time_array<datetime(2021,10,11,14,0))]
    pitch1=data1_pitch-give_def_val(IMU_name,"pitch",1)
    roll1=data1_roll-give_def_val(IMU_name,"roll",1)
    alpha=give_def_val(IMU_name,"yaw",1)/360*(2*np.pi)
    if IMU_name == "IMU":
        if alpha >= 0:
            pitch2=np.cos(alpha)*pitch1+np.sin(alpha)*roll1
            roll2=np.cos(alpha)*roll1+np.sin(alpha)*(-pitch1)
        elif alpha <0:
            pitch2=np.cos(alpha)*pitch1+np.sin(alpha)*(-roll1)
            roll2=np.cos(alpha)*roll1+np.sin(alpha)*pitch1
        else:
            print("Alpha Alarm!")
        
        roll_ret=roll2
        pitch_ret=-pitch2
        alpha_ret=-alpha/(2*np.pi)*360
                
    elif IMU_name == "extraIMU_0" or IMU_name == "extraIMU_1":
        print('extra IMU value does not exist')
    else:
        print(" Wat isn dat fuern Name? Dat kann ich nich lesen.")
        
    #print('data1',data1_pitch.shape[0])    
    ret_roll=roll_ret
    ret_pitch=pitch_ret    
    
    data2_pitch=pitch_array[np.where((time_array>=datetime(2021,10,11,14,0)) & (time_array<datetime(2022,9,15,0,0)))]
    data2_roll=roll_array[np.where((time_array>=datetime(2021,10,11,14,0)) & (time_array<datetime(2022,9,15,0,0)))]
    pitch1=data2_pitch-give_def_val(IMU_name,"pitch",2)
    roll1=data2_roll-give_def_val(IMU_name,"roll",2)
    alpha=give_def_val(IMU_name,"yaw",2)/360*(2*np.pi)
    if IMU_name == "IMU":
        if alpha >= 0:
            pitch2=np.cos(alpha)*pitch1+np.sin(alpha)*roll1
            roll2=np.cos(alpha)*roll1+np.sin(alpha)*(-pitch1)
        elif alpha <0:
            pitch2=np.cos(alpha)*pitch1+np.sin(alpha)*(-roll1)
            roll2=np.cos(alpha)*roll1+np.sin(alpha)*pitch1
        else:
            print("Alpha Alarm!")
        
        roll_ret=roll2
        pitch_ret=-pitch2
        alpha_ret=-alpha/(2*np.pi)*360
        
        
    elif IMU_name == "extraIMU_0" or IMU_name == "extraIMU_1":
        if alpha >= 0:
            pitch2=np.cos(alpha)*pitch1+np.sin(alpha)*roll1
            roll2=np.cos(alpha)*roll1+np.sin(alpha)*(-pitch1)
        elif alpha <0:
            pitch2=np.cos(alpha)*pitch1+np.sin(alpha)*(-roll1)
            roll2=np.cos(alpha)*roll1+np.sin(alpha)*pitch1
        else:
            print("Alpha Alarm!")
        roll_ret=roll2
        pitch_ret=pitch2
        alpha_ret=alpha/(2*np.pi)*360
                
    else:
        print(" Wat isn dat fuern Name? Dat kann ich nich lesen.")
    
    #print('data2',data2_pitch.shape[0],data2_pitch) 
    if data2_pitch.shape[0]>0:    
        ret_roll=np.append(ret_roll,roll_ret)
        ret_pitch=np.append(ret_pitch,pitch_ret) 
        #ret_pitch=np.append(ret_roll,roll_ret) #Fabians code 
       

    data3_pitch=pitch_array[np.where(time_array>=datetime(2022,9,15,0,0))]
    data3_roll=roll_array[np.where(time_array>=datetime(2022,9,15,0,0))]
    pitch1=data3_pitch-give_def_val(IMU_name,"pitch",3)
    roll1=data3_roll-give_def_val(IMU_name,"roll",3)
    alpha=give_def_val(IMU_name,"yaw",3)/360*(2*np.pi)
    if IMU_name == "IMU":
        if alpha >= 0:
            pitch2=np.cos(alpha)*pitch1+np.sin(alpha)*roll1
            roll2=np.cos(alpha)*roll1+np.sin(alpha)*(-pitch1)
        elif alpha <0:
            pitch2=np.cos(alpha)*pitch1+np.sin(alpha)*(-roll1)
            roll2=np.cos(alpha)*roll1+np.sin(alpha)*pitch1
        else:
            print("Alpha Alarm!")
        
        roll_ret=roll2
        pitch_ret=-pitch2
        alpha_ret=-alpha/(2*np.pi)*360

    elif IMU_name == "extraIMU_0" or IMU_name == "extraIMU_1":
        print('extra IMU value does not exist')
    else:
        print(" Wat isn dat fuern Name? Dat kann ich nich lesen.")    
        
    #print('data3',data3_pitch.shape[0],data3_pitch) 
    if data3_pitch.shape[0]>0: 
        ret_roll=np.append(ret_roll,roll_ret)
        ret_pitch=np.append(ret_pitch,pitch_ret) 
        
    return (ret_pitch,ret_roll)


def simple_trig_approx_ug1(p,r,z):
    """Approximation of ug1's vertical position from extraIMU1 by the VN100's data. c-values are obtained from pre_processing"""
    #old
    #c0=18.68128014
    #c1=13.61833253
    #c2=-13.57785234
    #c3=-1.21498561
    #c4=-0.03734729
    #new old period
    #c0=18.464643800905158
    #c1=-7.776360620918963
    #c2=7.859373950422994
    #c3=0.798907688317781
    #c4=-0.005620148761905527
    #new longer period
    c0=19.561963453672995
    c1=-13.746961156354672
    c2=13.817934839355296
    c3=0.4058005160519371
    c4=-0.028129167092660497
    ret=c0*np.sin(-p)+c1*np.cos(r)+c2*np.cos(p)+c3*np.sin(r)+c4*z
    return ret


def simple_trig_approx_ug2(p,r,z):
    """Approximation of ug2's and radar's vertical position from extraIMU1 by the VN100's data. c-values are obtained from pre_processing"""
    #old
    #c0=3.84804107
    #c1=5.25491401
    #c2=-5.24400568
    #c3=-0.20843634
    #c4=-0.01711533
    #new old period
    #c0=3.18331459386499
    #c1=-7.037361316136021
    #c2=7.076888668816003
    #c3=0.5874789630595725
    #c4=0.01207023138877749
    #new longer period
    c0=5.120822763485035
    c1=-2.6625318968750946
    c2=2.6752341301980054
    c3=0.06497426749423409
    c4=-0.05877697081238737
    ret=c0*np.sin(-p)+c1*np.cos(r)+c2*np.cos(p)+c3*np.sin(r)+c4*z
    return ret

def kts_to_ms(kts):
    ms = kts*1852/(60*60)
    return ms

def correct_freq_for_Doppler_effect(freq,kts,heading,wave_dir):
    """Corrects the measured frequency for the Doppler Effect with the positive value. This is the function that should be used for correction"""
    g=9.81
    v_ship=kts_to_ms(kts)
    theta=wave_dir-heading
    v_o=v_ship*np.cos(theta/360*2*np.pi)
    #print('freq',freq)
    corr_freq_plus=(g/(4*np.pi*v_o)+np.sqrt(g/(2*np.pi*v_o)*(g/(8*np.pi*v_o)-freq)))
    corr_freq_minus=(g/(4*np.pi*v_o)-np.sqrt(g/(2*np.pi*v_o)*(g/(8*np.pi*v_o)-freq)))
    return corr_freq_plus,corr_freq_minus

def load_ship_stats(filepath="/Users/judiththuolberg/Master/Pos/all_pos.csv"):
    """Loads the ships positions, wind_direction, heading, wind and speed in knts from a .csv file and returns it in a dict"""
    df = pd.read_csv(filepath)
    ret_dict={}
    ret_dict["lon"]=df.lon.values
    ret_dict["lat"]=df.lat.values
    ret_dict["wind"]=df.Wind.values
    ret_dict["wind_dir"]=df["Wind dir"].values
    ret_dict["heading"]=df.Heading.values
    ret_dict["kts"]=df.Speed.values
    #make datetime object without pandas method 
    ret_dict["datetime"]=[]
    for timestring in df.datetime:
        date=datetime.strptime(str(timestring),"%Y-%m-%d %H:%M:%S")
        ret_dict["datetime"].append(date)
    ret_dict["datetime"]=np.array((ret_dict["datetime"]))
    return ret_dict

def load_wave_direction(filepath="/Users/judiththuolberg/Master/Pos/wave_direction.csv"):
    """Loads the ECWAM wave direction from a .csv file and returns it in a dict"""
    df = pd.read_csv(filepath)
    ret_dict={}
    ret_dict["model_lons"]=df.model_lons.values
    ret_dict["model_lats"]=df.model_lats.values
    ret_dict["mdir"]=df.mdir.values
    #make datetime object without pandas method 
    ret_dict["datetime"]=[]
    for timestring in df.datetime:
        date=datetime.strptime(str(timestring),"%Y-%m-%d %H:%M:%S")
        ret_dict["datetime"].append(date)
    ret_dict["datetime"]=np.array((ret_dict["datetime"]))
    return ret_dict

def get_wave_data_one_file(dict_data,ship_dict,wave_dict):
    """Calculates the important wave parameters from dict_data obtained from a single data file via load_data_dict(filepath), 
    and return all the relevant information on wave parameters, ship's position, wind, speed and heading"""
    
    #needed data from data-dict is taken and corrected
    pitch=np.array((dict_data["list_time_series_interpolated"]["IMU"]["pitch"]))
    #print('pitch',pitch)
    roll=np.array((dict_data["list_time_series_interpolated"]["IMU"]["roll"]))
    #print('roll',roll)
    date_time=np.array((dict_data["list_time_series_interpolated"]["common_datetime"]))
    #print('date_time',date_time)
    accel_D=np.array((dict_data["list_time_series_interpolated"]["IMU"]["accel_D"]))
    pitch, roll = corr_pitch_n_roll_array("IMU",pitch,roll,date_time)
    #print('lengths',len(pitch),len(roll),len(date_time))
    pitch=pitch/360*(2*np.pi)
    roll=roll/360*(2*np.pi)
    ug1=np.array(dict_data["list_time_series_interpolated"]["gauges"]["gauge_ug1_meters"])
    ug2=np.array(dict_data["list_time_series_interpolated"]["gauges"]["gauge_ug2_meters"])
    radar=np.array(dict_data["list_time_series_interpolated"]["gauges"]["gauge_radar1_meters"])

    #the vertical positions of the VN100 by two times integration
    sample_freq=10
    accel_D_mean=np.mean(accel_D)
    z_tt_VN100=-(accel_D-accel_D_mean)
    z_VN100=fft_int_2(z_tt_VN100,sample_freq)
    #print('lenz',len(z_VN100))
    
    #Due to Fourier transformations the data is cut off in beginning and end
    #10 er akkuratt nok til å ta med datasett med 1min intervall, mindre intervall klarer ikke regne heading, kts og wind dir
    data_len=pitch.shape[0]
    cutoff=data_len//10
    date_time=date_time[cutoff:-cutoff]
    pitch=pitch[cutoff:-cutoff]
    roll=roll[cutoff:-cutoff]
    z_VN100=z_VN100[cutoff:-cutoff]
    ug1=ug1[cutoff:-cutoff]
    #disse to linjene manglet i Fabians kode
    ug2=ug2[cutoff:-cutoff]
    radar=radar[cutoff:-cutoff]                    
    
    #keep these to get lists with equal length when using wrong roll value
    # cut_off_height = 14.75
    # date_time=date_time[np.where(ug1<cut_off_height)]
    # pitch=pitch[np.where(ug1<cut_off_height)]
    # roll=roll[np.where(ug1<cut_off_height)]
    # z_VN100=z_VN100[np.where(ug1<cut_off_height)]
    # ug1=ug1[np.where(ug1<cut_off_height)]
    # ug2=ug2[np.where(ug1<cut_off_height)]
    # radar=radar[np.where(ug1<cut_off_height)]
    data_used = pitch.shape[0]
    
    #the vertical positions of the sensors are computed
    delta_x_IMU_ug1=0.5
    delta_y_IMU_ug2=0.5
    z_corr_ug1=simple_trig_approx_ug1(pitch,roll,z_VN100)
    z_ug1=(z_VN100+z_corr_ug1)-np.mean(z_VN100+z_corr_ug1)
    z_corr_ug2=simple_trig_approx_ug2(pitch,roll,z_VN100)
    z_ug2=(z_VN100+z_corr_ug2)-np.mean(z_VN100+z_corr_ug2)
    pos_z_ug1=z_ug1#+np.sin(-pitch)*delta_x_IMU_ug1 #sligth improvement from old pos, barely noticable
    pos_z_ug2=z_ug2#+np.cos(pitch)*np.cos(roll)*delta_y_IMU_ug2 #sligth improvement from old pos, barely noticable
    
    #the water suface elevation WH_ossc is computed
    WH_ug1=pos_z_ug1-ug1*np.cos(pitch)*np.cos(roll)
    WH_mean_ug1=np.mean(WH_ug1)    
    WH_ossc_ug1=WH_ug1-WH_mean_ug1
    
    WH_ug2=pos_z_ug2-ug2*np.cos(pitch)*np.cos(roll)
    WH_mean_ug2=np.mean(WH_ug2)    
    WH_ossc_ug2=WH_ug2-WH_mean_ug2
    
    WH_radar=pos_z_ug2-radar*np.cos(pitch)*np.cos(roll)
    WH_mean_radar=np.mean(WH_radar)    
    WH_ossc_radar=WH_radar-WH_mean_radar
    
    #significant wave heights SWH are computed from standard deviation
    SWH_ug1=4*np.std(WH_ossc_ug1)
    SWH_ug2=4*np.std(WH_ossc_ug2)    
    SWH_radar=4*np.std(WH_ossc_radar)
        
    #significant wave heights H1/3 are computed manually
    WH_ug1_filt_test=WH_ossc_ug1*1
    cut_index_ug1=[]
    for pnt in range(WH_ug1_filt_test.shape[0]-3):
        if pnt < 3:
            continue
        if np.abs(WH_ug1_filt_test)[pnt] >= 1*SWH_ug1 and (not (np.abs(WH_ug1_filt_test)[pnt-1] >=1*SWH_ug1 and np.abs(WH_ug1_filt_test)[pnt+1] >=1*SWH_ug1 and np.abs(WH_ug1_filt_test)[pnt-2]>=1*SWH_ug1 and np.abs(WH_ug1_filt_test)[pnt+2]>=1*SWH_ug1 and np.abs(WH_ug1_filt_test)[pnt-3]>=SWH_ug1 and np.abs(WH_ug1_filt_test)[pnt+3]>=SWH_ug1)):
            cut_index_ug1.append(pnt)
    WH_ug1_filt_test=np.delete(WH_ug1_filt_test,cut_index_ug1)
    
    SWH_ug1_filt_test=4*np.std(WH_ug1_filt_test[np.where(WH_ug1_filt_test>=0)])

    WHs_ug1=[]
    upcross_indx=[]
    for cnt in range(WH_ug1_filt_test.shape[0]-1):
        if (WH_ug1_filt_test[cnt]<=0 and WH_ug1_filt_test[cnt+1]>0):
            upcross_indx.append(cnt)
            
    for cnt in range(len(upcross_indx)-1):
        subdata=WH_ug1_filt_test[upcross_indx[cnt]:upcross_indx[cnt+1]]*1
        top=np.max(subdata)
        low=np.min(subdata)
        WHs_ug1.append(top-low)
    WHs_ug1s=np.array((WHs_ug1))
    WHs_ug1s_sorted=np.sort(WHs_ug1)
    third_len=WHs_ug1s_sorted.shape[0]//3
    WHs_ug1s_top_third=WHs_ug1s_sorted[-third_len:]*1
    Hs_ug1_manual=np.mean(WHs_ug1s_top_third)
    
    WH_ug2_filt_test=WH_ossc_ug2*1
    cut_index_ug2=[]
    for pnt in range(WH_ug2_filt_test.shape[0]-3):
        if pnt < 3:
            continue
        if np.abs(WH_ug2_filt_test)[pnt] >= 1*SWH_ug2 and (not (np.abs(WH_ug2_filt_test)[pnt-1] >=1*SWH_ug2 and np.abs(WH_ug2_filt_test)[pnt+1] >=1*SWH_ug2 and np.abs(WH_ug2_filt_test)[pnt-2]>=1*SWH_ug2 and np.abs(WH_ug2_filt_test)[pnt+2]>=1*SWH_ug2 and np.abs(WH_ug2_filt_test)[pnt-3]>=SWH_ug2 and np.abs(WH_ug2_filt_test)[pnt+3]>=SWH_ug2)):
            cut_index_ug2.append(pnt)
    WH_ug2_filt_test=np.delete(WH_ug2_filt_test,cut_index_ug2)
    
    WH_ug2s=[]
    upcross_indx=[]
    for cnt in range(WH_ug2_filt_test.shape[0]-1):
        if (WH_ug2_filt_test[cnt]<=0 and WH_ug2_filt_test[cnt+1]>0):
            upcross_indx.append(cnt)
            
    for cnt in range(len(upcross_indx)-1):
        subdata=WH_ug2_filt_test[upcross_indx[cnt]:upcross_indx[cnt+1]]*1
        top=np.max(subdata)
        low=np.min(subdata)
        WH_ug2s.append(top-low)
    WH_ug2s=np.array((WH_ug2s))
    WH_ug2s_sorted=np.sort(WH_ug2s)
    third_len=WH_ug2s_sorted.shape[0]//3
    WH_ug2s_top_third=WH_ug2s_sorted[-third_len:]*1
    Hs_ug2_manual=np.mean(WH_ug2s_top_third)
    
    WH_radar_filt_test=WH_ossc_radar*1
    cut_index_radar=[]
    for pnt in range(WH_radar_filt_test.shape[0]-3):
        if pnt < 3:
            continue
        if np.abs(WH_radar_filt_test)[pnt] >= 1*SWH_radar and (not (np.abs(WH_radar_filt_test)[pnt-1] >=1*SWH_radar and np.abs(WH_radar_filt_test)[pnt+1] >=1*SWH_radar and np.abs(WH_radar_filt_test)[pnt-2]>=1*SWH_radar and np.abs(WH_radar_filt_test)[pnt+2]>=1*SWH_radar and np.abs(WH_radar_filt_test)[pnt-3]>=SWH_radar and np.abs(WH_radar_filt_test)[pnt+3]>=SWH_radar)):
            cut_index_radar.append(pnt)
    WH_radar_filt_test=np.delete(WH_radar_filt_test,cut_index_radar)
    
    WH_radars=[]
    upcross_indx=[]
    for cnt in range(WH_radar_filt_test.shape[0]-1):
        if (WH_radar_filt_test[cnt]<=0 and WH_radar_filt_test[cnt+1]>0):
            upcross_indx.append(cnt)
            
    for cnt in range(len(upcross_indx)-1):
        subdata=WH_radar_filt_test[upcross_indx[cnt]:upcross_indx[cnt+1]]*1
        top=np.max(subdata)
        low=np.min(subdata)
        WH_radars.append(top-low)
    WH_radars=np.array((WH_radars))
    WH_radars_sorted=np.sort(WH_radars)
    third_len=WH_radars_sorted.shape[0]//3
    WH_radars_top_third=WH_radars_sorted[-third_len:]*1
    Hs_radar_manual=np.mean(WH_radars_top_third)
    
    #calculating one value within the 30minutes datasetwith minute resolution
    t_start=date_time[0]
    t_end=date_time[-1]
    #print("t_start:",t_start,"t_end:",t_end)
    kts_array=ship_dict["kts"][np.where(((ship_dict["datetime"]>t_start)&(ship_dict["datetime"]<t_end)))]
    heading_array=ship_dict["heading"][np.where(((ship_dict["datetime"]>t_start)&(ship_dict["datetime"]<t_end)))]
    wave_dir_array=ship_dict["wind_dir"][np.where(((ship_dict["datetime"]>t_start)&(ship_dict["datetime"]<t_end)))]
    wind_array=ship_dict["wind"][np.where(((ship_dict["datetime"]>t_start)&(ship_dict["datetime"]<t_end)))]
    lon_array=ship_dict["lon"][np.where(((ship_dict["datetime"]>t_start)&(ship_dict["datetime"]<t_end)))]
    lat_array=ship_dict["lat"][np.where(((ship_dict["datetime"]>t_start)&(ship_dict["datetime"]<t_end)))]

    kts=np.mean(kts_array)
    heading=np.median(heading_array)   #need to calculate this better for values betwen 0 and 360
    wind_dir=np.median(wave_dir_array) #need to calculate this better for values betwen 0 and 360
    wind=np.mean(wind_array)
    lon=np.mean(lon_array)
    lat=np.mean(lat_array)
    
    #extracting one value within from the dataset with 60min resolution
    delta_t_start=np.abs(wave_dict["datetime"]-t_start)
    closest_t_start=np.min(delta_t_start)
    delta_t_end=np.abs(wave_dict["datetime"]-t_end)
    closest_t_end=np.min(delta_t_end)
    if closest_t_start<=closest_t_end:
        ind_closest_t=np.argmin(delta_t_start)
        wave_dir=wave_dict["mdir"][ind_closest_t]
    else:
        ind_closest_t=np.argmin(delta_t_end)
        wave_dir=wave_dict["mdir"][ind_closest_t]
        
    if wave_dir == 0:
        wave_dir +=0.1    
    
    if heading == 0:
        heading+=0.1
    
    if kts == 0:
        kts += 0.1
        
    #print("heading=",heading)
    #print("kts=",kts)
    #print("wave_dir=",wave_dir)
        
    dt_sec=(t_end-t_start).seconds
    times_180=dt_sec//180 #often same as length of kts_array = 20/5 so seg_secs = 180 is used
    #print(times_180)
    if times_180 >= 3:
        seg_secs=180
    else:
        seg_secs=dt_sec/2
    #print(seg_secs)    

    freqs_ug1,Pxx_ug1=welch_spectrum(WH_ossc_ug1,segment_duration_seconds=seg_secs)
    freqs_ug2,Pxx_ug2=welch_spectrum(WH_ossc_ug2,segment_duration_seconds=seg_secs)
    freqs_radar,Pxx_radar=welch_spectrum(WH_ossc_radar,segment_duration_seconds=seg_secs)
    
    Pxx_ug1=Pxx_ug1[np.where((freqs_ug1 >= 0.05)&(freqs_ug1<= 0.5))]
    Pxx_ug2=Pxx_ug2[np.where((freqs_ug2 >= 0.05)&(freqs_ug2<= 0.5))]
    Pxx_radar=Pxx_radar[np.where((freqs_radar >= 0.05)&(freqs_radar<= 0.5))]
    
    freqs_ug1=freqs_ug1[np.where((freqs_ug1 >= 0.05)&(freqs_ug1<= 0.5))]
    freqs_ug2=freqs_ug2[np.where((freqs_ug2 >= 0.05)&(freqs_ug2<= 0.5))]
    freqs_radar=freqs_radar[np.where((freqs_radar >= 0.05)&(freqs_radar<= 0.5))]
    
    Pxx_DoppShift_ug1_pos=Pxx_ug1*1
    Pxx_DoppShift_ug2_pos=Pxx_ug2*1
    Pxx_DoppShift_radar_pos=Pxx_radar*1
    
    Pxx_DoppShift_ug1_neg=Pxx_ug1*1
    Pxx_DoppShift_ug2_neg=Pxx_ug2*1
    Pxx_DoppShift_radar_neg=Pxx_radar*1
    
    freqs_DoppShift_ug1_pos=[]
    freqs_DoppShift_ug2_pos=[]
    freqs_DoppShift_radar_pos=[]
    
    freqs_DoppShift_ug1_neg=[]
    freqs_DoppShift_ug2_neg=[]
    freqs_DoppShift_radar_neg=[]
    
    for freq in freqs_ug1:
        freq_corr_pos,freq_corr_neg=correct_freq_for_Doppler_effect(freq,kts,heading,wave_dir)
        freqs_DoppShift_ug1_pos.append(freq_corr_pos)
        freqs_DoppShift_ug1_neg.append(freq_corr_neg)
    
    for freq in freqs_ug2:
        freq_corr_pos,freq_corr_neg=correct_freq_for_Doppler_effect(freq,kts,heading,wave_dir)
        freqs_DoppShift_ug2_pos.append(freq_corr_pos)
        freqs_DoppShift_ug2_neg.append(freq_corr_neg)  
    
    for freq in freqs_radar:
        freq_corr_pos,freq_corr_neg=correct_freq_for_Doppler_effect(freq,kts,heading,wave_dir)
        freqs_DoppShift_radar_pos.append(freq_corr_pos)
        freqs_DoppShift_radar_neg.append(freq_corr_neg) 
    
    #keep only non-nan values    
    freqs_DoppShift_ug1_pos=np.array((freqs_DoppShift_ug1_pos))
    Pxx_DoppShift_ug1_pos=Pxx_DoppShift_ug1_pos[np.isfinite(freqs_DoppShift_ug1_pos)]
    freqs_DoppShift_ug1_pos=freqs_DoppShift_ug1_pos[np.isfinite(freqs_DoppShift_ug1_pos)]
    
    freqs_DoppShift_ug1_neg=np.array((freqs_DoppShift_ug1_neg))
    Pxx_DoppShift_ug1_neg=Pxx_DoppShift_ug1_neg[np.isfinite(freqs_DoppShift_ug1_neg)]
    freqs_DoppShift_ug1_neg=freqs_DoppShift_ug1_neg[np.isfinite(freqs_DoppShift_ug1_neg)]
    
    freqs_DoppShift_ug2_pos=np.array((freqs_DoppShift_ug2_pos))
    Pxx_DoppShift_ug2_pos=Pxx_DoppShift_ug2_pos[np.isfinite(freqs_DoppShift_ug2_pos)]
    freqs_DoppShift_ug2_pos=freqs_DoppShift_ug2_pos[np.isfinite(freqs_DoppShift_ug2_pos)]
    
    freqs_DoppShift_ug2_neg=np.array((freqs_DoppShift_ug2_neg))
    Pxx_DoppShift_ug2_neg=Pxx_DoppShift_ug2_neg[np.isfinite(freqs_DoppShift_ug2_neg)]
    freqs_DoppShift_ug2_neg=freqs_DoppShift_ug2_neg[np.isfinite(freqs_DoppShift_ug2_neg)]

    freqs_DoppShift_radar_pos=np.array((freqs_DoppShift_radar_pos)) 
    Pxx_DoppShift_radar_pos=Pxx_DoppShift_radar_pos[np.isfinite(freqs_DoppShift_radar_pos)]
    freqs_DoppShift_radar_pos=freqs_DoppShift_radar_pos[np.isfinite(freqs_DoppShift_radar_pos)]
    
    freqs_DoppShift_radar_neg=np.array((freqs_DoppShift_radar_neg)) 
    Pxx_DoppShift_radar_neg=Pxx_DoppShift_radar_neg[np.isfinite(freqs_DoppShift_radar_neg)]
    freqs_DoppShift_radar_neg=freqs_DoppShift_radar_neg[np.isfinite(freqs_DoppShift_radar_neg)]
    
    M0_ug1, M1_ug1, M2_ug1, M3_ug1, M4_ug1, MM1_ug1, MM2_ug1 = compute_wave_spectrum_moments(freqs_ug1,Pxx_ug1,min_freq=0.05)
    Hs_ug1= np.sqrt(M0_ug1) * 4.0
    Tp_ug1=1/(freqs_ug1[np.argmax(Pxx_ug1)])
    Tme_ug1=MM1_ug1/M0_ug1
    Tf_ug1=M0_ug1/M1_ug1
    Tzc_ug1=np.sqrt(M0_ug1/M2_ug1)
    
    M0_ug2, M1_ug2, M2_ug2, M3_ug2, M4_ug2, MM1_ug2, MM2_ug2 = compute_wave_spectrum_moments(freqs_ug2,Pxx_ug2,min_freq=0.05)
    Hs_ug2= np.sqrt(M0_ug2) * 4.0
    Tp_ug2=1/(freqs_ug2[np.argmax(Pxx_ug2)])
    Tme_ug2=MM1_ug2/M0_ug2
    Tf_ug2=M0_ug2/M1_ug2
    Tzc_ug2=np.sqrt(M0_ug2/M2_ug2)
    
    M0_radar, M1_radar, M2_radar, M3_radar, M4_radar, MM1_radar, MM2_radar = compute_wave_spectrum_moments(freqs_radar,Pxx_radar,min_freq=0.05)
    Hs_radar= np.sqrt(M0_radar) * 4.0
    Tp_radar=1/(freqs_radar[np.argmax(Pxx_radar)])
    Tme_radar=MM1_radar/M0_radar
    Tf_radar=M0_radar/M1_radar
    Tzc_radar=np.sqrt(M0_radar/M2_radar)
    
    #If wave direction is nan value near the shore, the dopplershift cannot be calculated
    if freqs_DoppShift_ug1_pos.shape[0]>0:
        M0_DoppShift_ug1, M1_DoppShift_ug1, M2_DoppShift_ug1, M3_DoppShift_ug1, M4_DoppShift_ug1, MM1_DoppShift_ug1, MM2_DoppShift_ug1 = compute_wave_spectrum_moments(freqs_DoppShift_ug1_pos,Pxx_DoppShift_ug1_pos,min_freq=0.05)
        Hs_DoppShift_ug1_pos= np.sqrt(M0_DoppShift_ug1) * 4.0
        Tp_DoppShift_ug1_pos=1/(freqs_DoppShift_ug1_pos[np.argmax(Pxx_DoppShift_ug1_pos)])
        Tme_DoppShift_ug1_pos=MM1_DoppShift_ug1/M0_DoppShift_ug1   
        Tf_DoppShift_ug1_pos=M0_DoppShift_ug1/M1_DoppShift_ug1
        Tzc_DoppShift_ug1_pos=np.sqrt(M0_DoppShift_ug1/M2_DoppShift_ug1)
    
        M0_DoppShift_ug2, M1_DoppShift_ug2, M2_DoppShift_ug2, M3_DoppShift_ug2, M4_DoppShift_ug2, MM1_DoppShift_ug2, MM2_DoppShift_ug2 = compute_wave_spectrum_moments(freqs_DoppShift_ug2_pos,Pxx_DoppShift_ug2_pos,min_freq=0.05)
        Hs_DoppShift_ug2_pos= np.sqrt(M0_DoppShift_ug2) * 4.0
        Tp_DoppShift_ug2_pos=1/(freqs_DoppShift_ug2_pos[np.argmax(Pxx_DoppShift_ug2_pos)])
        Tme_DoppShift_ug2_pos=MM1_DoppShift_ug2/M0_DoppShift_ug2    
        Tf_DoppShift_ug2_pos=M0_DoppShift_ug2/M1_DoppShift_ug2
        Tzc_DoppShift_ug2_pos=np.sqrt(M0_DoppShift_ug2/M2_DoppShift_ug2)
 
        M0_DoppShift_radar, M1_DoppShift_radar, M2_DoppShift_radar, M3_DoppShift_radar, M4_DoppShift_radar, MM1_DoppShift_radar, MM2_DoppShift_radar = compute_wave_spectrum_moments(freqs_DoppShift_radar_pos,Pxx_DoppShift_radar_pos,min_freq=0.05)
        Hs_DoppShift_radar_pos= np.sqrt(M0_DoppShift_radar) * 4.0
        Tp_DoppShift_radar_pos=1/(freqs_DoppShift_radar_pos[np.argmax(Pxx_DoppShift_radar_pos)])
        Tme_DoppShift_radar_pos=MM1_DoppShift_radar/M0_DoppShift_radar    
        Tf_DoppShift_radar_pos=M0_DoppShift_radar/M1_DoppShift_radar
        Tzc_DoppShift_radar_pos=np.sqrt(M0_DoppShift_radar/M2_DoppShift_radar)
        
    if freqs_DoppShift_ug1_neg.shape[0]>0:
        M0_DoppShift_ug1, M1_DoppShift_ug1, M2_DoppShift_ug1, M3_DoppShift_ug1, M4_DoppShift_ug1, MM1_DoppShift_ug1, MM2_DoppShift_ug1 = compute_wave_spectrum_moments(freqs_DoppShift_ug1_neg,Pxx_DoppShift_ug1_neg,min_freq=0.05)
        Hs_DoppShift_ug1_neg= np.sqrt(M0_DoppShift_ug1) * 4.0
        Tp_DoppShift_ug1_neg=1/(freqs_DoppShift_ug1_neg[np.argmax(Pxx_DoppShift_ug1_neg)])
        Tme_DoppShift_ug1_neg=MM1_DoppShift_ug1/M0_DoppShift_ug1   
        Tf_DoppShift_ug1_neg=M0_DoppShift_ug1/M1_DoppShift_ug1
        Tzc_DoppShift_ug1_neg=np.sqrt(M0_DoppShift_ug1/M2_DoppShift_ug1)
    
        M0_DoppShift_ug2, M1_DoppShift_ug2, M2_DoppShift_ug2, M3_DoppShift_ug2, M4_DoppShift_ug2, MM1_DoppShift_ug2, MM2_DoppShift_ug2 = compute_wave_spectrum_moments(freqs_DoppShift_ug2_neg,Pxx_DoppShift_ug2_neg,min_freq=0.05)
        Hs_DoppShift_ug2_neg= np.sqrt(M0_DoppShift_ug2) * 4.0
        Tp_DoppShift_ug2_neg=1/(freqs_DoppShift_ug2_neg[np.argmax(Pxx_DoppShift_ug2_neg)])
        Tme_DoppShift_ug2_neg=MM1_DoppShift_ug2/M0_DoppShift_ug2    
        Tf_DoppShift_ug2_neg=M0_DoppShift_ug2/M1_DoppShift_ug2
        Tzc_DoppShift_ug2_neg=np.sqrt(M0_DoppShift_ug2/M2_DoppShift_ug2)
 
        M0_DoppShift_radar, M1_DoppShift_radar, M2_DoppShift_radar, M3_DoppShift_radar, M4_DoppShift_radar, MM1_DoppShift_radar, MM2_DoppShift_radar = compute_wave_spectrum_moments(freqs_DoppShift_radar_neg,Pxx_DoppShift_radar_neg,min_freq=0.05)
        Hs_DoppShift_radar_neg= np.sqrt(M0_DoppShift_radar) * 4.0
        Tp_DoppShift_radar_neg=1/(freqs_DoppShift_radar_neg[np.argmax(Pxx_DoppShift_radar_neg)])
        Tme_DoppShift_radar_neg=MM1_DoppShift_radar/M0_DoppShift_radar    
        Tf_DoppShift_radar_neg=M0_DoppShift_radar/M1_DoppShift_radar
        Tzc_DoppShift_radar_neg=np.sqrt(M0_DoppShift_radar/M2_DoppShift_radar)        
    
    if freqs_DoppShift_ug1_pos.shape[0]>0 or freqs_DoppShift_ug1_neg.shape[0]>0:    
        Hs_DoppShift_ug1 = max(Hs_DoppShift_ug1_pos,Hs_DoppShift_ug1_neg)
        Tp_DoppShift_ug1 = max(Tp_DoppShift_ug1_pos,Tp_DoppShift_ug1_neg)
        Tme_DoppShift_ug1 = max(Tme_DoppShift_ug1_pos,Tme_DoppShift_ug1_neg)
        Tf_DoppShift_ug1 = max(Tf_DoppShift_ug1_pos,Tf_DoppShift_ug1_neg)
        Tzc_DoppShift_ug1 = max(Tzc_DoppShift_ug1_pos,Tzc_DoppShift_ug1_neg)
        
        Hs_DoppShift_ug2 = max(Hs_DoppShift_ug2_pos,Hs_DoppShift_ug2_neg)
        Tp_DoppShift_ug2 = max(Tp_DoppShift_ug2_pos,Tp_DoppShift_ug2_neg)
        Tme_DoppShift_ug2 = max(Tme_DoppShift_ug2_pos,Tme_DoppShift_ug2_neg)
        Tf_DoppShift_ug2 = max(Tf_DoppShift_ug2_pos,Tf_DoppShift_ug2_neg)
        Tzc_DoppShift_ug2 = max(Tzc_DoppShift_ug2_pos,Tzc_DoppShift_ug2_neg)
        
        Hs_DoppShift_radar = max(Hs_DoppShift_radar_pos,Hs_DoppShift_radar_neg)
        Tp_DoppShift_radar = max(Tp_DoppShift_radar_pos,Tp_DoppShift_radar_neg)
        Tme_DoppShift_radar = max(Tme_DoppShift_radar_pos,Tme_DoppShift_radar_neg)
        Tf_DoppShift_radar = max(Tf_DoppShift_radar_pos,Tf_DoppShift_radar_neg)
        Tzc_DoppShift_radar = max(Tf_DoppShift_radar_pos,Tf_DoppShift_radar_neg)   
     
    else:
        Hs_DoppShift_ug1 = np.nan
        Tp_DoppShift_ug1 = np.nan
        Tme_DoppShift_ug1 = np.nan
        Tf_DoppShift_ug1 = np.nan
        Tzc_DoppShift_ug1 = np.nan
        
        Hs_DoppShift_ug2 = np.nan
        Tp_DoppShift_ug2 = np.nan
        Tme_DoppShift_ug2 = np.nan
        Tf_DoppShift_ug2 = np.nan
        Tzc_DoppShift_ug2 = np.nan
        
        Hs_DoppShift_radar = np.nan
        Tp_DoppShift_radar = np.nan
        Tme_DoppShift_radar = np.nan
        Tf_DoppShift_radar = np.nan
        Tzc_DoppShift_radar = np.nan
    
    return data_len,data_used,SWH_ug1,SWH_ug2,SWH_radar,Hs_ug1_manual,Hs_ug2_manual,Hs_radar_manual,Hs_ug1,Hs_DoppShift_ug1,Tp_ug1,Tp_DoppShift_ug1,Tme_ug1,Tme_DoppShift_ug1,Tf_ug1,Tf_DoppShift_ug1,Tzc_ug1,Tzc_DoppShift_ug1,Hs_ug2,Hs_DoppShift_ug2,Tp_ug2,Tp_DoppShift_ug2,Tme_ug2,Tme_DoppShift_ug2,Tf_ug2,Tf_DoppShift_ug2,Tzc_ug2,Tzc_DoppShift_ug2,Hs_radar,Hs_DoppShift_radar,Tp_radar,Tp_DoppShift_radar,Tme_radar,Tme_DoppShift_radar,Tf_radar,Tf_DoppShift_radar,Tzc_radar,Tzc_DoppShift_radar,heading,kts,wave_dir,wind_dir,wind,lon,lat #,SWH_ug1_scaled,SWH_ug2_scaled,SWH_radar_scaled

def get_wave_data_all_files():
    """Goes through folders and files in base_path, computes the parameters with get_wave_data_one_file()
    and gathers all the relevant information on wave parameters, ship's position, wind, speed and heading in one dict"""
    
    base_path = "/Users/judiththuolberg/Master/data"
    
    ship_dict=load_ship_stats()
    wave_dict=load_wave_direction()
    
    SWH_ug1=[]
    SWH_ug2=[]
    SWH_radar=[]    
    Hs_ug1_manual=[]
    Hs_ug2_manual=[]
    Hs_radar_manual=[]
    Hs_ug1=[]
    Hs_ug2=[]
    Hs_radar=[]
    Hs_DoppShift_ug1=[]
    Hs_DoppShift_ug2=[]
    Hs_DoppShift_radar=[]
    Tp_ug1=[]
    Tp_DoppShift_ug1=[]
    Tme_ug1=[]
    Tme_DoppShift_ug1=[]
    Tf_ug1=[]
    Tf_DoppShift_ug1=[]
    Tzc_ug1=[]
    Tzc_DoppShift_ug1=[]
    Tp_ug2=[]
    Tp_DoppShift_ug2=[]
    Tme_ug2=[]
    Tme_DoppShift_ug2=[]
    Tf_ug2=[]
    Tf_DoppShift_ug2=[]
    Tzc_ug2=[]
    Tzc_DoppShift_ug2=[]
    Tp_radar=[]
    Tp_DoppShift_radar=[]
    Tme_radar=[]
    Tme_DoppShift_radar=[]
    Tf_radar=[]
    Tf_DoppShift_radar=[]
    Tzc_radar=[]
    Tzc_DoppShift_radar=[]    
    date_time=[]
    heading=[]
    kts=[]
    wave_dir=[]
    wind_dir=[]
    wind=[]
    lon=[]
    lat=[]
    data_len = []
    data_used =[]

    dir_cnt=0
    for (root,dirs,files) in os.walk(base_path):
        dirs.sort()
        files.sort()
        dir_cnt+=1
        file_cnt=0
    
        for file_name in files:
            if not(file_name.endswith(".lzma")):
                continue
            file_cnt+=1
            print('\x1b[1K' ,end="\r")
            print("folder counter: ",dir_cnt," file counter: ",file_cnt, end="\r")
            file_path=root+"/"+file_name
            dict_data=load_data_dump(file_path)

            curr_data_len,curr_data_used,curr_SWH_ug1,curr_SWH_ug2,curr_SWH_radar,curr_Hs_ug1_manual,curr_Hs_ug2_manual,curr_Hs_radar_manual,curr_Hs_ug1,curr_Hs_DoppShift_ug1,curr_Tp_ug1,curr_Tp_DoppShift_ug1,curr_Tme_ug1,curr_Tme_DoppShift_ug1,curr_Tf_ug1,curr_Tf_DoppShift_ug1,curr_Tzc_ug1,curr_Tzc_DoppShift_ug1,curr_Hs_ug2,curr_Hs_DoppShift_ug2,curr_Tp_ug2,curr_Tp_DoppShift_ug2,curr_Tme_ug2,curr_Tme_DoppShift_ug2,curr_Tf_ug2,curr_Tf_DoppShift_ug2,curr_Tzc_ug2,curr_Tzc_DoppShift_ug2,curr_Hs_radar,curr_Hs_DoppShift_radar,curr_Tp_radar,curr_Tp_DoppShift_radar,curr_Tme_radar,curr_Tme_DoppShift_radar,curr_Tf_radar,curr_Tf_DoppShift_radar,curr_Tzc_radar,curr_Tzc_DoppShift_radar,curr_heading,curr_kts,curr_wave_dir,curr_wind_dir,curr_wind,curr_lon,curr_lat = get_wave_data_one_file(dict_data,ship_dict,wave_dict) #curr_SWH_ug1_scaled,curr_SWH_ug2_scaled,curr_SWH_radar_scaled
                
            SWH_ug1.append(curr_SWH_ug1)
            SWH_ug2.append(curr_SWH_ug2)
            SWH_radar.append(curr_SWH_radar)
            Hs_ug1_manual.append(curr_Hs_ug1_manual)
            Hs_ug2_manual.append(curr_Hs_ug2_manual)
            Hs_radar_manual.append(curr_Hs_radar_manual)
            
            Hs_ug1.append(curr_Hs_ug1)
            Hs_DoppShift_ug1.append(curr_Hs_DoppShift_ug1)
            Tp_ug1.append(curr_Tp_ug1)
            Tp_DoppShift_ug1.append(curr_Tp_DoppShift_ug1)
            Tme_ug1.append(curr_Tme_ug1)
            Tme_DoppShift_ug1.append(curr_Tme_DoppShift_ug1)
            Tf_ug1.append(curr_Tf_ug1)
            Tf_DoppShift_ug1.append(curr_Tf_DoppShift_ug1)
            Tzc_ug1.append(curr_Tzc_ug1)
            Tzc_DoppShift_ug1.append(curr_Tzc_DoppShift_ug1)
            
            Hs_ug2.append(curr_Hs_ug2)
            Hs_DoppShift_ug2.append(curr_Hs_DoppShift_ug2)
            Tp_ug2.append(curr_Tp_ug2)
            Tp_DoppShift_ug2.append(curr_Tp_DoppShift_ug2)
            Tme_ug2.append(curr_Tme_ug2)
            Tme_DoppShift_ug2.append(curr_Tme_DoppShift_ug2)
            Tf_ug2.append(curr_Tf_ug2)
            Tf_DoppShift_ug2.append(curr_Tf_DoppShift_ug2)
            Tzc_ug2.append(curr_Tzc_ug2)
            Tzc_DoppShift_ug2.append(curr_Tzc_DoppShift_ug2)
            
            Hs_radar.append(curr_Hs_radar)
            Hs_DoppShift_radar.append(curr_Hs_DoppShift_radar)
            Tp_radar.append(curr_Tp_radar)
            Tp_DoppShift_radar.append(curr_Tp_DoppShift_radar)
            Tme_radar.append(curr_Tme_radar)
            Tme_DoppShift_radar.append(curr_Tme_DoppShift_radar)
            Tf_radar.append(curr_Tf_radar)
            Tf_DoppShift_radar.append(curr_Tf_DoppShift_radar)
            Tzc_radar.append(curr_Tzc_radar)
            Tzc_DoppShift_radar.append(curr_Tzc_DoppShift_radar)

            heading.append(curr_heading)
            kts.append(curr_kts)
            wave_dir.append(curr_wave_dir)
            wind_dir.append(curr_wind_dir)
            wind.append(curr_wind)
            lon.append(curr_lon)
            lat.append(curr_lat)
            date_time.append(dict_data["list_time_series_interpolated"]["common_datetime"][0])
            data_len.append(curr_data_len)
            data_used.append(curr_data_used)
    
    SWH_ug1=np.array((SWH_ug1))
    SWH_ug2=np.array((SWH_ug2))
    SWH_radar=np.array((SWH_radar))
    Hs_ug1_manual=np.array((Hs_ug1_manual))
    Hs_ug2_manual=np.array((Hs_ug2_manual))
    Hs_radar_manual=np.array((Hs_radar_manual))
    
    Hs_ug1=np.array((Hs_ug1))
    Hs_DoppShift_ug1=np.array((Hs_DoppShift_ug1))
    Tp_ug1=np.array((Tp_ug1))
    Tp_DoppShift_ug1=np.array((Tp_DoppShift_ug1))
    Tme_ug1=np.array((Tme_ug1))
    Tme_DoppShift_ug1=np.array((Tme_DoppShift_ug1))
    Tf_ug1=np.array((Tf_ug1))
    Tf_DoppShift_ug1=np.array((Tf_DoppShift_ug1))
    Tzc_ug1=np.array((Tzc_ug1))
    Tzc_DoppShift_ug1=np.array((Tzc_DoppShift_ug1))
    
    Hs_ug2=np.array((Hs_ug2))
    Hs_DoppShift_ug2=np.array((Hs_DoppShift_ug2))
    Tp_ug2=np.array((Tp_ug2))
    Tp_DoppShift_ug2=np.array((Tp_DoppShift_ug2))
    Tme_ug2=np.array((Tme_ug2))
    Tme_DoppShift_ug2=np.array((Tme_DoppShift_ug2))
    Tf_ug2=np.array((Tf_ug2))
    Tf_DoppShift_ug2=np.array((Tf_DoppShift_ug2))
    Tzc_ug2=np.array((Tzc_ug2))
    Tzc_DoppShift_ug2=np.array((Tzc_DoppShift_ug2))
    
    Hs_radar=np.array((Hs_radar))
    Hs_DoppShift_radar=np.array((Hs_DoppShift_radar))
    Tp_radar=np.array((Tp_radar))
    Tp_DoppShift_radar=np.array((Tp_DoppShift_radar))
    Tme_radar=np.array((Tme_radar))
    Tme_DoppShift_radar=np.array((Tme_DoppShift_radar))
    Tf_radar=np.array((Tf_radar))
    Tf_DoppShift_radar=np.array((Tf_DoppShift_radar))
    Tzc_radar=np.array((Tzc_radar))
    Tzc_DoppShift_radar=np.array((Tzc_DoppShift_radar))

    heading=np.array((heading))
    kts=np.array((kts))
    wave_dir=np.array((wave_dir))
    wind_dir=np.array((wind_dir))
    wind=np.array((wind))
    lon=np.array((lon))
    lat=np.array((lat))    
    date_time=np.array((date_time))
    data_len=np.array(data_len)
    data_used=np.array(data_used)
    
    ret_dict={}
    ret_dict["datetime"]=date_time
    ret_dict["heading"]=heading
    ret_dict["wave_dir"]=wave_dir
    ret_dict["wind_dir"]=wind_dir
    ret_dict["kts"]=kts
    ret_dict["wind"]=wind
    ret_dict["lon"]=lon
    ret_dict["lat"]=lat
    ret_dict['data_len']=data_len
    ret_dict['data_used']=data_used
    ret_dict["ug1"]={}
    ret_dict["ug2"]={}
    ret_dict["radar"]={}
    
    ret_dict["ug1"]["Hs"]={}
    ret_dict["ug1"]["Hs"]["SWH"]=SWH_ug1
    ret_dict["ug1"]["Hs"]["H_1/3"]=Hs_ug1_manual
    ret_dict["ug1"]["Hs"]["Hs"]=Hs_ug1
    ret_dict["ug1"]["Hs"]["Hs_DoppShift"]=Hs_DoppShift_ug1
    ret_dict["ug1"]["Tp"]={}
    ret_dict["ug1"]["Tp"]["Tp"]=Tp_ug1
    ret_dict["ug1"]["Tp"]["Tp_DoppShift"]=Tp_DoppShift_ug1
    ret_dict["ug1"]["Tme"]={}
    ret_dict["ug1"]["Tme"]["Tme"]=Tme_ug1
    ret_dict["ug1"]["Tme"]["Tme_DoppShift"]=Tme_DoppShift_ug1
    ret_dict["ug1"]["Tf"]={}
    ret_dict["ug1"]["Tf"]["Tf"]=Tf_ug1
    ret_dict["ug1"]["Tf"]["Tf_DoppShift"]=Tf_DoppShift_ug1
    ret_dict["ug1"]["Tzc"]={}
    ret_dict["ug1"]["Tzc"]["Tzc"]=Tzc_ug1
    ret_dict["ug1"]["Tzc"]["Tzc_DoppShift"]=Tzc_DoppShift_ug1
    
    ret_dict["ug2"]["Hs"]={}
    ret_dict["ug2"]["Hs"]["SWH"]=SWH_ug2
    ret_dict["ug2"]["Hs"]["H_1/3"]=Hs_ug2_manual
    ret_dict["ug2"]["Hs"]["Hs"]=Hs_ug2
    ret_dict["ug2"]["Hs"]["Hs_DoppShift"]=Hs_DoppShift_ug2
    ret_dict["ug2"]["Tp"]={}    
    ret_dict["ug2"]["Tp"]["Tp"]=Tp_ug2
    ret_dict["ug2"]["Tp"]["Tp_DoppShift"]=Tp_DoppShift_ug2
    ret_dict["ug2"]["Tme"]={}
    ret_dict["ug2"]["Tme"]["Tme"]=Tme_ug2
    ret_dict["ug2"]["Tme"]["Tme_DoppShift"]=Tme_DoppShift_ug2
    ret_dict["ug2"]["Tf"]={}
    ret_dict["ug2"]["Tf"]["Tf"]=Tf_ug2
    ret_dict["ug2"]["Tf"]["Tf_DoppShift"]=Tf_DoppShift_ug2
    ret_dict["ug2"]["Tzc"]={}
    ret_dict["ug2"]["Tzc"]["Tzc"]=Tzc_ug2
    ret_dict["ug2"]["Tzc"]["Tzc_DoppShift"]=Tzc_DoppShift_ug2
    
    ret_dict["radar"]["Hs"]={}
    ret_dict["radar"]["Hs"]["SWH"]=SWH_radar
    ret_dict["radar"]["Hs"]["H_1/3"]=Hs_radar_manual
    ret_dict["radar"]["Hs"]["Hs"]=Hs_radar
    ret_dict["radar"]["Hs"]["Hs_DoppShift"]=Hs_DoppShift_radar
    ret_dict["radar"]["Tp"]={}    
    ret_dict["radar"]["Tp"]["Tp"]=Tp_radar
    ret_dict["radar"]["Tp"]["Tp_DoppShift"]=Tp_DoppShift_radar
    ret_dict["radar"]["Tme"]={}
    ret_dict["radar"]["Tme"]["Tme"]=Tme_radar
    ret_dict["radar"]["Tme"]["Tme_DoppShift"]=Tme_DoppShift_radar
    ret_dict["radar"]["Tf"]={}
    ret_dict["radar"]["Tf"]["Tf"]=Tf_radar
    ret_dict["radar"]["Tf"]["Tf_DoppShift"]=Tf_DoppShift_radar
    ret_dict["radar"]["Tzc"]={}
    ret_dict["radar"]["Tzc"]["Tzc"]=Tzc_radar
    ret_dict["radar"]["Tzc"]["Tzc_DoppShift"]=Tzc_DoppShift_radar

    return ret_dict


def PSD_over_some_hours(probe_name,starttime=10,min_per_PSD=30,number_of_files=6):
    #plots the PSDs over some hours for some files in the same folder
    base_path = "/Users/judiththuolberg/master/data/2021/11/20"
    
    ship_dict=load_ship_stats()
    wave_dict=load_wave_direction()
    
    pitch=[]
    roll=[]
    accel_D=[]
    ug1=[]
    ug2=[]
    radar=[]
    date_time=[]
     
    dir_cnt=0
    for (root,dirs,files) in os.walk(base_path):
        dirs.sort()
        files.sort()
        dir_cnt+=1
        file_cnt=0
    
        for file_name in files[starttime*2:]:
            if not(file_name.endswith(".lzma")):
                continue
            file_cnt+=1
            print('\x1b[1K' ,end="\r")
            print("folder counter: ",dir_cnt," file counter: ",file_cnt, end="\r")
            file_path=root+"/"+file_name
            dict_data=load_data_dump(file_path)
    
            pitch.extend(dict_data["list_time_series_interpolated"]["IMU"]["pitch"])
            roll.extend(dict_data["list_time_series_interpolated"]["IMU"]["roll"])
            accel_D.extend(dict_data["list_time_series_interpolated"]["IMU"]["accel_D"])
            ug1.extend(dict_data["list_time_series_interpolated"]["gauges"]["gauge_ug1_meters"])
            ug2.extend(dict_data["list_time_series_interpolated"]["gauges"]["gauge_ug1_meters"])
            radar.extend(dict_data["list_time_series_interpolated"]["gauges"]["gauge_radar1_meters"])
            date_time.extend(dict_data["list_time_series_interpolated"]["common_datetime"])
            
            if file_cnt >= number_of_files:
                break
        if file_cnt >= number_of_files:
            break
    
    pitch=np.array((pitch))
    roll=np.array((roll))
    accel_D=np.array((accel_D))
    date_time=np.array((date_time))
    pitch, roll = corr_pitch_n_roll_array("IMU",pitch,roll,date_time)
    pitch=pitch/360*(2*np.pi)
    roll=roll/360*(2*np.pi)
      
    #the vertical positions of the VN100 by two times integration
    sample_freq=10
    accel_D_mean=np.mean(accel_D)
    z_tt_VN100=-(accel_D-accel_D_mean)
    z_VN100=fft_int_2(z_tt_VN100,sample_freq)
    
    #Due to Fourier transformations the data is cut off in beginning and end
    data_len=pitch.shape[0]
    cutoff=data_len//10
    date_time=date_time[cutoff:-cutoff]
    pitch=pitch[cutoff:-cutoff]
    roll=roll[cutoff:-cutoff]
    z_VN100=z_VN100[cutoff:-cutoff]
    ug1=ug1[cutoff:-cutoff]
    ug2=ug2[cutoff:-cutoff]
    radar=radar[cutoff:-cutoff]                    
    
    #the vertical positions of the sensors are computed
    z_corr_ug1=simple_trig_approx_ug1(pitch,roll,z_VN100)
    z_ug1=(z_VN100+z_corr_ug1)-np.mean(z_VN100+z_corr_ug1)
    z_corr_ug2=simple_trig_approx_ug2(pitch,roll,z_VN100)
    z_ug2=(z_VN100+z_corr_ug2)-np.mean(z_VN100+z_corr_ug2)
    pos_z_ug1=z_ug1
    pos_z_ug2=z_ug2
    
    #the water suface elevation WH_ossc is computed
    WH_ug1=pos_z_ug1-ug1*np.cos(pitch)*np.cos(roll)
    WH_mean_ug1=np.mean(WH_ug1)    
    WH_ossc_ug1=WH_ug1-WH_mean_ug1
    
    WH_ug2=pos_z_ug2-ug2*np.cos(pitch)*np.cos(roll)
    WH_mean_ug2=np.mean(WH_ug2)    
    WH_ossc_ug2=WH_ug2-WH_mean_ug2
    
    WH_radar=pos_z_ug2-radar*np.cos(pitch)*np.cos(roll)
    WH_mean_radar=np.mean(WH_radar)    
    WH_ossc_radar=WH_radar-WH_mean_radar

    #calculating one value within the 30minutes datasetwith minute resolution
    t_start=date_time[0]
    t_end=date_time[-1]
    kts_array=ship_dict["kts"][np.where(((ship_dict["datetime"]>t_start)&(ship_dict["datetime"]<t_end)))]
    heading_array=ship_dict["heading"][np.where(((ship_dict["datetime"]>t_start)&(ship_dict["datetime"]<t_end)))]

    kts=np.mean(kts_array)
    heading=np.median(heading_array)   #need to calculate this better for values betwen 0 and 360
    
    if heading == 0:
        heading+=0.1
    
    if kts == 0:
        kts += 0.1
    
    #wave direction
    delta_t_start=np.abs(wave_dict["datetime"]-t_start)
    closest_t_start=np.min(delta_t_start)
    delta_t_end=np.abs(wave_dict["datetime"]-t_end)
    closest_t_end=np.min(delta_t_end)
    if closest_t_start<=closest_t_end:
        ind_closest_t=np.argmin(delta_t_start)
        wave_dir=wave_dict["mdir"][ind_closest_t]
    else:
        ind_closest_t=np.argmin(delta_t_end)
        wave_dir=wave_dict["mdir"][ind_closest_t]
        
    if wave_dir == 0:
        wave_dir +=0.1    
    
    seg_secs=400#180 is smoother than 400
    
    starttimes=[]
    SWHs=[]
    power_law=[]
    PSD_origs=[]
    freqs_origs=[]
    PSD_Dopplers=[]
    freqs_Dopplers=[]
    Tp_Dopplers=[]
    
    if probe_name=='ug1':
        WH_ossc = WH_ossc_ug1
    elif probe_name=='ug2':
        WH_ossc = WH_ossc_ug2        
    elif probe_name=='radar':
        WH_ossc = WH_ossc_radar
    
    data_len=(t_end-t_start).seconds
    num_PSDs=data_len//(min_per_PSD*60)
    for cnt in range(num_PSDs):
        crr_WH_ossc=WH_ossc[cnt*min_per_PSD*600:(cnt+1)*min_per_PSD*600]
        starttime=date_time[cnt*min_per_PSD*600]
        starttimes.append(starttime)
        
        SWH=4*np.std(crr_WH_ossc)
        
        freqs,Pxx=welch_spectrum(crr_WH_ossc,segment_duration_seconds=seg_secs)
        Pxx_orig=Pxx
        freqs_orig=freqs
        
        Pxx_DoppShift_pos=Pxx_orig*1
        Pxx_DoppShift_neg=Pxx_orig*1
        freqs_DoppShift_pos=[]
        freqs_DoppShift_neg=[]
        for freq in freqs_orig:
            freq_corr_pos,freq_corr_neg=correct_freq_for_Doppler_effect(freq,kts,heading,wave_dir)
            freqs_DoppShift_pos.append(freq_corr_pos)
            freqs_DoppShift_neg.append(freq_corr_neg)
            
        #keep only non-nan values    
        freqs_DoppShift_pos=np.array((freqs_DoppShift_pos))
        Pxx_DoppShift_pos=Pxx_DoppShift_pos[np.isfinite(freqs_DoppShift_pos)]
        freqs_DoppShift_pos=freqs_DoppShift_pos[np.isfinite(freqs_DoppShift_pos)]
        
        freqs_DoppShift_neg=np.array((freqs_DoppShift_neg))
        Pxx_DoppShift_neg=Pxx_DoppShift_neg[np.isfinite(freqs_DoppShift_neg)]
        freqs_DoppShift_neg=freqs_DoppShift_neg[np.isfinite(freqs_DoppShift_neg)]
        
        Tp_orig=1/(freqs_orig[np.argmax(Pxx_orig)])
        Tp_DoppShift_pos=1/(freqs_DoppShift_pos[np.argmax(Pxx_DoppShift_pos)])
        Tp_DoppShift_neg=1/(freqs_DoppShift_neg[np.argmax(Pxx_DoppShift_neg)])
        
        #Test the fit of the power law freqs_orig**(slope)*10**(y-axis interception at lower right corner)
        tail=freqs_orig**(-2.3)*10**(-2.5) 
        
        if Tp_DoppShift_pos<5:
            print('check doppler shifted frequency value')
        
        SWHs.append(SWH)
        PSD_origs.append(Pxx_orig)
        freqs_origs.append(freqs_orig)
        Tp_Dopplers.append(Tp_DoppShift_pos)
        freqs_Dopplers.append(freqs_DoppShift_pos)
        PSD_Dopplers.append(Pxx_DoppShift_pos)
        power_law.append(tail) #power law for the tail of the power spectrum

    fig,ax=plt.subplots(3,2,figsize=(17,21))
    #original frequency and spectum
    for cnt in range(len(starttimes)):
        lbl=(starttimes[cnt]).strftime("%H:%M") #"PSD uncorrected "+str(min_per_PSD)+"min from "+
        ax[0,0].plot(freqs_origs[cnt],PSD_origs[cnt],label=lbl)
    ax[0,0].legend()
    ax[0,0].set_xlim((0.05,0.6))
    ax[0,0].set_ylabel("PSD uncorrected")
    ax[0,0].set_xlabel("Frequency [Hz]")
    plot0_title="Uncorrected PSD for "+starttimes[0].strftime("%d.%m.%Y")# %H:%M:%S")+", "+"Segment duration = "+str(seg_secs)+" seconds"
    ax[0,0].set_title(plot0_title)
    
    #doppler shifted frequency and spectrum
    for cnt in range(len(starttimes)):
        lbl=(starttimes[cnt]).strftime("%H:%M") #"PSD Doppler effect corrected "+str(min_per_PSD)+"min from "+
        ax[0,1].plot(freqs_Dopplers[cnt],PSD_Dopplers[cnt],label=lbl)
    ax[0,1].legend()
    ax[0,1].set_xlim((0.05,0.6))
    ax[0,1].set_ylabel("PSD Doppler effect corrected")
    ax[0,1].set_xlabel("Frequency [Hz]")
    plot1_title="Doppler effect corrected PSD, "+"Segment duration = "+str(seg_secs)+" seconds"
    ax[0,1].set_title(plot1_title)
   
    #original frequency and spectum logaritmic   
    for cnt in range(len(starttimes)):
        lbl=(starttimes[cnt]).strftime("%H:%M") #"PSD uncorrected "+str(min_per_PSD)+"min from "+
        ax[1,0].loglog(freqs_origs[cnt],PSD_origs[cnt],label=lbl)
    ax[1,0].loglog(freqs_origs[0],power_law[0],color='k')
    ax[1,0].legend()
    ax[1,0].set_xlim((0.05,0.6))
    ax[1,0].set_ylim((10**(-3),10**(0.5)))
    ax[1,0].set_ylabel("Log (PSD uncorrected)")
    ax[1,0].set_xlabel("Log (Frequency [Hz])")
    plot2_title="Uncorrected PSD for "+starttimes[0].strftime("%d.%m.%Y ")#%H:%M:%S")+", "+"Segment duration = "+str(seg_secs)+" seconds"
    ax[1,0].set_title(plot2_title)
    
    #doppler shifted frequency and spectrum logaritmic
    for cnt in range(len(starttimes)):
        lbl=(starttimes[cnt]).strftime("%H:%M") #"PSD Doppler effect corrected "+str(min_per_PSD)+"min from "+
        ax[1,1].loglog(freqs_Dopplers[cnt],PSD_Dopplers[cnt],label=lbl)
    ax[1,1].loglog(freqs_Dopplers[0],power_law[0],color='k')
    ax[1,1].legend()
    ax[1,1].set_xlim((0.05,0.6))
    ax[1,1].set_ylim((10**(-3),10**(0.5)))
    ax[1,1].set_ylabel("Log (PSD Doppler effect corrected)")
    ax[1,1].set_xlabel("Log (Frequency [Hz])")
    plot3_title="Doppler effect corrected PSD, "+"Segment duration = "+str(seg_secs)+" seconds"
    ax[1,1].set_title(plot3_title)
    
    plottimes=[]
    for time in starttimes:
        plottimes.append(time.strftime("%H:%M"))
        
    ax[2,0].plot(starttimes,SWHs,"b*",label="SWH")
    #ax[2,0].legend()
    ax[2,0].set_ylabel("SWH [m]")
    #ax[2,0].set_xlabel("time")
    #ax[2,0].set_ylim((0,6))
    plot4_title="Significant wave height over time"#" for 2.5 h from: "+ t_start.strftime("%d.%m.%Y %H:%M:%S")
    ax[2,0].set_title(plot4_title)
    
    ax[2,1].plot(starttimes,Tp_Dopplers,"r*",label="Tp Doppler effect corrected pos")
    ax[2,1].set_ylabel("Tp [s]")
    #ax[2,1].legend()
    #ax[1,1].set_ylim((0,15))
    #ax[2,1].set_xlabel("time")
    plot5_title="Peak period over time"#" for 2.5 h from: "+ t_start.strftime("%d.%m.%Y %H:%M:%S")
    ax[2,1].set_title(plot5_title)
