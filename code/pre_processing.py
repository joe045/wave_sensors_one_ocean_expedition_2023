#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 08:35:21 2023

@author: judiththuolberg
"""
import os
import math
import numpy as np
from sklearn.decomposition import PCA
from scipy import interpolate
from scipy.optimize import curve_fit

from wave_processing import *
#script that loads the raw data (needs arduino_logging.py) and does wave processing
from run_you_fools import *

def get_def_pitch_and_roll():
    """Goes trough all files in the folder for harbour days and returns the default value of pitch and roll for each IMU"""
    base_path="/Users/judiththuolberg/master/old_harbour/"
    
    IMU_pitch=[]
    IMU_roll=[]
    extraIMU_0_pitch=[]
    extraIMU_0_roll=[]
    extraIMU_1_pitch=[]
    extraIMU_1_roll=[]
    
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
            print("folder counter: ",dir_cnt,"file counter: ",file_cnt, end="\r")
            file_path=root+"/"+file_name
            dict_data=load_data_dump(file_path)
                
            IMU_pitch.extend(dict_data["list_time_series_interpolated"]["IMU"]["pitch"])
            IMU_roll.extend(dict_data["list_time_series_interpolated"]["IMU"]["roll"])
            extraIMU_0_pitch.extend(dict_data["list_time_series_interpolated"]["extraIMU_0"]["pitch"])
            extraIMU_0_roll.extend(dict_data["list_time_series_interpolated"]["extraIMU_0"]["roll"])
            extraIMU_1_pitch.extend(dict_data["list_time_series_interpolated"]["extraIMU_1"]["pitch"])
            extraIMU_1_roll.extend(dict_data["list_time_series_interpolated"]["extraIMU_1"]["roll"])
                
    IMU_def_pitch=np.nanmean(IMU_pitch)
    IMU_def_roll=np.nanmean(IMU_roll)
    extraIMU_0_def_pitch=np.nanmean(extraIMU_0_pitch)
    extraIMU_0_def_roll=np.nanmean(extraIMU_0_roll)
    extraIMU_1_def_pitch=np.nanmean(extraIMU_1_pitch)
    extraIMU_1_def_roll=np.nanmean(extraIMU_1_roll)
    
    return (IMU_def_pitch,IMU_def_roll,extraIMU_0_def_pitch,extraIMU_0_def_roll,extraIMU_1_def_pitch,extraIMU_1_def_roll)

def two_PCs_daily_scatter():
    """Goes trough all the files in the folder for sailing days and returns the principal components for pitch and roll for each IMU"""
    base_path="/Users/judiththuolberg/master/new_data/"
    
    IMU_PC1s_roll=[]
    IMU_PC1s_pitch=[]
    IMU_PC2s_roll=[]
    IMU_PC2s_pitch=[]
    
    extraIMU_0_PC1s_roll=[]
    extraIMU_0_PC1s_pitch=[]
    extraIMU_0_PC2s_roll=[]
    extraIMU_0_PC2s_pitch=[]
    
    extraIMU_1_PC1s_roll=[]
    extraIMU_1_PC1s_pitch=[]
    extraIMU_1_PC2s_roll=[]
    extraIMU_1_PC2s_pitch=[]
    
    dir_cnt=0
    for (root,dirs,files) in os.walk(base_path):
        dirs.sort()
        files.sort()
        dir_cnt+=1
        file_cnt=0

        pitch_IMU_day=[]
        roll_IMU_day=[]
        pitch_extraIMU_0_day=[]
        roll_extraIMU_0_day=[]
        pitch_extraIMU_1_day=[]
        roll_extraIMU_1_day=[]
        
        for file_name in files:
            if not(file_name.endswith(".lzma")):
                continue
            file_cnt+=1
            print('\x1b[1K' ,end="\r")
            print("folder counter: ",dir_cnt," file counter: ",file_cnt, end="\r")
            file_path=root+"/"+file_name
            dict_data=load_data_dump(file_path)
            
            pitch_IMU_day.extend(dict_data["list_time_series_interpolated"]["IMU"]["pitch"])
            roll_IMU_day.extend(dict_data["list_time_series_interpolated"]["IMU"]["roll"])
            pitch_extraIMU_0_day.extend(dict_data["list_time_series_interpolated"]["extraIMU_0"]["pitch"])
            roll_extraIMU_0_day.extend(dict_data["list_time_series_interpolated"]["extraIMU_0"]["roll"])
            pitch_extraIMU_1_day.extend(dict_data["list_time_series_interpolated"]["extraIMU_1"]["pitch"])
            roll_extraIMU_1_day.extend(dict_data["list_time_series_interpolated"]["extraIMU_1"]["roll"])
        
        pitch_IMU_day = [item for item in pitch_IMU_day if not(math.isnan(item)) == True]
        roll_IMU_day = [item for item in roll_IMU_day if not(math.isnan(item)) == True]
        pitch_extraIMU_0_day =[item for item in pitch_extraIMU_0_day if not(math.isnan(item)) == True]
        roll_extraIMU_0_day = [item for item in roll_extraIMU_0_day if not(math.isnan(item)) == True]
        pitch_extraIMU_1_day = [item for item in pitch_extraIMU_1_day if not(math.isnan(item)) == True]
        roll_extraIMU_1_day = [item for item in roll_extraIMU_1_day if not(math.isnan(item)) == True]
        
        if file_cnt < 1:
            continue

        if np.shape(roll_IMU_day)[0]>0 and np.shape(pitch_IMU_day)[0]>0:
            pr_IMU_array=np.array(((roll_IMU_day),(pitch_IMU_day)))
            pr_IMU_array_t=np.transpose(pr_IMU_array)
            pca_IMU=PCA().fit(pr_IMU_array_t)
            
            if pca_IMU.components_[0,0]>0:
                IMU_PC1s_roll.append(pca_IMU.components_[0,0])
                IMU_PC1s_pitch.append(pca_IMU.components_[0,1])
            else:
                IMU_PC1s_roll.append(-pca_IMU.components_[0,0])
                IMU_PC1s_pitch.append(-pca_IMU.components_[0,1])
            
            if pca_IMU.components_[1,1]>0:
                IMU_PC2s_roll.append(pca_IMU.components_[1,0])
                IMU_PC2s_pitch.append(pca_IMU.components_[1,1])
            else:
                IMU_PC2s_roll.append(-pca_IMU.components_[1,0])
                IMU_PC2s_pitch.append(-pca_IMU.components_[1,1])
                  
        if np.shape(roll_extraIMU_0_day)[0]>0 and np.shape(pitch_extraIMU_0_day)[0]>0:
            pr_extraIMU_0_array=np.array(((roll_extraIMU_0_day),(pitch_extraIMU_0_day)))
            pr_extraIMU_0_array_t=np.transpose(pr_extraIMU_0_array)
            pca_extraIMU_0=PCA().fit(pr_extraIMU_0_array_t)
            
            if pca_extraIMU_0.components_[0,0]>0:
                extraIMU_0_PC1s_roll.append(pca_extraIMU_0.components_[0,0])
                extraIMU_0_PC1s_pitch.append(pca_extraIMU_0.components_[0,1])
            else:
                extraIMU_0_PC1s_roll.append(-pca_extraIMU_0.components_[0,0])
                extraIMU_0_PC1s_pitch.append(-pca_extraIMU_0.components_[0,1])
            
            if pca_extraIMU_0.components_[1,1]>0:
                extraIMU_0_PC2s_roll.append(pca_extraIMU_0.components_[1,0])
                extraIMU_0_PC2s_pitch.append(pca_extraIMU_0.components_[1,1])
            else:
                extraIMU_0_PC2s_roll.append(-pca_extraIMU_0.components_[1,0])
                extraIMU_0_PC2s_pitch.append(-pca_extraIMU_0.components_[1,1])   
        
        if np.shape(roll_extraIMU_1_day)[0]>0 and np.shape(pitch_extraIMU_1_day)[0]>0:
            pr_extraIMU_1_array=np.array(((roll_extraIMU_1_day),(pitch_extraIMU_1_day)))
            pr_extraIMU_1_array_t=np.transpose(pr_extraIMU_1_array)
            pca_extraIMU_1=PCA().fit(pr_extraIMU_1_array_t)
        
            if pca_extraIMU_1.components_[0,0]>0:
                extraIMU_1_PC1s_roll.append(pca_extraIMU_1.components_[0,0])
                extraIMU_1_PC1s_pitch.append(pca_extraIMU_1.components_[0,1])
            else:
                extraIMU_1_PC1s_roll.append(-pca_extraIMU_1.components_[0,0])
                extraIMU_1_PC1s_pitch.append(-pca_extraIMU_1.components_[0,1]) 
        
            if pca_extraIMU_1.components_[1,1] > 0:
                extraIMU_1_PC2s_roll.append(pca_extraIMU_1.components_[1,0])
                extraIMU_1_PC2s_pitch.append(pca_extraIMU_1.components_[1,1])
            else:
                extraIMU_1_PC2s_roll.append(-pca_extraIMU_1.components_[1,0])
                extraIMU_1_PC2s_pitch.append(-pca_extraIMU_1.components_[1,1])
        
    fig,ax=plt.subplots(1,3,figsize=(15,5))
    
    ax[0].plot(IMU_PC1s_roll,IMU_PC1s_pitch,"b+",label="PC1")
    ax[0].plot(IMU_PC2s_roll,IMU_PC2s_pitch,"r+",label="PC2")
    ax[0].plot([0],[0],"k+",label="(0,0)")
    ax[0].set_xlabel("roll")
    ax[0].set_ylabel("pitch")
    ax[0].legend()
    ax[0].set_xlim((-1.1,1.1))
    ax[0].set_ylim((-1.1,1.1))
    ax[0].set_title("IMU PCs") #from 22.10. to 27.11.")
    
    ax[1].plot(extraIMU_0_PC1s_roll,extraIMU_0_PC1s_pitch,"b+",label="PC1")
    ax[1].plot(extraIMU_0_PC2s_roll,extraIMU_0_PC2s_pitch,"r+",label="PC2")
    ax[1].plot([0],[0],"k+",label="(0,0)")
    ax[1].set_xlabel("roll")
    ax[1].set_ylabel("pitch")
    ax[1].legend()
    ax[1].set_xlim((-1.1,1.1))
    ax[1].set_ylim((-1.1,1.1))
    ax[1].set_title("extraIMU_0 PCs")# from 22.10. to 27.11.")
    
    ax[2].plot(extraIMU_1_PC1s_roll,extraIMU_1_PC1s_pitch,"b+",label="PC1")
    ax[2].plot(extraIMU_1_PC2s_roll,extraIMU_1_PC2s_pitch,"r+",label="PC2")
    ax[2].plot([0],[0],"k+",label="(0,0)")
    ax[2].set_xlabel("roll")
    ax[2].set_ylabel("pitch")
    ax[2].legend()
    ax[2].set_xlim((-1.1,1.1))
    ax[2].set_ylim((-1.1,1.1))
    ax[2].set_title("extraIMU_1 PCs")#" from 22.10. to 27.11.")
    
    PC_dict={}
    PC_dict["IMU"]={}
    PC_dict["extraIMU_0"]={}
    PC_dict["extraIMU_1"]={}
    PC_dict["IMU"]["PC1s"]={}
    PC_dict["IMU"]["PC2s"]={}
    PC_dict["extraIMU_0"]["PC1s"]={}
    PC_dict["extraIMU_0"]["PC2s"]={}
    PC_dict["extraIMU_1"]["PC1s"]={}
    PC_dict["extraIMU_1"]["PC2s"]={}
    
    PC_dict["IMU"]["PC1s"]["pitch"]=IMU_PC1s_pitch
    PC_dict["IMU"]["PC1s"]["roll"]=IMU_PC1s_roll
    PC_dict["IMU"]["PC2s"]["pitch"]=IMU_PC2s_pitch
    PC_dict["IMU"]["PC2s"]["roll"]=IMU_PC2s_roll

    PC_dict["extraIMU_0"]["PC1s"]["pitch"]=extraIMU_0_PC1s_pitch
    PC_dict["extraIMU_0"]["PC1s"]["roll"]=extraIMU_0_PC1s_roll
    PC_dict["extraIMU_0"]["PC2s"]["pitch"]=extraIMU_0_PC2s_pitch
    PC_dict["extraIMU_0"]["PC2s"]["roll"]=extraIMU_0_PC2s_roll

    PC_dict["extraIMU_1"]["PC1s"]["pitch"]=extraIMU_1_PC1s_pitch
    PC_dict["extraIMU_1"]["PC1s"]["roll"]=extraIMU_1_PC1s_roll
    PC_dict["extraIMU_1"]["PC2s"]["pitch"]=extraIMU_1_PC2s_pitch
    PC_dict["extraIMU_1"]["PC2s"]["roll"]=extraIMU_1_PC2s_roll
    
    return PC_dict

def Empirical_Distribution_Function(PC_dict,binwidth=6): 
    """Takes in the PC_dict from two_PCs_daily_scatter() and returns the default yaw value for each IMU"""
    IMUs=["IMU","extraIMU_0","extraIMU_1"]
    x_axis_ship=np.array((0,1))
    y_axis_ship=np.array((1,0))
    dist_dict={}
    fig,ax=plt.subplots(2,3,figsize=(35,15))
    plt_cnt=0
    yaw_offset = []
    peak_val = []
    for imu in IMUs:
        dist_dict[imu]={}
        dist_dict[imu]["PC1_angles"]=[]
        dist_dict[imu]["PC2_angles"]=[]
        for cnt in range(len(PC_dict[imu]["PC1s"]["roll"])):
            PC1v=np.array((PC_dict[imu]["PC1s"]["roll"][cnt],PC_dict[imu]["PC1s"]["pitch"][cnt]))
            PC2v=np.array((PC_dict[imu]["PC2s"]["roll"][cnt],PC_dict[imu]["PC2s"]["pitch"][cnt]))
            PC1_dot=np.dot(x_axis_ship,PC1v)
            PC2_dot=np.dot(-y_axis_ship,PC2v)
            PC1_angle=np.arccos(PC1_dot)/np.pi*180
            PC2_angle=np.arccos(PC2_dot)/np.pi*180
            dist_dict[imu]["PC1_angles"].append(PC1_angle-90)
            dist_dict[imu]["PC2_angles"].append(PC2_angle-90)
            
        bin_vals=np.arange(-90,90,binwidth)
        bin_mids=np.arange(-90+0.5*binwidth,90-0.5*binwidth,binwidth)
        PC2_PDF,whatever=np.histogram(dist_dict[imu]["PC2_angles"],bins=bin_vals)
        PC1_PDF,whatever=np.histogram(dist_dict[imu]["PC1_angles"],bins=bin_vals)
        PC2_interpol=interpolate.interp1d(bin_mids,PC2_PDF/len(dist_dict[imu]["PC2_angles"]),kind="quadratic")
        PC1_interpol=interpolate.interp1d(bin_mids,PC1_PDF/len(dist_dict[imu]["PC1_angles"]),kind="quadratic")
        plt_x=np.arange(-90+2*binwidth,90-2*binwidth,0.1)
        PC1_plt_y=PC1_interpol(plt_x)
        PC2_plt_y=PC2_interpol(plt_x)
        
        ax[0,plt_cnt].plot(plt_x,PC1_plt_y,"b")
        ax[0,plt_cnt].set_xlabel("angle in degrees perpendicular midships")
        ax[0,plt_cnt].set_ylabel("probabiliy")
        title_str="Empirical Distribution Fuction of PC1 from "+imu
        ax[0,plt_cnt].set_title(title_str)
        
        ax[1,plt_cnt].plot(plt_x,PC2_plt_y,"b")
        ax[1,plt_cnt].set_xlabel("angle in degrees from midships")
        ax[1,plt_cnt].set_ylabel("probabiliy")
        title_str="Empirical Distribution Fuction of PC2 from "+imu
        ax[1,plt_cnt].set_title(title_str)
       
        peak_val=np.argmax(PC2_plt_y)
        
        if peak_val>0:
            peak_angle=plt_x[peak_val]
            yaw_offset.append(peak_angle)
        else:
            peak_angle=np.nan
            yaw_offset.append(peak_angle)

        plt_cnt+=1
        
    yaw_dict={}
    yaw_dict["IMU"]={}
    yaw_dict["extraIMU_0"]={}
    yaw_dict["extraIMU_1"]={}
    yaw_dict["IMU"]["yaw_offset"]=yaw_offset[0]
    yaw_dict["extraIMU_0"]["yaw_offset"]=yaw_offset[1]
    yaw_dict["extraIMU_1"]["yaw_offset"]=yaw_offset[2]
        
    return yaw_dict

def get_fit_coefficients():
    """Goes trough all files in the folder base_path and returns the coefficients to fit extra IMU data to main IMU data"""
    base_path = "/Users/judiththuolberg/master/data"
    c0 = []
    c1 = []
    c2 = []
    c3 = []
    c4 = []
    
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
            extraIMU_coeffs = simple_trig_curve_fit(dict_data,'extraIMU_0',False)
            if len(extraIMU_coeffs)>0:
                c0.append(extraIMU_coeffs[0])
                c1.append(extraIMU_coeffs[1])
                c2.append(extraIMU_coeffs[2])
                c3.append(extraIMU_coeffs[3])
                c4.append(extraIMU_coeffs[4])   
    
    c0_val=np.nanmean(c0)
    c1_val=np.nanmean(c1)
    c2_val=np.nanmean(c2)
    c3_val=np.nanmean(c3)
    c4_val=np.nanmean(c4)
    
    return (c0_val,c1_val,c2_val,c3_val,c4_val)

def simple_trig_curve_fit(dict_data,IMU_name,plot):
    """Simple trigonometric function used to fit extra IMU data to main IMU data""" 
    IMU_accel_D=np.array((dict_data["list_time_series_interpolated"]["IMU"]["accel_D"]))
    extraIMU_0_accel_D=np.array((dict_data["list_time_series_interpolated"]["extraIMU_0"]["accel_D"]))
    extraIMU_1_accel_D=np.array((dict_data["list_time_series_interpolated"]["extraIMU_1"]["accel_D"]))
    pitch=np.array((dict_data["list_time_series_interpolated"]["IMU"]["pitch"]))
    roll=np.array((dict_data["list_time_series_interpolated"]["IMU"]["roll"]))
    date_time=np.array((dict_data["list_time_series_interpolated"]["common_datetime"]))
    pitch, roll = corr_pitch_n_roll_array("IMU",pitch,roll,date_time)
    roll=roll/360*(2*np.pi)
    pitch=pitch/360*(2*np.pi)
    
    if IMU_name == "extraIMU_0":
        extraIMU_accel_D = extraIMU_0_accel_D
    elif IMU_name == "extraIMU_1":   
        extraIMU_accel_D = extraIMU_1_accel_D
        
    IMU_accel_D_mean=np.mean(IMU_accel_D)
    extraIMU_accel_D_mean=np.mean(extraIMU_accel_D)
    
    IMU_z_tt=-(IMU_accel_D-IMU_accel_D_mean)
    extraIMU_z_tt=(extraIMU_accel_D-extraIMU_accel_D_mean)
    sample_freq=10
                
    IMU_z=fft_int_2(IMU_z_tt,sample_freq)
    extraIMU_z=fft_int_2(extraIMU_z_tt,sample_freq)
    
    data_len=pitch.shape[0]
    cutoff=data_len//10
    IMU_z=IMU_z[cutoff:-cutoff]
    extraIMU_z=extraIMU_z[cutoff:-cutoff]
    pitch=pitch[cutoff:-cutoff]
    roll=roll[cutoff:-cutoff]
    
    extraIMU_z_corr=extraIMU_z-IMU_z
    extraIMU_z_corr=extraIMU_z_corr[np.isfinite(extraIMU_z_corr)]

    def simple_trig_func(prz,c0,c1,c2,c3,c4):
        p,r,z=prz
        ret=c0*np.sin(-p)+c1*np.cos(r)+c2*np.cos(p)+c3*np.sin(r)+c4*z
        return ret

    if extraIMU_z_corr.shape[0]>0:
        extraIMU_coeffs,waste=curve_fit(simple_trig_func,(pitch,roll,IMU_z),extraIMU_z_corr)
        extraIMU_z_comp=simple_trig_func((pitch,roll,IMU_z),extraIMU_coeffs[0],extraIMU_coeffs[1],extraIMU_coeffs[2],extraIMU_coeffs[3],extraIMU_coeffs[4])
    
        extraIMU_sigma=np.sqrt(np.sum((extraIMU_z_corr-extraIMU_z_comp)**2)/extraIMU_z_corr.shape[0])
        extraIMU_z_comp_save=extraIMU_z_comp[np.where(np.abs(extraIMU_z_corr)>0.1)]
        extraIMU_z_corr_save=extraIMU_z_corr[np.where(np.abs(extraIMU_z_corr)>0.1)]
        extraIMU_z_save2=extraIMU_z[np.where(np.abs(extraIMU_z>0.1))]
        extraIMU_z_corr_save2=extraIMU_z_corr[np.where(np.abs(extraIMU_z>0.1))]
        extraIMU_z_comp_save2=extraIMU_z_comp[np.where(np.abs(extraIMU_z>0.1))]
        extraIMU_relative_mean_error=np.sum(np.abs(extraIMU_z_corr_save-extraIMU_z_comp_save)/np.abs(extraIMU_z_corr_save))/extraIMU_z_corr_save.shape[0]
        extraIMU_sigma_rel=extraIMU_sigma/np.mean(np.abs(extraIMU_z_corr))
        extraIMU_sigma_rel_ztot=extraIMU_sigma/np.mean(np.abs(extraIMU_z))
        extraIMU_sigma_err_rel=np.sqrt(np.sum(((extraIMU_z_corr_save-extraIMU_z_comp_save)/extraIMU_z_corr_save)**2)/extraIMU_z_corr_save.shape[0])
        extraIMU_sigma_err_rel_ztot=np.sqrt(np.sum(((extraIMU_z_corr_save2-extraIMU_z_comp_save2)/extraIMU_z_save2)**2)/extraIMU_z_corr_save2.shape[0])
        #print("extraIMU_1_sigma:",extraIMU_sigma, "relative mean error:",extraIMU_relative_mean_error,"sigma rel:",extraIMU_sigma_rel,"sigma rel to total z:",extraIMU_sigma_rel_ztot,"sigma of rel error:",extraIMU_sigma_err_rel,"sigma of err rel to total z:",extraIMU_sigma_err_rel_ztot)

        if plot==True:
            plt_len=3000
            plt.figure(0)
            plt.plot(np.arange(0,plt_len,1),extraIMU_z_corr[0:plt_len],"b",label="z measured")
            plt.plot(np.arange(0,plt_len,1),extraIMU_z_comp[0:plt_len],"r",label="z approximated")
            plt.title("z at beginning of data")
            plt.legend()
    
            plt_len=1000
            plt.figure(1)
            plt.plot(np.arange(0,plt_len,1),extraIMU_z_corr[0:plt_len],"b",label="z measured")
            plt.plot(np.arange(0,plt_len,1),extraIMU_z_comp[0:plt_len],"r",label="z approximated")
            plt.title("z at beginning of data")
            plt.legend()
    
            diff_err=np.abs(extraIMU_z_corr-extraIMU_z_comp)/(np.abs(extraIMU_z_corr)+1)
            min_ind=np.argmax(diff_err)
            plt_len=1000
            plt.figure(2)
            plt.plot(np.arange(0,plt_len,1),extraIMU_z_corr[min_ind-plt_len//2:min_ind+plt_len//2],"b",label="z measured")
            plt.plot(np.arange(0,plt_len,1),extraIMU_z_comp[min_ind-plt_len//2:min_ind+plt_len//2],"r",label="z approximated")
            plt.title("z at highest error area")
            plt.legend()
    
            plt_len=1000
            plt.figure(3)
            plt.plot(np.arange(0,plt_len,1),extraIMU_z[0:plt_len],"k",label="total z at extraIMU_1")
            plt.plot(np.arange(0,plt_len,1),extraIMU_z_comp[0:plt_len]+IMU_z[0:plt_len],"g",label="approximated z + z at VN100")
            plt.title("Comparison of total zs at beginning of data")
            plt.legend()
        
            plt_len=1000
            plt.figure(4)
            plt.plot(np.arange(0,plt_len,1),extraIMU_z[min_ind-plt_len//2:min_ind+plt_len//2],"k",label="z measured")
            plt.plot(np.arange(0,plt_len,1),extraIMU_z_comp[min_ind-plt_len//2:min_ind+plt_len//2]+IMU_z[min_ind-plt_len//2:min_ind+plt_len//2],"g",label="z approximated")
            plt.title("total zs at highest error area")
            plt.legend()
    else:
        print("extra IMU value does not exist")
        extraIMU_coeffs = []
    
    return extraIMU_coeffs

def read_raw_data():
    """Goes trough all files in folder base_path and returns the raw values in a dictionary"""
    pitch_IMU=[]
    roll_IMU=[]
    yaw_IMU=[]
    accel_D_IMU=[]
    pitch_extraIMU_0=[]
    roll_extraIMU_0=[]
    yaw_extraIMU_0=[]
    accel_D_extraIMU_0=[]
    pitch_extraIMU_1=[]
    roll_extraIMU_1=[]
    yaw_extraIMU_1=[]
    accel_D_extraIMU_1=[]
    ug_1=[]
    ug_2=[]
    radar=[]
    date_time=[]

    dir_cnt=0
    base_path = "/Users/judiththuolberg/master/data"
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
        
            pitch_IMU.extend(dict_data["list_time_series_interpolated"]["IMU"]["pitch"])
            roll_IMU.extend(dict_data["list_time_series_interpolated"]["IMU"]["roll"])
            yaw_IMU.extend(dict_data["list_time_series_interpolated"]["IMU"]["yaw"])
            accel_D_IMU.extend(dict_data["list_time_series_interpolated"]["IMU"]["accel_D"])
            pitch_extraIMU_0.extend(dict_data["list_time_series_interpolated"]["extraIMU_0"]["pitch"])
            roll_extraIMU_0.extend(dict_data["list_time_series_interpolated"]["extraIMU_0"]["roll"])
            yaw_extraIMU_0.extend(dict_data["list_time_series_interpolated"]["extraIMU_0"]["yaw"])
            accel_D_extraIMU_0.extend(dict_data["list_time_series_interpolated"]["extraIMU_0"]["accel_D"])
            pitch_extraIMU_1.extend(dict_data["list_time_series_interpolated"]["extraIMU_1"]["pitch"])
            roll_extraIMU_1.extend(dict_data["list_time_series_interpolated"]["extraIMU_1"]["roll"])
            yaw_extraIMU_1.extend(dict_data["list_time_series_interpolated"]["extraIMU_1"]["yaw"])
            accel_D_extraIMU_1.extend(dict_data["list_time_series_interpolated"]["extraIMU_1"]["accel_D"])
            ug_1.extend(dict_data["list_time_series_interpolated"]["gauges"]["gauge_ug1_meters"])
            ug_2.extend(dict_data["list_time_series_interpolated"]["gauges"]["gauge_ug2_meters"])
            radar.extend(dict_data["list_time_series_interpolated"]["gauges"]["gauge_radar1_meters"])
            date_time.extend(dict_data["list_time_series_interpolated"]["common_datetime"])

    ret_dict={}
    ret_dict["IMU"]={}
    ret_dict["IMU"]["pitch"]=pitch_IMU
    ret_dict["IMU"]["roll"]=roll_IMU
    ret_dict["IMU"]["yaw"]=yaw_IMU
    ret_dict["IMU"]["accel_D"]=accel_D_IMU
    ret_dict["extraIMU_0"]={}
    ret_dict["extraIMU_0"]["pitch"]=pitch_extraIMU_0
    ret_dict["extraIMU_0"]["roll"]=roll_extraIMU_0
    ret_dict["extraIMU_0"]["yaw"]=yaw_extraIMU_0
    ret_dict["extraIMU_0"]["accel_D"]=accel_D_extraIMU_0
    ret_dict["extraIMU_1"]={}
    ret_dict["extraIMU_1"]["pitch"]=pitch_extraIMU_1
    ret_dict["extraIMU_1"]["roll"]=roll_extraIMU_1
    ret_dict["extraIMU_1"]["yaw"]=yaw_extraIMU_1
    ret_dict["extraIMU_1"]["accel_D"]=accel_D_extraIMU_1
    ret_dict["gauges"]={}
    ret_dict["gauges"]["ug1"]=ug_1
    ret_dict["gauges"]["ug2"]=ug_2
    ret_dict["gauges"]["radar"]=radar
    ret_dict["datetime"]=date_time
            
    return ret_dict

def read_raw_data_means():
    """Goes trough all files in folder base_path and returns the mean value of each file in a dictionary"""
    pitch_IMU=[]
    roll_IMU=[]
    yaw_IMU=[]
    accel_D_IMU=[]
    pitch_extraIMU_0=[]
    roll_extraIMU_0=[]
    yaw_extraIMU_0=[]
    accel_D_extraIMU_0=[]
    pitch_extraIMU_1=[]
    roll_extraIMU_1=[]
    yaw_extraIMU_1=[]
    accel_D_extraIMU_1=[]
    ug_1=[]
    ug_2=[]
    radar=[]
    date_time=[]

    dir_cnt=0
    base_path = "/Users/judiththuolberg/master/data"
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
        
            pitch_IMU.append(np.nanmean(dict_data["list_time_series_interpolated"]["IMU"]["pitch"]))
            roll_IMU.append(np.nanmean(dict_data["list_time_series_interpolated"]["IMU"]["roll"]))
            yaw_IMU.append(np.nanmean(dict_data["list_time_series_interpolated"]["IMU"]["yaw"]))
            accel_D_IMU.append(np.nanmean(dict_data["list_time_series_interpolated"]["IMU"]["accel_D"]))
            pitch_extraIMU_0.append(np.nanmean(dict_data["list_time_series_interpolated"]["extraIMU_0"]["pitch"]))
            roll_extraIMU_0.append(np.nanmean(dict_data["list_time_series_interpolated"]["extraIMU_0"]["roll"]))
            yaw_extraIMU_0.append(np.nanmean(dict_data["list_time_series_interpolated"]["extraIMU_0"]["yaw"]))
            accel_D_extraIMU_0.append(np.nanmean(dict_data["list_time_series_interpolated"]["extraIMU_0"]["accel_D"]))
            pitch_extraIMU_1.append(np.nanmean(dict_data["list_time_series_interpolated"]["extraIMU_1"]["pitch"]))
            roll_extraIMU_1.append(np.nanmean(dict_data["list_time_series_interpolated"]["extraIMU_1"]["roll"]))
            yaw_extraIMU_1.append(np.nanmean(dict_data["list_time_series_interpolated"]["extraIMU_1"]["yaw"]))
            accel_D_extraIMU_1.append(np.nanmean(dict_data["list_time_series_interpolated"]["extraIMU_1"]["accel_D"]))
            ug_1.append(np.nanmean(dict_data["list_time_series_interpolated"]["gauges"]["gauge_ug1_meters"]))
            ug_2.append(np.nanmean(dict_data["list_time_series_interpolated"]["gauges"]["gauge_ug2_meters"]))
            radar.append(np.nanmean(dict_data["list_time_series_interpolated"]["gauges"]["gauge_radar1_meters"]))
            date_time.append(dict_data["list_time_series_interpolated"]["common_datetime"][0])

    ret_dict={}
    ret_dict["IMU"]={}
    ret_dict["IMU"]["pitch"]=pitch_IMU
    ret_dict["IMU"]["roll"]=roll_IMU
    ret_dict["IMU"]["yaw"]=yaw_IMU
    ret_dict["IMU"]["accel_D"]=accel_D_IMU
    ret_dict["extraIMU_0"]={}
    ret_dict["extraIMU_0"]["pitch"]=pitch_extraIMU_0
    ret_dict["extraIMU_0"]["roll"]=roll_extraIMU_0
    ret_dict["extraIMU_0"]["yaw"]=yaw_extraIMU_0
    ret_dict["extraIMU_0"]["accel_D"]=accel_D_extraIMU_0
    ret_dict["extraIMU_1"]={}
    ret_dict["extraIMU_1"]["pitch"]=pitch_extraIMU_1
    ret_dict["extraIMU_1"]["roll"]=roll_extraIMU_1
    ret_dict["extraIMU_1"]["yaw"]=yaw_extraIMU_1
    ret_dict["extraIMU_1"]["accel_D"]=accel_D_extraIMU_1
    ret_dict["gauges"]={}
    ret_dict["gauges"]["ug1"]=ug_1
    ret_dict["gauges"]["ug2"]=ug_2
    ret_dict["gauges"]["radar"]=radar
    ret_dict["datetime"]=date_time
            
    return ret_dict

