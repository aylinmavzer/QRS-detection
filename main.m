clear all
close all
clc 
rng(0, 'twister') %seed
data1= load('data.mat');
sig1 = data1.ecg_data;
fs1  = data1.fs;

pan_tompkins(sig1, fs1, 0.25, 0)

sig2 = awgn(sig1 , 18); %Add white Gaussian Noise. I added white gaussian noise because of the efficiency. It is between -0.5,0.5 in the time domain as you can see in the plot.

pan_tompkins(sig2, fs1, 0.25, 1)

function []=pan_tompkins(sig, fs, th, n)
    tensecofdata=sig(1:1000);%fs =100, 100x10sec=1000 sample
    subplot(2,2,1+n);
    plot((1:length(tensecofdata))./fs,tensecofdata);
    if n==0
    title('First 10s of raw ECG data')
    xlabel('Time (second)')
    ylabel('Magnitude')
    else
        title('First 10s of ECG data with noise')
        xlabel('Time (second)')
        ylabel('Magnitude')
    end
    
    %FILTERS
    %Low Pass Filter 
    B=[1 0 0 0 0 0 -2 0 0 0 0 0 1]; %Transfer Function's numerator array according to article
    A=[1 -2 1];
    sig_L=filter(B,A,sig(1:540000));
    sig_L = sig_L./max(sig_L);
    
    % High Pass Filter
    B=[-1/32 zeros(1,15) 1 -1 zeros(1,14) 1/32];
    A=[1 -1];
    sig_h=filter(B,A,sig_L);
    
    
    %Derivative
    B=1/8*[2 1 0 -1 -2];
    A=[1];
    sig_d=filter(B,A,sig_h);
    
    
    % Squaring Function
    sig_s=sig_d.^2;
    
   
    % Moving-Window Integration
    N=round(0.15*fs);
    B=1/N*ones(1,N);
    A=[1];
    sig_w=filter(B,A,sig_s);
    sig_w = sig_w./max(sig_w);
    
    %Threshold

    sig_th = arrayfun(@(x)(filter_th(x,th)), sig_w);
   
    
    %Determining MAX value of the each signal and creating diracs 
    last_n = 4;
    sig_dirac = zeros(length(sig_th),1);
    in_block = false;
    max_val = 0;
    max_pos = 0;
    for i = (last_n+1):length(sig_th)
        if(all(sig_th(i-last_n:i-1) == 0) && sig_th(i) ~= 0) 
            in_block = true; 
            max_val = sig_th(i);
            max_pos = i;
        elseif(all(sig_th(i-last_n:i-1) ~= 0) && sig_th(i) == 0)
            in_block = false;
            sig_dirac(max_pos) = 1;
        else
           if(sig_th(i) > max_val)
              max_val = sig_th(i);
              max_pos = i;
           end
        end
    end

    
   %Number of beat per minute
    min =0;
    for t=1:1:540000
        if(mod(t,6000) == 0)
            min = min + 1;
            hr(min) = sum(sig_dirac((min-1) * 6000 + 1: min * 6000));
            
        end       
    end
    subplot(2,2,3+n)
    plot(hr)
     if n==0
    title('Heart Rate')
    xlabel('Time (minutes)')
    ylabel('Number of Heart beats')
     else
        title('Heart Rate')
        xlabel('Time (minutes)')
        ylabel('Number of Heart beats')
     end
    %number of QRS in interval of 45 to 50 mins
    if n == 0
        disp(['number of QRS in interval of 45 to 50 mins on original ECG signal= ',num2str(sum(sig_dirac(270000:300000)))])
    else
        disp(['number of QRS in interval of 45 to 50 mins on ECG signal with noise= ',num2str(sum(sig_dirac(270000:300000)))])
    end
    
   
end

    
function [y]=filter_th(x,th)
    if(x > th)
        y = x;
    else
        y = 0;
    end
end



    