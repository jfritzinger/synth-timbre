function calc_DV = generateDV(cond, DV_type, CFs, avBS, BS, Fs)
% Calculates the decision variable for several different cases
% J. Fritzinger, updated 1/6/2021
%
% Inputs: cond: pitch or timbre
%         DV_type: type of DV 
%         CFs: array of CFs
%         avBS: average BS response
%         BS: BS response
%         Fs: sampling rate, Hz

switch cond
    case 1 % timbre
        log_CFs = log10(CFs);
        if DV_type == 1 % DV 1
            % peaks
            
        elseif DV_type == 2  % DV 2, center of mass
            calc_DV = (sum(CFs.*avBS))/sum(avBS);
            
        elseif DV_type == 3 % DV 3, smoothed and subtracted center of mass
            smoothing = supsmu(log_CFs, avBS, 'Span', 0.5);
            y_d = avBS - smoothing;
            y_d(y_d<0)=0;
            calc_DV = 10^((sum(log_CFs.*y_d))/(sum(y_d)));
            
        elseif DV_type == 4 % DV 4, Braden's DV
            calc_DV = new_approach(avBS, CFs);
            
        elseif DV_type == 5 % DV 5, subtracting average center of mass
            mean_sub_BS = max(avBS - 28, 0);
            calc_DV = 10^((sum(log_CFs.*mean_sub_BS))/sum(mean_sub_BS));
        end
        
    case 2 % pitch
        %             new = -average_ic_sout_BS(inoise,:);
        %             [~, indices] = findpeaks(new, 'MinPeakProminence', 2);
        %             BE_peaks = CFs(indices);
        %             below_Fc = find(BE_peaks<1000);
        %             calc_DV(trial) = BE_peaks(below_Fc(end))-BE_peaks(below_Fc(end-1));
        mean_ic_sout_BS = mean(BS);
        Y = fft(mean_ic_sout_BS);
        P2 = abs(Y/length(mean_ic_sout_BS));
        P1 = P2(1:length(mean_ic_sout_BS)/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(length(mean_ic_sout_BS)/2))/length(mean_ic_sout_BS);
        
%         figure;
%         plot(f, P1)
%         title('Single-Sided Amplitude Spectrum of S(t)')
%         xlabel('f (Hz)')
%         ylabel('|P1(f)|')
        [~, I] = max(P1(5:end));
        calc_DV = f(I);
end
end