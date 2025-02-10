function plotIntermediateConfigurations(stim_params,n, ic_bs_ex,ic_bs_inh,ic_inh_lo,ic_inh_hi,ic, Fs, ic_bs_ex_hi, ic_bs_inh_hi)
% 
% figure('Position',[52,164,560,1036])
% tiledlayout(7, 1)
% 
% nexttile
% t = linspace(0, length(stim_params.stim(n,:))./Fs, length(stim_params.stim(n,:)));
% plot(t, stim_params.stim(n,:))
% title('Stimulus')
% xlim([0 1])
% 
% nexttile
% t = linspace(0, length(ic_bs_ex)./Fs, length(ic_bs_ex));
% plot(t, ic_bs_ex, "LineWidth",1.5)
% title('Exc On')
% xlim([0 1])
% 
% nexttile
% t = linspace(0, length(ic_bs_inh)./Fs, length(ic_bs_inh));
% plot(t, ic_bs_inh, "LineWidth",1.5)
% title('Inh On')
% xlim([0 1])
% 
% nexttile
% t = linspace(0, length(ic_inh_lo)./Fs, length(ic_inh_lo));
% plot(t, ic_inh_lo, "LineWidth",1.5)
% title('Inh Low')
% xlim([0 1])
% 
% nexttile
% t = linspace(0, length(ic_inh_hi)./Fs, length(ic_inh_hi));
% plot(t, ic_inh_hi, "LineWidth",1.5)
% title('Inh High')
% xlim([0 1])
% 
% nexttile
% t = linspace(0, length(ic)./Fs, length(ic));
% plot(t, ic, "LineWidth",1.5)
% title('Final')
% xlim([0 1])
% 
% nexttile
% plot(t, (ic + abs(ic))/2, "LineWidth",1.5)
% title('Half-Wave Rectified')
% xlim([0 1])


%% 

figure('Position',[10,179,521,702])
tiledlayout(4, 1)

a(1) = nexttile;
t = linspace(0, length(stim_params.stim(n,:))./Fs, length(stim_params.stim(n,:)));
plot(t, stim_params.stim(n,:))
%title(['Stimulus: ' num2str(stim_params.all_fms(n)) ' Hz'])
xlim([0 1])

t_on = min([length(ic_bs_inh) length(ic_inh_hi) length(ic_inh_lo)]);
t = linspace(0, length(ic_bs_ex)./Fs, length(ic_bs_ex));
t_inh = linspace(0, t_on./Fs, t_on);

% a(2) = nexttile;
% hold on
% plot(t, ic_bs_ex, "LineWidth",1.5)
% plot(t_inh, ic_bs_inh(1:t_on), "LineWidth",1.5)
% xlim([0 1])
% legend('On-CF Excitation', 'On-CF Inhibition')
% ylim([0 250])

a(2) = nexttile(2, [2 1]);
hold on
%plot(t, ic_bs_ex-ic_bs_inh, "LineWidth",1.5)
yline(600, 'k')
plot(t, ic_bs_ex+600, "LineWidth",1.5)
dc_component(1) = getACDCratio(ic_bs_ex);

yline(200, 'k')
plot(t_inh, ic_bs_inh(1:t_on)+200, "LineWidth",1.5)
dc_component(2) = getACDCratio(ic_bs_inh(1:t_on));

yline(100, 'k')
plot(t_inh, ic_inh_lo(1:t_on)+100, "LineWidth",1.5)
dc_component(3) = getACDCratio(ic_inh_lo(1:t_on));

yline(0, 'k')
plot(t_inh, ic_inh_hi(1:t_on), "LineWidth",1.5)
dc_component(4) = getACDCratio(ic_inh_hi(1:t_on));
xlim([0 1])


legend('', ['Exc On, ' num2str(dc_component(1))],...
	'', ['On Inh, ' num2str(dc_component(2))],...
	'', ['Low Inh, ' num2str(dc_component(3))],...
	'', ['High Inh, ' num2str(dc_component(4))])
%ylim([0 100])


a(3) = nexttile;
t = linspace(0, length(ic)./Fs, length(ic));
plot(t, (ic + abs(ic))/2, "LineWidth",1.5)
title('Half-Wave Rectified, Broad')
xlim([0 1])
avg_rate = mean((ic + abs(ic))/2);
legend(['Mean = ' num2str(avg_rate)])
ylim([0 110])

linkaxes(a, 'x')


	function dc_component = getACDCratio(signal)
		dc_component = mean(signal);
		ac_component = signal - dc_component;
		dc_power = dc_component^2;
		ac_power = sum(ac_component.^2) / length(ac_component);
		ac_dc_ratio = ac_power / dc_power;
	end




%% Plot intermediate another way 
% 
% figure
% tiledlayout(5, 1)
% 
% a(1) = nexttile;
% t = linspace(0, length(stim_params.stim(n,:))./Fs, length(stim_params.stim(n,:)));
% plot(t, stim_params.stim(n,:))
% %title(['Stimulus: ' num2str(stim_params.all_fms(n)) ' Hz'])
% xlim([0 1])
% 
% t_on = min([length(ic_bs_inh) length(ic_inh_hi) length(ic_inh_lo)]);
% t = linspace(0, length(ic_bs_ex)./Fs, length(ic_bs_ex));
% t_inh = linspace(0, t_on./Fs, t_on);
% 
% a(2) = nexttile;
% hold on
% plot(t, ic_bs_ex, "LineWidth",1.5)
% plot(t_inh, ic_bs_inh(1:t_on), "LineWidth",1.5)
% xlim([0 1])
% legend('On-CF Excitation', 'On-CF Inhibition')
% ylim([0 250])
% 
% a(3) = nexttile;
% hold on
% plot(t, ic_bs_ex_hi(1:t_on), "LineWidth",1.5)
% plot(t_inh, ic_bs_inh_hi(1:t_on), "LineWidth",1.5)
% plot(t_inh,  ic_bs_ex_hi(1:t_on)-ic_bs_inh_hi(1:t_on), "LineWidth",1.5)
% xlim([0 1])
% legend('High-CF Excitation', 'High-CF Inhibition', 'Total')
% ylim([0 250])
% 
% a(4) = nexttile;
% hold on
% %plot(t, ic_bs_ex-ic_bs_inh, "LineWidth",1.5)
% %plot(t_inh, ic_inh_lo(1:t_on), "LineWidth",1.5)
% plot(t_inh, ic_inh_hi(1:t_on), "LineWidth",1.5)
% plot(t_inh, ic_bs_ex_hi(1:t_on)-ic_bs_inh_hi(1:t_on), "LineWidth",1.5)
% 
% xlim([0 1])
% legend('BS On', 'Low Inh', 'High Inh')
% %ylim([0 100])
% 
% 
% a(5) = nexttile;
% t = linspace(0, length(ic)./Fs, length(ic));
% plot(t, (ic + abs(ic))/2, "LineWidth",1.5)
% title('Half-Wave Rectified, Broad')
% xlim([0 1])
% avg_rate = mean((ic + abs(ic))/2);
% legend(['Mean = ' num2str(avg_rate)])
% ylim([0 110])
% 
% linkaxes(a, 'x')

end