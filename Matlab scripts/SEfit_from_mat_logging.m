clf;
addpath('mat\')
% N3, NN7, NN4, 16/64, four dates each
n=0;
n=n+1;file_in{n}='QL40_denoised_16_fev22_NN4.mat';
n=n+1;file_in{n}='QL40_denoised_16_fev22_NN7.mat';
n=n+1;file_in{n}='QL40_denoised_16_july22_NN3.mat';
n=n+1;file_in{n}='QL40_denoised_16_july22_NN4.mat';
n=n+1;file_in{n}='QL40_denoised_16_july22_NN7.mat';
n=n+1;file_in{n}='QL40_denoised_16_mars21_NN3.mat';
n=n+1;file_in{n}='QL40_denoised_16_mars21_NN4.mat';
n=n+1;file_in{n}='QL40_denoised_16_mars21_NN7.mat';
n=n+1;file_in{n}='QL40_denoised_16_NN7_28Sep20_FW_complete.mat';
n=n+1;file_in{n}='QL40_denoised_16_nov21_NN3.mat';
n=n+1;file_in{n}='QL40_denoised_16_sep20_NN3.mat';
n=n+1;file_in{n}='QL40_denoised_16_sep20_NN4.mat';
% n=n+1;file_in{n}='QL40_denoised_64_fev22_NN4.mat';
n=n+1;file_in{n}='QL40_denoised_64_fev22_NN7.mat';
n=n+1;file_in{n}='QL40_denoised_64_july22_NN3.mat';
n=n+1;file_in{n}='QL40_denoised_64_july22_NN4.mat';
n=n+1;file_in{n}='QL40_denoised_64_july22_NN7.mat';
n=n+1;file_in{n}='QL40_denoised_64_mars21_NN3.mat';
n=n+1;file_in{n}='QL40_denoised_64_mars21_NN4.mat';
n=n+1;file_in{n}='QL40_denoised_64_mars21_NN7.mat';
n=n+1;file_in{n}='QL40_denoised_64_NN7_28Sep20_FW_complete.mat';
n=n+1;file_in{n}='QL40_denoised_64_nov21_NN3.mat';
n=n+1;file_in{n}='QL40_denoised_64_sep20_NN3.mat';
n=n+1;file_in{n}='QL40_denoised_64_sep20_NN4.mat';

depths_cell=cell(n,1);
Mapp_cell=cell(n,1);
IP_flag_cell=cell(n,1);
t_center=cell(n,1);
ydata=cell(n,1);
xdata=cell(n,1);
depth_data=cell(n,1);
current_cell=cell(n,1);
num_dp=NaN(n,1);
num_gate=NaN(n,1);
% col=turbo(n);

num_depth_ref =1643; %max num_depths sur les 24 datasets importés
fitresult_se=cell(n,num_depth_ref);
gof_se=cell(n,num_depth_ref);
output_se=cell(n,num_depth_ref);
se_fit_mapp=cell(n,num_depth_ref);

fitresult_se_2=cell(n,num_depth_ref);
gof_se_2=cell(n,num_depth_ref);
output_se_2=cell(n,num_depth_ref);
se_fit_mapp_2=cell(n,num_depth_ref);

save_ydata_rloess=cell(n,num_depth_ref);
save_depthdata=NaN(n,num_depth_ref);
weights_SEfit=cell(n,num_depth_ref);
save_xdata=cell(n,num_depth_ref);
save_ydata=cell(n,num_depth_ref);

rmse_data=NaN(n,num_depth_ref);

Rho_save=cell(n,1);
Res_save=cell(n,1);
save_current=cell(n,1);

coeff_alpha_GEMGAS=NaN(n,num_depth_ref);
coeff_tauKWW_GEMGAS=NaN(n,num_depth_ref);
coeff_beta_GEMGAS=NaN(n,num_depth_ref);
rsquare_GEMGAS=NaN(n,num_depth_ref);

coeff_alpha=NaN(n,num_depth_ref);
coeff_tauKWW=NaN(n,num_depth_ref);
coeff_beta=NaN(n,num_depth_ref);

coeff_alpha_2=NaN(n,num_depth_ref);
coeff_tauKWW_2=NaN(n,num_depth_ref);
coeff_beta_2=NaN(n,num_depth_ref);

NGates=36; % including the first one
GateWidth=NaN(n,NGates);

%% launch the loop, start at chosen index 
name_current=cell(n,1);
for v=1:n
    if contains(file_in{v},'NN3')
        if contains(file_in{v},'sep20')
            name_current{v}='NN3_21Sep20_FW';
        elseif contains(file_in{v},'mars21')
            name_current{v}='NN3_10mars21_FW';
        elseif contains(file_in{v},'nov21')
            name_current{v}='NN3_15nov21_FW';
        elseif contains(file_in{v},'july22')
            name_current{v}='NN3_5july22_FW';
        end
    elseif contains(file_in{v},'NN4')
        if contains(file_in{v},'sep20')
            name_current{v}='NN4_22Sep20_FW';
        elseif contains(file_in{v},'mars21')
            name_current{v}='NN4_9mars21_FW';
        elseif contains(file_in{v},'fev22')
            name_current{v}='NN4_18fev22_FW';
        elseif contains(file_in{v},'july22')
            name_current{v}='NN4_7july22_FW';
        end
    elseif contains(file_in{v},'NN7')
        if contains(file_in{v},'Sep20')
            name_current{v}='NN7_28Sep20_FW_complete';
        elseif contains(file_in{v},'mars21')
            name_current{v}='NN7_11mars21_FW';
        elseif contains(file_in{v},'fev22')
            name_current{v}='NN7_18fev22_FW';
        elseif contains(file_in{v},'july22')
            name_current{v}='NN7_6july22_FW';
        end
    end
end
for k=1:n
    inch=str2double(file_in{k}(15:16));
    load(file_in{k})
    depths_cell{k}=depths;
    Mapp_cell{k}=Mapp;
    t_center{k}=t_center_pos;
    num_dp(k)=numel(depths_cell{k});
    depth_data{k}=NaN(1,num_dp(k));
    xdata{k}=NaN(NGates,num_dp(k));
    ydata{k}=NaN(NGates,num_dp(k));

    % name_var=['QL40_raw_' name_current{k} '.mat'];
    % load(name_var,'QL40_voltage','QL40_time','QL40_current');
    name_var=['QL40_ready_' name_current{k} '.mat'];
    load(name_var,'QL40_voltage','QL40_time','QL40_current','QL40_val');
    % store some values from QL40_val // different if input file is tx2
    if inch == 16
        Rho_save{k}=QL40_val{1};
        Res_save{k}=QL40_val{13};
    elseif inch == 64
        Rho_save{k}=QL40_val{2};
        Res_save{k}=QL40_val{14};
    end

    % remove low current first once and for all?
    % But only important for NN4since we dont keep the february 2022 in NN3
    depth_current=QL40_current; %1 is depth, 2 is current
    save_current{k}=depth_current;
    depth_vector_current=depth_current(:,1); %should be the same as for IP but maybe there is a shift
    current_vector=depth_current(:,2);
    current_cell{k,1}=NaN(num_dp(k),1);
    for i=1:num_dp(k)
        test_dp=1;
        depth_test=depths_cell{k}(i);
        depth_id_current= ismembertol(depth_vector_current,depth_test,0.125,'DataScale', 1);
        if sum(depth_id_current)==0
            warning(['depth=' num2str(depth_test) ' not found in current for file' num2str(k)])
        else
            current_cell{k,1}(i)=current_vector(depth_id_current);
            if current_vector(depth_id_current)<200
                warning(['depth=' num2str(depth_test) ' current lower than 200 mA for file' num2str(k)])
                Mapp_cell{k}(i,:)=NaN;
                test_dp=0;
            end
        end
        if test_dp==1
            xdata0_tmp = t_center{k}(1:end); 
            ydata0_tmp = Mapp_cell{k}(i,1:end)';
    
            ind_gate_noneg=find(ydata0_tmp>0);
            xdata0_noneg = xdata0_tmp(ind_gate_noneg);
            ydata0_noneg = ydata0_tmp(ind_gate_noneg);
    
            n_gate_start = 1;
%             xdata0=xdata0_noneg(n_gate_start:end); %removing first gate only (0 gates in the other...)
%             ydata0=ydata0_noneg(n_gate_start:end);
            
            xdata0=log(xdata0_noneg(n_gate_start:end)); %removing first gate? not now
            ydata0=log(ydata0_noneg(n_gate_start:end));

            if numel(xdata0)>15
                [xData_exp, yData_exp] = prepareCurveData( xdata0, ydata0);
                ft = fittype( 'a1 - exp(c*x-b1)', 'independent', 'x', 'dependent', 'y' );
%                 ft = fittype( 'a*exp(-(x/b)^c)', 'independent', 'x', 'dependent', 'y' );
                opts = fitoptions( 'Method', 'NonlinearLeastSquares');
                opts.Display = 'Off';
                opts.Robust = 'LAR';
                opts.Lower = [0 -Inf 0.1];
                opts.StartPoint = [1 1 0.5];
                opts.Upper = [10 Inf 1];
                [fitresult_se{k,i}, gof_se{k,i}, output_se{k,i}] = fit( xData_exp, yData_exp, ft, opts );
    
                if gof_se{k,i}.rsquare<0.99
                    warning(['Bad SE fit at DP=' num2str(i) ' for file ' file_in{k}(1:end-4) '; try new starting point'])
                    opts.StartPoint = [5 5 0.5]
                    [fitresult_se{k,i}, gof_se{k,i}, output_se{k,i}] = fit( xData_exp, yData_exp, ft, opts );
                    if gof_se{k,i}.rsquare<0.99
                        warning(['(again) Bad SE fit at DP=' num2str(i) ' for file ' file_in{k}(1:end-4) '; try new starting point'])
                        opts.StartPoint = [0.6 0.2 0.7];
                        [fitresult_se{k,i}, gof_se{k,i}, output_se{k,i}] = fit( xData_exp, yData_exp, ft, opts );
%                         start_point=[1 0.4 0.1];
                    else
%                         start_point=[10 400 0.8];
                    end
                else
%                     start_point=[20 100 0.5];
                end
                start_point=opts.StartPoint;
                coeff_alpha(k,i)=exp(fitresult_se{k,i}.a1);
                coeff_tauKWW(k,i)=exp(fitresult_se{k,i}.b1/fitresult_se{k,i}.c);
                coeff_beta(k,i)=fitresult_se{k,i}.c;
            
                se_fit_mapp{k,i}=coeff_alpha(k,i)*exp(-(xdata0_noneg/coeff_tauKWW(k,i)).^coeff_beta(k,i));   
%                 se_fit_mapp{k,i}=fitresult_se{k,i}.a*exp(-(xdata0/fitresult_se{k,i}.b).^fitresult_se{k,i}.c); 

                plot(xdata0_noneg,ydata0_noneg,'o','Color',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',9);
                hold on
                plot (xdata0_noneg,se_fit_mapp{k,i},'-','Color','k','LineWidth',2)
%                 plot(xdata0,ydata0,'o','Color',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',9);
%                 hold on
%                 plot (xdata0,se_fit_mapp{k,i},'-','Color','k','LineWidth',2)
                
                set(gca, 'YScale', 'log')
                set(gca, 'XScale', 'log')
    
                hXLabel=xlabel('Time (msec)');
                hYLabel=ylabel('Polarizability \eta (mV/V)');
                xlim([1e0,1e4])
                ylim([1e-2,1e3])
                grid on;

                % try to adjust the fit based on residuals    
                res_thres=0.2;
                residuals_GEMGAS_tmp = output_se{k, i}.residuals;
    
%                 opts = fitoptions( 'Method', 'NonlinearLeastSquares');
%                 opts.Robust = 'LAR';
%                 opts.StartPoint = start_point;
%                 opts.Lower = [0 50 0];
%                 opts.Upper = [Inf 5000 1];
                points_to_keep=find(abs(residuals_GEMGAS_tmp)<res_thres);
                opts.Exclude=find(abs(residuals_GEMGAS_tmp)>res_thres);
                % opts.Weights=W;
%                 xdata0_excl=xdata0(points_to_keep);
                xdata0_excl=xdata0_noneg(points_to_keep);
                ydata0_excl=ydata0_noneg(points_to_keep);

                if numel(xdata0_excl)>5
                    [fitresult_se_2{k,i}, gof_se_2{k,i}, output_se_2{k,i}] = fit( xData_exp, yData_exp, ft, opts );
                    coeff_alpha_2(k,i)=exp(fitresult_se_2{k,i}.a1);
                    coeff_tauKWW_2(k,i)=exp(fitresult_se_2{k,i}.b1/fitresult_se_2{k,i}.c);
                    coeff_beta_2(k,i)=fitresult_se_2{k,i}.c;

                    se_fit_mapp_2{k,i}=coeff_alpha_2(k,i)*exp(-(xdata0_noneg/coeff_tauKWW_2(k,i)).^coeff_beta_2(k,i)); 
%                     se_fit_mapp_2{k,i}=fitresult_se_2{k,i}.a*exp(-(xdata0/fitresult_se_2{k,i}.b).^fitresult_se_2{k,i}.c);   
        
                    % plot
                    plot(xdata0_excl,ydata0_excl,'o','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'MarkerSize',5);
                    plot(xdata0_noneg,se_fit_mapp_2{k,i},'-','Color','b','LineWidth',1)
                    txt_title=['Depth:' num2str(depths_cell{k}(i)) 'm ; I =' num2str(round(current_cell{k,1}(i),0)) ' mA' '; R^2_1 =' num2str(gof_se{k,i}.rsquare,3) '; R^2_2 =' num2str(gof_se_2{k, i}.rsquare,3) ];
                    hLegend=legend('All data','Fit 1','Data 2','Fit 2','Location','Best');
                else
                    plot(xdata0_excl,ydata0_excl,'o','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'MarkerSize',5);
                    % plot (xdata0,se_fit_mapp_2{k,i},'-','Color','b','LineWidth',1)
                    % plot(xdata0_noflag,ydata0_noflag,'*','Color',[1 0 0],'MarkerSize',7);
                    txt_title=['Depth:' num2str(depths_cell{k}(i)) 'm ; I =' num2str(round(current_cell{k,1}(i),0)) ' mA' '; R^2_1 =' num2str(gof_se{k, i}.rsquare,3)];
                    hLegend=legend('All data','Fit 1','Data 2','Location','Best');
                end
                hold off
    
                % txt_title=['DP:' num2str(i) '; I =' num2str(round(1000*current_cell{k,1}(i),0)) ' mA' '; R^2_1 =' num2str(gof_se{k, i}.rsquare,3) '; R^2_2 =' num2str(gof_se_2{k, i}.rsquare,3) ];
                title(txt_title)
                set(gca,'FontName','Times New Roman','FontSize',14);
                set([hXLabel,hYLabel,hLegend],'FontName','Times New Roman','FontSize',14);
                % name_out=['png\SE_fit_GEMGAS\new_fit\' file_in{k}(1:end-4) '_DP' num2str(i) '_SEfit.png'];
                name_out=['png\SE_fit_GEMGAS\new_fit\' file_in{k}(1:end-4) '_Depth' num2str(depths_cell{k}(i)) '_SEfit.png'];
                exportgraphics(gcf,name_out,'Resolution',300)

            % store final fit results
                save_xdata{k,i}=xdata0_tmp;
                save_ydata{k,i}=ydata0_tmp;     
                if ~isempty(gof_se_2{k, i})
                    if gof_se{k, i}.rsquare>gof_se_2{k, i}.rsquare
                        coeff_alpha_GEMGAS(k,i) = coeff_alpha(k,i);
                        coeff_tauKWW_GEMGAS(k,i) = coeff_tauKWW(k,i);
                        coeff_beta_GEMGAS(k,i) = coeff_beta(k,i);
                        rsquare_GEMGAS(k,i) = gof_se{k, i}.rsquare;
                        first_or_second=1;
                    else
                        coeff_alpha_GEMGAS(k,i) = coeff_alpha_2(k,i);
                        coeff_tauKWW_GEMGAS(k,i) = coeff_tauKWW_2(k,i);
                        coeff_beta_GEMGAS(k,i) = coeff_beta_2(k,i);
                        rsquare_GEMGAS(k,i) = gof_se_2{k, i}.rsquare;
                        first_or_second=2;
                    end
                end
            else
                warning(['Not enough data points left after flags removed for dp=' num2str(i) 'in file' file_in{k}(1:end-4)])
            end
        end
    end
    % keyboard
end

% save mat file
name_final='mat\GEMGAS_SEfit_log.mat';
save(name_final,'first_or_second','fitresult_se','gof_se','output_se','se_fit_mapp','fitresult_se_2','gof_se_2','output_se_2','se_fit_mapp_2', 'save_ydata','save_xdata','Mapp_cell','IP_flag_cell','t_center','depths_cell','GateWidth','coeff_alpha_GEMGAS','coeff_tauKWW_GEMGAS','coeff_beta_GEMGAS','rsquare_GEMGAS','current_cell','Rho_save','Res_save','file_in')