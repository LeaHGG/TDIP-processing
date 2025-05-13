addpath('DCIP_data\')
clf;
n=0;

%% fetch data
% Krafla 2017
% % Around KS1
n=n+1;file_in{n}='ISL1.tx2';
n=n+1;file_in{n}='ISL2.tx2';
n=n+1;file_in{n}='ISL3.tx2';
n=n+1;file_in{n}='ISL4.tx2';
% % % Around KH3
% n=n+1;file_in{n}='ISL5.tx2';
% n=n+1;file_in{n}='ISL6.tx2';
% n=n+1;file_in{n}='ISL7.tx2';
% n=n+1;file_in{n}='ISL8.tx2';
% n=n+1;file_in{n}='ISL9.tx2';
% % % Around KH1
n=n+1;file_in{n}='ISL10.tx2';
n=n+1;file_in{n}='ISL11.tx2';
n=n+1;file_in{n}='ISL12.tx2';

% Hvedemarken 2019
% n=n+1;file_in{n}='HXB_E6789_R3_M15_IPproc7g_AutoProc_noiseGF_reproc030720.tx2';
% n=n+1;file_in{n}='HXB_E6789_R4_IPproc7g_Res15_c_AutoProc.tx2';
% n=n+1;file_in{n}='HXB_E6789_R5_7g_fObStd_M_15_AutoProc_noiseGF.tx2';

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
col=turbo(n);

num_depth_ref =1678; %max num_dp sur les 12 profils Krafla
% num_depth_ref =2232; %max num_dp sur les 3 rounds E6789 Hvedemarken
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

coeff_alpha_Krafla=NaN(n,num_depth_ref);
coeff_tauKWW_Krafla=NaN(n,num_depth_ref);
coeff_beta_Krafla=NaN(n,num_depth_ref);
rsquare_Krafla=NaN(n,num_depth_ref);

coeff_alpha=NaN(n,num_depth_ref);
coeff_tauKWW=NaN(n,num_depth_ref);
coeff_beta=NaN(n,num_depth_ref);

coeff_alpha_2=NaN(n,num_depth_ref);
coeff_tauKWW_2=NaN(n,num_depth_ref);
coeff_beta_2=NaN(n,num_depth_ref);

%% store data
for k=1:n
    clf;
    data=importdata(file_in{k});
    data=data.data;
    Ndata=size(data,1);
    Rho=data(:,22);
    Res=data(:,21);
    flagdc=data(:,24);
    NGates=data(1,25);
    flagIP=data(:,end-8-NGates:end-9);
    delay=data(1,26+NGates);
    GateWidth=NaN(n,NGates);
    GateWidth(k,:) = data(1,26+NGates+1:26+2*NGates);
    GateStart=NaN(1,NGates);
    GateTime=NaN(1,NGates);
    GateStart(1)=delay;
    for i=2:NGates
        GateStart(i)=GateStart(i-1)+GateWidth(k,i-1);
    end
    for i=1:NGates
        % if GateWidth(i)>0 %some gate widths are at 0 or -1
            GateTime(i)=GateStart(i)+0.5*GateWidth(k,i);
        % else
        %     flagIP(i,:)<1
        % end
    end
    IP_Mapp = data(:,26:25+NGates);
    off=(data(end,end-3))/1000;
    on=(data(end,end-4)+data(end,end-5))/1000;
    Ax=data(:,1);
    Bx=data(:,2);
    Mx=data(:,3);
    Nx=data(:,4);
    Az=data(:,5);
    Bz=data(:,6);
    Mz=data(:,8);
    Nz=data(:,4);
    std_dc=data(:,23);
    std_ip=data(:,end-2*NGates-8:end-9-NGates);
    Current=data(:,end-7);
    
    current_cell{k}=Current;

    % test plotting
    % clf;
    % for i=1:20:Ndata
    %     ind_gate_non_flagged=find(flagIP(i,:)<1);
    %     plot(GateTime(ind_gate_non_flagged),IP_Mapp(i,ind_gate_non_flagged),'o-')
    %     hold on
    % end
    % hold off
    % xlim([1e0,1e4])
    % ylim([1e-1,1e3])
    % set(gca,'XScale','log');
    % set(gca,'YScale','log');
    % hXLabel=xlabel('Time (msec)');
    % hYLabel=ylabel('\eta (mV/V)');
    % grid on;
    % grid minor
    % axis square
    % set(gca,'FontName','Times New Roman','FontSize',14);
    % set([hXLabel,hYLabel],'FontName','Times New Roman','FontSize',14);

    % fileout=['png\Krafla_polarizability_processed_' file_in{k}(1:end-4) '.png'];
    % exportgraphics(gcf,fileout,'Resolution',300)

    % store values
    depths_cell{k}=[Ax Bx Mx Nx Az Bz Mz Nz];
    Mapp_cell{k}=IP_Mapp;
    IP_flag_cell{k}=flagIP;
    t_center{k}=GateTime;
    num_dp(k)=Ndata;
    num_gate(k)=NGates;
    xdata{k}=NaN(num_gate(k),num_dp(k));
    ydata{k}=NaN(num_gate(k),num_dp(k));
    Rho_save{k}=Rho;
    Res_save{k}=Res;

    for i=1:num_dp(k)
        xdata0_tmp = t_center{k}(1:end); 
        ydata0_tmp = Mapp_cell{k}(i,1:end)';
        
        ind_gate_noneg=find(ydata0_tmp>0);
        xdata0_noneg = xdata0_tmp(ind_gate_noneg);
        ydata0_noneg = ydata0_tmp(ind_gate_noneg);

        ind_gate_non_flagged=find(flagIP(i,:)<1);
        xdata0_noflag = xdata0_tmp(ind_gate_non_flagged);
        ydata0_noflag = ydata0_tmp(ind_gate_non_flagged);

        % manually remove some gates or add weights
        % % id_20msec=13; %Krafla remove
        % id_20msec=3; %Krafla remove before 3
        % xdata0 = xdata0_tmp(id_20msec:end);
        % ydata0 = ydata0_tmp(id_20msec:end);
        % % id_50msec=5; %13-17 Krafla
        % id_50msec=5; %3-7 Krafla

        % add weights
        xdata0=log(xdata0_noneg(1:end)); %removing first gates always, not now
        ydata0=log(ydata0_noneg(1:end));
        % weights=ones(numel(xdata0),1);
        % weights(1:id_50msec)=0.1; %Krafla
        % weights(1:id_50msec)=0.5; %Krafla
        % autre option ici au lieu d'enlever les points avec des flags:
        % leur mettre des weights très faibles.
        % Pb : je teste pas encore le SE fit en data processing method.
        % xdata0=xdata0_noflag;
        % ydata0=ydata0_noflag;
        

        if numel(xdata0)>15
            [xData_exp, yData_exp] = prepareCurveData( xdata0, ydata0);
            % W = ones(length(ydata0),1)*30;
            % W(1:4)=10;
            % ft = fittype( 'a*exp(-(x/b)^c)', 'independent', 'x', 'dependent', 'y' );
            ft = fittype( 'a1 - exp(c*x-b1)', 'independent', 'x', 'dependent', 'y' );
            % opts = fitoptions( 'Method', 'NonlinearLeastSquares','Weights', W);
            % opts = fitoptions( 'Method', 'NonlinearLeastSquares','Weights', weights);
            opts = fitoptions( 'Method', 'NonlinearLeastSquares');
            opts.Display = 'Off';
            opts.Robust = 'LAR';
            opts.Lower = [0 -Inf 0.1];
            opts.StartPoint = [1 1 0.5];
            opts.Upper = [10 Inf 1];
            % opts.Lower = [0 50 0];
            % opts.StartPoint = [0.0597795429471558 0.234779913372406 0.353158571222071];
            % opts.StartPoint = [20 100 0.5];
            % opts.Upper = [Inf 5000 1];
            % opts.Weights=W;
            [fitresult_se{k,i}, gof_se{k,i}, output_se{k,i}] = fit( xData_exp, yData_exp, ft, opts );
    
            if gof_se{k, i}.rsquare<0.99
                warning(['Bad SE fit at DP=' num2str(i) ' for file ' file_in{k}(1:end-4) '; try new starting point'])
                opts.StartPoint = [5 5 0.5];
                [fitresult_se{k,i}, gof_se{k,i}, output_se{k,i}] = fit( xData_exp, yData_exp, ft, opts );
                if gof_se{k, i}.rsquare<0.99
                    warning(['(again) Bad SE fit at DP=' num2str(i) ' for file ' file_in{k}(1:end-4) '; try new starting point'])
                    opts.StartPoint = [0.6 0.2 0.7];
                    [fitresult_se{k,i}, gof_se{k,i}, output_se{k,i}] = fit( xData_exp, yData_exp, ft, opts );
                    % start_point=[1 0.4 0.1];
                else
                    % start_point=[10 400 0.8];
                end
            else
                % start_point=[20 100 0.5];
            end
            start_point=opts.StartPoint;
            coeff_alpha(k,i)=exp(fitresult_se{k,i}.a1);
            coeff_tauKWW(k,i)=exp(fitresult_se{k,i}.b1/fitresult_se{k,i}.c);
            coeff_beta(k,i)=fitresult_se{k,i}.c;

            se_fit_mapp{k,i}=coeff_alpha(k,i)*exp(-(xdata0_noneg/coeff_tauKWW(k,i)).^coeff_beta(k,i));   

             % plot
            plot(xdata0_noneg,ydata0_noneg,'o','Color',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',9);
            hold on
            plot (xdata0_noneg,se_fit_mapp{k,i},'-','Color','k','LineWidth',2)
            
            set(gca, 'YScale', 'log')
            set(gca, 'XScale', 'log')

            hXLabel=xlabel('Time (msec)');
            hYLabel=ylabel('Polarizability \eta (mV/V)');
            xlim([1e0,1e4])
            % ylim([1e-2,1e3])
            grid on;

            % try to adjust the fit based on residuals
            % rel_res_thres=0.2;
            res_thres=0.2;
            % W = zeros(length(ydata0),1);
            residuals_Krafla_tmp = output_se{k, i}.residuals;
            % residuals_Krafla_tmp = output_se{k, i}.residuals ./ ydata0;
            % recipvar = abs(1./residuals_Krafla_tmp);
            
            % W=W+recipvar;
            % W(1:4)=10;
            % W(W>60)=60;
            % W(1:6)=max(0,W(1:6)-50);
            % W(1:6)=min(2,W(1:6));
            % weights(abs(residuals_Krafla_tmp)>5)=0.1;
            % opts = fitoptions( 'Method', 'NonlinearLeastSquares','Weights', W);
            % % opts = fitoptions( 'Method', 'NonlinearLeastSquares');
            % % opts.Robust = 'LAR';
            % % opts.StartPoint = start_point;
            % % opts.Lower = [0 50 0];
            % % opts.Upper = [Inf 5000 1];
            points_to_keep=find(abs(residuals_Krafla_tmp)<res_thres);
            opts.Exclude=find(abs(residuals_Krafla_tmp)>res_thres);
            % opts.Weights=W;
            xdata0_excl=xdata0_noneg(points_to_keep);
            ydata0_excl=ydata0_noneg(points_to_keep);
            if numel(xdata0_excl)>15
                [fitresult_se_2{k,i}, gof_se_2{k,i}, output_se_2{k,i}] = fit( xData_exp, yData_exp, ft, opts );
                coeff_alpha_2(k,i)=exp(fitresult_se_2{k,i}.a1);
                coeff_tauKWW_2(k,i)=exp(fitresult_se_2{k,i}.b1/fitresult_se_2{k,i}.c);
                coeff_beta_2(k,i)=fitresult_se_2{k,i}.c;

                se_fit_mapp_2{k,i}=coeff_alpha_2(k,i)*exp(-(xdata0_noneg/coeff_tauKWW_2(k,i)).^coeff_beta_2(k,i)); 

    
                % plot
                plot(xdata0_excl,ydata0_excl,'o','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'MarkerSize',5);
                plot (xdata0_noneg,se_fit_mapp_2{k,i},'-','Color','b','LineWidth',1)
                plot(xdata0_noflag,ydata0_noflag,'*','Color',[1 0 0],'MarkerSize',7);
                txt_title=[file_in{k}(1:end-4) '-DP:' num2str(i) '; I =' num2str(round(1000*current_cell{k,1}(i),0)) ' mA' '; R^2_1 =' num2str(gof_se{k, i}.rsquare,3) '; R^2_2 =' num2str(gof_se_2{k, i}.rsquare,3) ];
                hLegend=legend('All data','Fit 1','Data 2','Fit 2','Data not flagged','Location','Best');
            else
                plot(xdata0_excl,ydata0_excl,'o','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'MarkerSize',5);
                % plot (xdata0,se_fit_mapp_2{k,i},'-','Color','b','LineWidth',1)
                plot(xdata0_noflag,ydata0_noflag,'*','Color',[1 0 0],'MarkerSize',7);
                txt_title=[file_in{k}(1:end-4) '-DP:' num2str(i) '; I =' num2str(round(1000*current_cell{k,1}(i),0)) ' mA' '; R^2_1 =' num2str(gof_se{k, i}.rsquare,3)];
                hLegend=legend('All data','Fit 1','Data 2','Data not flagged','Location','Best');
            end
            hold off

            % txt_title=['DP:' num2str(i) '; I =' num2str(round(1000*current_cell{k,1}(i),0)) ' mA' '; R^2_1 =' num2str(gof_se{k, i}.rsquare,3) '; R^2_2 =' num2str(gof_se_2{k, i}.rsquare,3) ];
            title(txt_title)
            set(gca,'FontName','Times New Roman','FontSize',14);
            set([hXLabel,hYLabel,hLegend],'FontName','Times New Roman','FontSize',14);
            name_out=['png\SE_fit_Krafla\new_fit_log\' file_in{k}(1:end-4) '_DP' num2str(i) '_SEfit_test_residuals.png'];
            exportgraphics(gcf,name_out,'Resolution',300)

            % store final fit results
            
            save_xdata{k,i}=xdata0_tmp;
            save_ydata{k,i}=ydata0_tmp;     
            if ~isempty(gof_se_2{k, i})
                if gof_se{k, i}.rsquare>gof_se_2{k, i}.rsquare
                    coeff_alpha_Krafla(k,i) = coeff_alpha(k,i);
                    coeff_tauKWW_Krafla(k,i) = coeff_tauKWW(k,i);
                    coeff_beta_Krafla(k,i) = coeff_beta(k,i);
                    rsquare_Krafla(k,i) = gof_se{k, i}.rsquare;
                    first_or_second=1;
                else
                    coeff_alpha_Krafla(k,i) = coeff_alpha_2(k,i);
                    coeff_tauKWW_Krafla(k,i) = coeff_tauKWW_2(k,i);
                    coeff_beta_Krafla(k,i) = coeff_beta_2(k,i);
                    rsquare_Krafla(k,i) = gof_se_2{k, i}.rsquare;
                    first_or_second=2;
                end
            end
        else
            warning(['Not enough data points left after flags removed for DP=' num2str(i) 'in file' file_in{k}(1:end-4)])
        end
    end
    % keyboard
end

% save mat file
name_final='mat\Krafla_SEfit_log.mat';
save(name_final,'first_or_second','fitresult_se','gof_se','output_se','se_fit_mapp','fitresult_se_2','gof_se_2','output_se_2','se_fit_mapp_2', 'save_ydata','save_xdata','Mapp_cell','IP_flag_cell','t_center','depths_cell','GateWidth','coeff_alpha_Krafla','coeff_tauKWW_Krafla','coeff_beta_Krafla','rsquare_Krafla','current_cell','Rho_save','Res_save','file_in')