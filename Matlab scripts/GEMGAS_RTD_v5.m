%% initialize
addpath('mat\')

name_SEfit_GEMGAS='GEMGAS_SEfit_log_saveOK_v2.mat';
load(name_SEfit_GEMGAS,'first_or_second','fitresult_se','gof_se','output_se','se_fit_mapp','fitresult_se_2','gof_se_2','output_se_2','se_fit_mapp_2', 'save_ydata','save_xdata','Mapp_cell','IP_flag_cell','t_center','depths_cell','GateWidth','coeff_alpha_GEMGAS','coeff_tauKWW_GEMGAS','coeff_beta_GEMGAS','rsquare_GEMGAS','current_cell','Rho_save','Res_save','file_in')

%% suite
n=numel(file_in);
col=turbo(n);
num_depth_ref= 1643;
rsquare_fit=NaN(n,num_depth_ref);
rmse_fit=NaN(n,num_depth_ref);
residuals_fit=cell(n,num_depth_ref);
Mapp_std_rel=NaN(n,num_depth_ref);
ydata_rms_tot=NaN(n,num_depth_ref);
misfit_sum=NaN(n,num_depth_ref);
sum_x=NaN(n,num_depth_ref);
weightedTau_x=NaN(n,num_depth_ref);
sum_x_SEfit=NaN(n,num_depth_ref);
weightedTau_x_SEfit=NaN(n,num_depth_ref);
coeff_x=cell(n,num_depth_ref);
resnorms_RTD=NaN(n,num_depth_ref);
residuals_RTD=cell(n,num_depth_ref);

num_dp=size(fitresult_se,2);

first_gate=2; % chosen from CDF
% gof_se_to_use=cell(n,num_dp);
% output_se_to_use=cell(n,num_dp);

%% calculate RTD and plots
for k=1:n
    for i=1:num_dp
        if first_or_second(k,i)==1
            % fitresult_se=fitresult_se;
            gof_se_to_use=gof_se{k,i};
            output_se_to_use=output_se{k,i};
            % se_fit_mapp_to_use{k,i}=se_fit_mapp{k,i};
        elseif first_or_second(k,i)==2
            % fitresult_se=fitresult_se_2;
            gof_se_to_use=gof_se_2{k,i};
            output_se_to_use=output_se_2{k,i};
            % se_fit_mapp_to_use=se_fit_mapp_2;
        end
        % %xdata = real gates
        xdata0=save_xdata{k,i};
        xdata=xdata0(first_gate:end); % on skip la première gate pour avoir l'essentiel des SE fit qui passent la condition residual
        % %xdata = tau interval 1e0 à 1e4 avec 300 points?
        % xdata=logspace(0,4,300);
        ydata=[];
        clf;
        if isempty(fitresult_se{k,i})
            warning(['no SE fit for DP =' num2str(i) 'm for file' file_in{k}])
        else
            rmse_fit(k,i)=gof_se_to_use.rmse;
            rsquare_fit(k,i)=gof_se_to_use.rsquare;   
            residuals_fit{k,i}=output_se_to_use.residuals;
            numobs_fit=output_se_to_use.numobs;
            if rsquare_fit(k,i)<0.99 %0.998
                warning(['Bad SE fit for DP =' num2str(i) 'm for file' file_in{k}])
                % here we dont give value to ydata
            else
                ydata_0=coeff_alpha_GEMGAS(k,i)*exp(-(xdata/coeff_tauKWW_GEMGAS(k,i)).^coeff_beta_GEMGAS(k,i));   
                ydata=ydata_0;

            end
        end
        if isempty(ydata)
            warning(['no data for DP =' num2str(i) 'm for file' file_in{k}])
        else
            % count_file=count_file+1;
            % lgd{2*count_file-1}=suff_date{k};
            % lgd{2*count_file}="";

            % % coeff tau large interval
            % coeff_tau=logspace(0,4,300);
            % % coeff tau short interval
            tau_min=log10(min(xdata));
            tau_max=log10(max(xdata));
            coeff_tau=logspace(tau_min,tau_max,300);

            num_param=numel(coeff_tau);
            num_data=numel(xdata);
            coeff_C=NaN(num_data,num_param);
            for l=1:num_data
                for j=1:num_param
                    coeff_C(l,j)=exp(-xdata(l)/coeff_tau(j));
                end
            end
            [x,resnorm,residual,exitflag,output] = lsqnonneg(coeff_C,ydata);
            resnorms_RTD(k,i)=resnorm; % one value, L2 norm
            residuals_RTD{k,i}=residual; % vector with size = number of gates
            coeff_x{k,i}=x;
            % max_x(k)=max(x);
            sum_x(k,i)=sum(x);
            prod_tmp=x'.*log(coeff_tau);
            weightedTau_x(k,i)=exp(sum(prod_tmp)/sum(x));
            y_test_3 = coeff_C*x;

            % subplot(2,1,1)
            % plot(xdata, ydata,'o','Color',col(k,:),'MarkerFaceColor',col(k,:),'MarkerSize',7)
            % hold on
            % plot(xdata, y_test_3,'-','Color','k','LineWidth',2);
            % 
            % set(gca, 'YScale', 'log')
            % set(gca, 'XScale', 'log')
            % grid on;
            % hXLabel=xlabel('Time (ms)');
            % hYLabel=ylabel('Polarizability \eta (mV/V)');
            % txt=['M_{tot}=' num2str(round(sum(x),0)) 'mV/V'];
            % hText=text(1.5,200,txt);
            % output_title=['File: ' file_in{k}(5:12) ' - DP = ' num2str(i)];
            % hTitle=title(output_title,'Interpreter','none');
            % xlim([1e0,1e4])
            % ylim([1e-2,1e3])
            % yticks_vector=logspace(-1,3,5);
            % yticks(yticks_vector);
            % yticklabels(yticks_vector);
            % 
            % set(gca,'FontName','Times New Roman','FontSize',14);
            % set([hXLabel,hYLabel,hText],'FontName','Times New Roman','FontSize',14);
            % set(hTitle,'FontName','Times New Roman','FontSize',14);
            % hold off;
            % 
            % subplot(2,1,2)
            % semilogx(coeff_tau,x,'o-','Color','k','MarkerFaceColor',col(k,:),'MarkerEdgeColor',col(k,:),'MarkerSize',5)
            % 
            % if k>1
            %     hXLabel=xlabel('Relaxation time distribution (msec)');
            %     set(hXLabel,'FontName','Times New Roman','FontSize',14);
            % end
            % hYLabel2=ylabel('m_k (mV/V)');
            % grid on;
            % 
            % set(gca,'FontName','Times New Roman','FontSize',14);
            % set(hYLabel2,'FontName','Times New Roman','FontSize',14);
            % 
            % fileout=['png\RTD_GEMGAS\new_fit\newfit_RTD_' file_in{k}(1:end-4) '_DP' num2str(i) '.png'];
            % exportgraphics(gcf,fileout,'Resolution',300)
        end
    end
    sprintf([file_in{k} ' done'])
end

% save variables
name_final='mat\GEMGAS_RTD_099_saveXk_fev25_short_interval.mat';
save(name_final,'sum_x','weightedTau_x','coeff_tau', 'coeff_x', 'save_ydata','save_xdata','coeff_alpha_GEMGAS','coeff_tauKWW_GEMGAS','coeff_beta_GEMGAS','rsquare_GEMGAS','file_in','current_cell','Rho_save','Res_save','rsquare_fit','rmse_fit','residuals_fit','resnorms_RTD','residuals_RTD','t_center','first_gate')
