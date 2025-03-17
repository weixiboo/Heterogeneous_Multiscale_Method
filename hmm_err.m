
function hmm_err(N,experiment_name,delx_sol,dely_sol,delt_sol,X_sol,Y_sol,H_sol,...
    delx_full,err_full,rho_constant,rho_div,m_Nx)

    if rho_constant == 0
        my_filepath = strcat("Numerical Experiment/",experiment_name,...
            '/rho_shrinking/',strcat("mNx_",num2str(m_Nx)),...
            strcat("/rho_delx_div_",num2str(rho_div)));
    else
        my_filepath = strcat("Numerical Experiment/",experiment_name,...
            '/rho_constant/',strcat("mNx_",num2str(m_Nx)),...
            strcat("/rho_",num2str(rho_constant/rho_div, '%.3e')));
    end

    mkdir(fullfile(pwd,my_filepath))

    delx_hmm = zeros(size(N));
    dely_hmm = zeros(size(N));
    delt_hmm = zeros(size(N));
    X_hmm = {};
    Y_hmm = {};
    H_hmm = {};
    
    guide_line_1 = 1./N;
    guide_line_2 = 1./(N.^2);
    
    err_hmm = zeros(size(N));
    
    
    for i = 1:length(N)
        [delx_hmm(i),dely_hmm(i),delt_hmm(i),...
            m_delx_hmm(i),m_dely_hmm(i),...
            X_hmm{i},Y_hmm{i},H_hmm{i}] = hmm_first_order_2d(N(i),rho_constant,rho_div,m_Nx);
    
        interp_H = interp2(X_hmm{i},Y_hmm{i},H_hmm{i},X_sol,Y_sol,'spline');
        
        err_hmm(i) = (norm(interp_H - H_sol))*sqrt(delx_sol)...
            *sqrt(dely_sol);
    end
    
    %%
    line_diff_1 = log(err_hmm(1)) - log(guide_line_1(1))+ 0.5;
    line_diff_2 = log(err_hmm(1)) - log(guide_line_2(1))- 0.5;
    
    figure(6)
    plot(log(delx_hmm),log(err_hmm),'-o')
    hold on
    plot(log(delx_full),log(err_full),'-*')
    plot(log(delx_hmm),log(guide_line_1) + line_diff_1,':',Color="k")
    plot(log(delx_hmm),log(guide_line_2) + line_diff_2,'--',Color="k")
    title('Log-Log Plot of Final Snapshot Error')
    xlabel('log(\Delta x)')
    ylabel('log(L_2 Error)')
    legend('HMM Error','Full Res Error','O(h) guideline','O(h^2) guideline','Location','southeast')
    hold off
    disp("HMM Final Snapshot order of convergence")
    disp((log(err_hmm(end)/err_hmm(end-1)))/(log(delx_hmm(end)/delx_hmm(end-1))))
    disp("Full Final Snapshot order of convergence")
    disp((log(err_full(end)/err_full(end-1)))/(log(delx_full(end)/delx_full(end-1))))
    saveas(figure(6),fullfile(pwd,my_filepath,'err_plot_hmm.fig'));
    saveas(figure(6),fullfile(pwd,my_filepath,'err_plot_hmm.png'));
    
    %%
    
    if any(delx_hmm ~= dely_hmm)
        warning('Delx != Dely')
    end
    
    out_matrix = zeros(length(N),6);
    out_matrix(:,1) = delx_hmm;
    out_matrix(:,2) = delt_hmm;
    
    out_matrix(:,3) = m_delx_hmm;
    
    out_matrix(:,4) = err_hmm;
    out_matrix(2:end,5) = err_hmm(1:end-1)./err_hmm(2:end);
    out_matrix(2:end,6) = log(out_matrix(2:end,5))./log(delx_hmm(1:end-1)./delx_hmm(2:end))';
    out_matrix(1,5) = nan;
    out_matrix(1,6) = nan;
    

    if rho_constant == 0
        matrix2latex(out_matrix,fullfile(pwd,my_filepath,'err_table_hmm.tex'),'columnLabels', ...
            {'$\Delta x = \Delta y$', '$\Delta t$', ...
            '$\delta x = \delta y$',...
            '$L_2$ Error', 'Ratio', 'Rate'},...
            'caption',sprintf("$\\rho = \\frac{\\Delta x}{%i}$,$N^{\\text{micro}}_x = N^{\\text{micro}}_y = %i$",rho_div,m_Nx))
    else
        matrix2latex(out_matrix,fullfile(pwd,my_filepath,'err_table_hmm.tex'),'columnLabels', ...
            {'$\Delta x = \Delta y$', '$\Delta t$', ...
            '$\delta x = \delta y$',...
            '$L_2$ Error', 'Ratio', 'Rate'},...
            'caption',sprintf("$\\rho = %s $,$ N^{\\text{micro}}_x = N^{\\text{micro}}_y = %i $",my_num2str(rho_constant/rho_div),m_Nx))
    end


end

function latexStr = my_num2str(number)
    if number >= 1e4 || number < 1e-2
        scientificStr = num2str(number, '%.3e'); % Example: '1.234568e+05'
        latexStr = strrep(scientificStr, 'e+0', ' \times 10^{');
        latexStr = strrep(latexStr, 'e-0', ' \times 10^{-');
        
        latexStr = strrep(latexStr, 'e+', ' \times 10^{');
        latexStr = strrep(latexStr, 'e-', ' \times 10^{-');
    
        latexStr = strcat(latexStr,'}');
    else
        number = round(number,3);
        latexStr = num2str(number);
    end

end
