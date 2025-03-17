
addpath('~/Desktop/Research/2D code/free space/Yee/');
addpath('~/Desktop/Research/2D code/')

experiment_name = "Example";
mkdir(fullfile(pwd,"Numerical Experiment",experiment_name))

%%
K = 1;
N = [1,2,4,8].*K;

[delx_sol,dely_sol,delt_sol,X_sol,Y_sol,H_sol] = yee_2d_TE(128*K);

delx_full = zeros(size(N));
dely_full = zeros(size(N));
delt_full = zeros(size(N));
X_full = {};
Y_full = {};
H_full = {};

err_full = zeros(size(N));

for i = 1:length(N)
    [delx_full(i),dely_full(i),delt_full(i),...
        X_full{i},Y_full{i},H_full{i}] = yee_2d_TE(N(i));
    
    interp_H = interp2(X_full{i},Y_full{i},H_full{i},X_sol,Y_sol,'spline');
    
    err_full(i) = (norm(interp_H - H_sol))*sqrt(X_sol(end,2) - X_sol(end,1))...
        *sqrt(Y_sol(2,end) - Y_sol(1,end));
end

%%

out_matrix_full = zeros(length(N),5);
out_matrix_full(:,1) = delx_full;
out_matrix_full(:,2) = delt_full;

out_matrix_full(:,3) = err_full;
out_matrix_full(2:end,4) = err_full(1:end-1)./err_full(2:end);
out_matrix_full(2:end,5) = log(out_matrix_full(2:end,4))./log(delx_full(1:end-1)./delx_full(2:end))';
out_matrix_full(1,4) = nan;
out_matrix_full(1,5) = nan;

matrix2latex(out_matrix_full,fullfile(pwd,"Numerical Experiment",...
    experiment_name,'err_table_full.tex'),'columnLabels', ...
    {'$\Delta x = \Delta y$', '$\Delta t$', ...
    '$L_2$ Error', 'Ratio', 'Rate'})
%%


hmm_err(N,experiment_name,delx_sol,dely_sol,delt_sol,X_sol,Y_sol,H_sol,...
    delx_full,err_full,1e-3,1,40)


