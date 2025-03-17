function [Ex_t,Ey_t] = hmm_E_flux(macro_data,macro_grids,micro_grids,steps_size,material_parameters)
    
    eps0 = 1;
    Ex_x_ind = 0;
    Ex_y_ind = 0;
    eps_inf_mat_Ex = 0;
    Ey_x_ind = 0;
    Ey_y_ind = 0;
    eps_inf_mat_Ey = 0;
    a = 1;

    x = 0;
    x_dual = 0;
    y = 0;
    y_dual = 0;

    Ex_mx_dual = 0;
    Ex_my_dual = 0;
    m_dely = 0;

    Ey_mx = 0;
    Ey_mx_dual = 0;
    Ey_my_dual = 0;
    m_delx = 0;

    unpack(macro_data);
    unpack(macro_grids);
    unpack(micro_grids);
    unpack(steps_size);
    unpack(material_parameters);

    Ex_t = zeros(size(X_dual_Y_main));
    Ey_t = zeros(size(X_main_Y_dual));

    Hz_interpolator = griddedInterpolant(X_dual',Y_dual',H_z','linear');

    J = size(X_dual_Y_main,2);
    parfor i = 2:size(X_dual_Y_main,1)-1
        for j = 1:J
            Ex_t(i,j) = mean(diff(Hz_interpolator({Ex_mx_dual(j,:),Ex_my_dual(i,:)'})',1,1)...
                ./(eps0*eps_inf_mat_Ex(:,:,Ex_x_ind(x_dual(j)),Ex_y_ind(y(i)))),'all')./m_dely;

            
        end
    end
        
    J = size(X_main_Y_dual,2)-1;
    parfor i = 1:size(X_main_Y_dual,1)
        for j = 2:J
            Ey_t(i,j) = mean(diff(Hz_interpolator({Ey_mx_dual(j,:),Ey_my_dual(i,:)'})',1,2)...
                ./(eps0*eps_inf_mat_Ey(:,:,Ey_x_ind(x(j)),Ey_y_ind(y_dual(i)))),'all')./m_delx;
        end
    end


end

function out = epsilon_inf(x)
    
    out = zeros(size(x));

    out(x<=(0.5)) = 2.7;
    out(x>(0.5)) = 1.03;
    out(x==0.5) = mean([2.7,1.03]);
    
end


function unpack(inputs)
    keys = inputs.keys;
    val = inputs.values;
    for i = 1:length(keys)
        assignin('caller',keys{i},val{i});
    end
end


