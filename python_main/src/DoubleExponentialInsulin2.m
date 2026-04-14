function [del] = DoubleExponentialInsulin2(event_times,bolus_insulin,c_gamma,c,d,x,z)

% % f1 = @(d) log(d/c)/(d-c) - y;
% % 
% % % Provide an initial guess for b (should not be equal to a)
% % d0 = c + 10;
% % 
% % % Use fsolve (Optimization Toolbox)
% % d = fsolve(f1, d0);

t = event_times;
bolus_ins_time = bolus_insulin(:,1);
n = length(t);

const = 1/((exp(-d*(z-x))-1)/d - (exp(-c*(z-x))-1)/c);

del = zeros(n,1);

% f = @(param, ins_time, delay, next_event_time, integral_time) exp(param*(ins_time+delay-integral_time)-c_gamma*(next_event_time-integral_time));

for k = 1:n-1
    
    hk = t(k+1)-t(k);
    ind1 = find(bolus_ins_time <= t(k));
    t_vec = bolus_ins_time(ind1);

    if isempty(ind1)
        
        del(k) = 0;
        
    else
        
        sum = 0;
        
        for i = 1:length(t_vec)
            if t_vec(i) <= t(k) && t(k) < t_vec(i)+x && t_vec(i) < t(k+1) && t(k+1) <= t_vec(i)+x
                sum = sum;
            elseif t_vec(i)+x <= t(k) && t(k) < t_vec(i)+z && t_vec(i)+x < t(k+1) && t(k+1) <= t_vec(i)+z
                sum = sum + (const * bolus_insulin(ind1(i),2) * (((exp(c*(t_vec(i)+x-t(k+1)))-exp((-c_gamma*hk)+c*(t_vec(i)+x-t(k))))/(c_gamma-c))-((exp(d*(t_vec(i)+x-t(k+1)))-exp((-c_gamma*hk)+d*(t_vec(i)+x-t(k))))/(c_gamma-d))));
            elseif t_vec(i)+z <= t(k)
                sum = sum;
            else
                warning('There should not be another case!')
            end
        end

        del(k) = sum;

    end
    
end

end


%(((exp(c*(t_vec(i)+x-tend))-exp((-c_gamma*(tend-tstart))+c*(t_vec(i)+x-tstart)))/(c_gamma-c))-((exp(d*(t_vec(i)+x-tend))-exp((-c_gamma*(tend-tstart))+d*(t_vec(i)+x-tstart)))/(c_gamma-d))));

% % %             if t_vec(i) <= t(k) && t(k) < t_vec(i)+x && t_vec(i) <= t(k+1) && t(k+1) < t_vec(i)+x
% % %                 tstart = 0;
% % %                 tend = 0;
% % %             elseif t_vec(i) <= t(k) && t(k) < t_vec(i)+x && t_vec(i)+x <= t(k+1) && t(k+1) < t_vec(i)+z
% % %                 tstart = t_vec(i) + x;
% % %                 tend = t(k+1);
% % %             elseif t_vec(i) <= t(k) && t(k) < t_vec(i)+x && t_vec(i)+z <= t(k+1)
% % %                 tstart = t_vec(i) + x;
% % %                 tend = t_vec(i) + z;
% % %             elseif t_vec(i) + x <= t(k) && t(k) < t_vec(i)+z && t_vec(i) + x <= t(k+1) && t(k+1) < t_vec(i)+z
% % %                 tstart = t(k);
% % %                 tend = t(k+1);
% % %             elseif t_vec(i) + x <= t(k) && t(k) < t_vec(i)+z && t_vec(i)+z <= t(k+1)
% % %                 tstart = t(k);
% % %                 tend = t_vec(i) + z;
% % %             else
% % %                 tstart = 0;
% % %                 tend = 0;
% % %             end