%the output of this function is the integrated version of the DE meal
%function to be used in the T1D version of the model developed for the
%generalized version with Gb+nk, which requires dealing with the
%normalizing constants.

function [del] = DoubleExponentialMeals4(event_times,nutrition,c_gamma,a,b)

% function [del] = DoubleExponentialMeals2(event_times,nutrition,c_gamma,a,b,T)

t = event_times;
nut_time = nutrition(:,1);
n = length(t);

c = (a*b)/(b-a); %normalizing constant (unit: 1/min)
% c = (a*b)/(b*(1-exp(-a*T))-a*(1-exp(-b*T)));

del = zeros(n,1);

for k = 1:n-1
    
    hk = t(k+1)-t(k);
    ind1 = find(nut_time <= t(k));
    t_vec = nut_time(ind1);
    
    if isempty(ind1)
        
        del(k) = 0;
        
    else
        
        sum = 0;
        
        for i = 1:length(t_vec)
            sum = sum + 
                (c * nutrition(ind1(i),2) * (
                    (
                        (exp(a*(t_vec(i)-t(k+1))) - exp((-c_gamma*hk)+a*(t_vec(i)-t(k))))/(c_gamma-a)
                        )-(
                            (exp(b*(t_vec(i)-t(k+1))) - exp((-c_gamma*hk)+b*(t_vec(i)-t(k))))/(c_gamma-b)
                            )
                    ));
        end
        
        del(k) = sum;

    end
    
end

end