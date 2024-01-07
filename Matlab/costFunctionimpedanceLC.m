%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Written by Juan A. Vasquez-Peralvo
%Luxembourg-Luxembourg
%1/1/2024
%Modified 7/1/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Description%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code allows to calculate the cost function that compares two S21
%parameters and obtain the mean square error to further bein processed by
%an optimizer 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Syntaxys %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input values 
% x is an array of two elements L and C 
% w is the angular velocity of free space
% s21S is the transmission coefficient that will be the reference for the
% optimiztion.
% output values
% cost is the calculated mean square error 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cost = costFunctionimpedanceLC(x, w, s21S)
  
    L = x(1)*1e-12 % here we multiply by the factor to obtain pH
    C = x(2)*1e-15 % here we multiply by the factor to obtain fC
    z0 = 377;

    % This the formula that calculates the S21 parameters based on the L
    % and C values 
    %yFssC = abs(  (  (1j.*w*L./(1j*w*C)) ./  (1j.*w*L +1./(1j*w*C)) ).^(-1) );   
    s21C = 20*log10(abs((2*(  (1j.*w*L./(1j*w*C)) ./  (1j.*w*L +1./(1j*w*C)) ))./ ...
        (z0+2*(  (1j.*w*L./(1j*w*C)) ./  (1j.*w*L +1./(1j*w*C)) ))));
    cost = sum(abs(s21C - s21S).^2); % Compute the cost function

end