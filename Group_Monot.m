function [Monot_indicator, Thr_sub] = Group_Monot(weight, alpha, xi, rho, Psi, Index_begin, Index_end)
Monot_indicator=0;                                  
%Monot_indicator = 1 represents increasing，
%Monot_indicator = -1 represents decreasing，
%Monot_indicator = 0 represents a threshold

if Index_begin==1 || weight(Index_end)>=weight(Index_begin-1)
    Monot_indicator = 1;
    Thr_sub = Inf;
else
    Thr_sub=(weight(Index_end)*alpha(Index_end)*xi(Index_begin,Index_end)-weight(Index_begin-1)*alpha(Index_begin-1)...
        +weight(Index_begin-1)*alpha(Index_begin-1)*alpha(Index_end)*rho(Index_begin,Index_end))/...
        ((weight(Index_begin-1)-weight(Index_end))*alpha(Index_begin-1)*alpha(Index_end)*xi(Index_begin,Index_end));
    if Thr_sub - Psi(Index_begin) > 10^-12
        Monot_indicator = 0;
    else
        Monot_indicator = -1;
        Thr_sub = Psi(Index_begin);
    end
end
end