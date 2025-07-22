function [Conv_indicator, Merge_index] = Conv_check(weight, alpha, xi, rho, Psi, Monot_list, Thr, Group_begin)
global N;
Conv_indicator = 1;
Merge_index =-1;
for i=2:length(Group_begin)
    if Monot_list(i)==1
        Conv_indicator = 0;
        Merge_index = i;
        break;
    end
    if Monot_list(i) == -1 && i<length(Group_begin)
        Conv_indicator = 0;
        Merge_index = i+1;
        break;
    end
    if i<length(Group_begin) && xi(Group_begin(i),Group_begin(i+1))*Thr(i)-rho(Group_begin(i),Group_begin(i+1))<Thr(i+1)
        Conv_indicator = 0;
        Merge_index = i+1;
        break;
    end
end