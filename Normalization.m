function [ Weight ] = Normalization( Weight )     
global N
   Weight=Weight*N/sum(Weight);
end

