function [ Weight ] = Gen_Weight_Rnd( Channel )
global M
global N

Weight = zeros(N,M);
for m=1:M
    Weight_tmp = rand( N,1 );
    Weight(:,m) = Normalization( Weight_tmp );
end

end

