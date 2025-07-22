function [ Distance ] = Rand_Circle( N,M,Radius_max,Radius_min)       
i=1;
Distance=zeros(N,M);                                 
Coordinate_temporary=zeros(1,2);
while i<=M*N
    Coordinate_temporary=Radius_max*(rand(2,1)*2-1);
    if norm(Coordinate_temporary, 2)<Radius_max && norm(Coordinate_temporary, 2)>Radius_min
       n = ceil(i/M);
       m = i-(n-1)*M;
       Distance(n,m) = norm(Coordinate_temporary, 2);
       i=i+1;
    end        
end    
end

