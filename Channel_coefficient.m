function [ Channel,Distance ] = Channel_coefficient( N,M,Radius_max,Radius_min,Bandwidth )    
    %Function to generate the channel
    Noise_dB=-144+10*log10(Bandwidth);                                           
    Variance_Noise=10^(Noise_dB/10);
    Channel=sqrt(0.5)*randn(N,M)+sqrt(0.5)*j*randn(N,M);
    Distance=Rand_Circle( N,M,Radius_max,Radius_min);
    L=-128-37.6*log10(Distance/1000);

    Channel = abs(Channel);
    Channel = Channel.*Channel;
    Channel = Channel.*(10.^(L/10));
    Channel = Channel/Variance_Noise;
    
    for m=1:M        
        Channel(:,m) = sort( Channel(:,m) );
        Distance(:,m) = sort( Distance(:,m), 'descend');
    end
end

