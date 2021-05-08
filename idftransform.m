function [X]=idftransform(values, t, f)
    k=0;
    X=zeros(1,length(t));
    for ti = t % time range
        k = k+1;
        X(k) = trapz(f, values.*exp(2i*pi*ti*f)); % integrate 
    end
end