function [res]=mu(k)
    absk=abs(k);
    if absk<1/3
        res=0;
    elseif absk<2/3 && 1/3<=absk
        res=3*(absk-1/3);
    elseif absk<4/3 && 2/3 <= absk
        res=1 + 3/2*(2/3-absk);
    else
        res=0;
    end
end