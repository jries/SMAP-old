function indg=brightestinregion(x,y,int,R)
%%

[xs,xind]=sort(x);
ys=y(xind);
ints=int(xind);
[~,iinds]=sort(ints,'descend');

for k=1:length(ints)
    indh=iinds(k);
    if ints(indh)==0
        continue
    end
    xh=xs(indh);
    l=1;
    while xs(indh+l)<xh+R
        d2=(xs(indh+l)-xh).^2+(ys(indh+l)-yh).^2;
        if d2<R^2 
            if ints(indh)>=ints(indh+l)
                ints(indh+l)=0;
            else
                ints(indh)=0;
            end
        end
        l=l+1;
    end
end
