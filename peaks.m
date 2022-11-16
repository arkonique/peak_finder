function [peak,k] = peaks(rates,times,ts,te,alimit,slimit,aevtime,sevtime,background,j)
    k=j+1;
    r=rates(times>=ts & times<=te);
    t=times(times>=ts & times<=te);
    peak=[];
    [c,i]=max(r);
    tp=t(i);
    %plot(t/3.154e7,log10(r));
    %saveas(gcf,['tests/part_' num2str(k) '.png']);
    if (numel(c)==0)
        c=background;
    end
    if (c>=slimit)
        evtime = sevtime;
    elseif (c>=alimit && c<slimit) 
        evtime=aevtime;
    end
    if (c>=alimit)
        if (i==1)
            r2=r(2:end);
            t2=t(2:end);
            [peak,k]=peaks(r2,t2,ts,te,alimit,slimit,aevtime,sevtime,background,k);
        elseif (i==numel(r))
            r2=r(1:end-1);
            t2=t(1:end-1);
            [peak,k]=peaks(r2,t2,ts,te,alimit,slimit,aevtime,sevtime,background,k);
        else
            if (c>=1e-3)
                l=1;
                
            else
                l=0;
            end
            peak=[peak;tp,c,l];
            % fprintf("Current peak: %f,%f,%d\n",tp/3.154e7,log10(c),l);
            t1=t(t>=tp+evtime);
            r1=r(t>=tp+evtime);
            t2=t(t<=tp-evtime);
            r2=r(t<=tp-evtime);
            [p1,k]=peaks(r1,t1,ts,te,alimit,slimit,aevtime,sevtime,background,k);
            [p2,k]=peaks(r2,t2,ts,te,alimit,slimit,aevtime,sevtime,background,k);
            peak=[peak;p1;p2];

        end
    else
        peak=[];
    end

end