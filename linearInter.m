function SInterp=linearInter(t,t_,P)
        if find(t_==t)
            SInterp=P(find(t_==t),:);
        else
            smaller=find(t_<t);
            greater=find(t_>t);
            
            ti=smaller(end);
            tiplus=greater(1);
            
            if abs(ti-t)<abs(tiplus-t)
                SInterp=P(ti,:);
            else
                SInterp=P(tiplus,:);
            end
        end
    end
