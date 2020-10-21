% Sahan Kanchnana Madduma Gamage 
% Undergrad @ Monash University
% MATLAB Code that creates a vector of sugar consumption with time steps of
% 1 minute for a chosen number of meals at a chosen consumption rate.
clear all;
nomeals = floor(rand*3)+1;
sugarintake=(floor((rand*150))+50)*10^-3;
sugarintakerate=(rand*5)*10^-3;
a =  round(sugarintake/sugarintakerate);
b = rand(nomeals,1);
b  = b/sum(b)*a;
out=round(b);
out(1) = out(1) - (a-sum(out));
possiblepositions=zeros(nomeals);
okay=0;
error=0;
unacceptable=0;
while (okay==0)
    startimes=floor(rand(length(out),1)*1440);
    startimes=sort(startimes);
    for i2=1:length(out)
        endtimes=startimes+out(i2);
        if((sum(endtimes>1440))>=nomeals)
           unacceptable=1;
           break
        else
            possiblepositions(:,i2)=[~(endtimes(1:length(endtimes)-1)>=startimes(2:length(endtimes))).*(~(endtimes(1:length(endtimes)-1)>1440));~(endtimes(length(endtimes))>1440)]';
        end
    end
    if(unacceptable==1)
        unacceptable=0;
        continue
    else
        concatonated=[out';sum(possiblepositions);possiblepositions];
        [temp, order] = sort(concatonated(2,:));
        concatonated = concatonated(:,order);
        out=concatonated(1,:)';
        possiblepositions=concatonated(3:end,:);
        storage=zeros(length(possiblepositions),1);
        for i3=1:length(possiblepositions)
            if(isempty(min(find(possiblepositions(:,i3)))))
                error=1;
                break;
            end
            storage(min(find(possiblepositions(:,i3))))=i3;
            possiblepositions(min(find(possiblepositions(:,i3))),:)=0;
        end
        if(error==1)
            unacceptable=0;
            error=0;
        else
            okay=1;
        end
    end
end
food_vec=zeros(1,1440);
for i=1:length(storage)
    food_vec(startimes(i):startimes(i)+out(storage(i)))=1;
end
food_vec=[0,food_vec]*sugarintakerate;
plot(food_vec);