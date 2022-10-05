function [estmean] = emergence(x,y,num)

len = 100000;

maxTime = 24*60*9;
warmupTime = 24*60*2;
firstArrivalAfterWarmUp = 0;

ambarr = ones(1,num);
norarr = ones(1,num);
ent = ones(1,num);
rec = ones(1,num);
lab = ones(1,num);
exam = ones(1,num);
reexam = ones(1,num);
treat = ones(1,num);
emer = ones(1,num);

ambarrtime = exprnd(y(1),num,len);
norarrtime = exprnd(y(2),num,len);
enttime = exprnd(1,num,len);
rectime = exprnd(5,num,len);
labtime = exprnd(30,num,len);
examtime = exprnd(5,num,len);
reexamtime = exprnd(3,num,len);
treattime = exprnd(60,num,len);
emertime = exprnd(120,num,len);

for i=1:num

    entqueue = zeros(3,1);
    recqueue = zeros(3,1);
    examqueue = zeros(4,1);
    labqueue = zeros(3,1);
    treatqueue = zeros(3,1);
    emerqueue = zeros(3,1);

    Events = zeros(4,2);

    next = 0;

    Events(:,1) = [norarrtime(i,norarr(i));1;1;1];
    Events(:,2) = [ambarrtime(i,ambarr(i));2;2;2];
    norarr(i) = norarr(i) + 1;
    ambarr(i) = ambarr(i) + 1;

    nArrivals = 2;

    if (next >=warmupTime && firstArrivalAfterWarmUp ==0)
        firstArrivalAfterWarmUp = 1;
    elseif (Events(1,2) >= warmupTime && firstArrivalAfterWarmUp ==0)
        firstArrivalAfterWarmUp = 2;
    end

    [~,nextEvent] = min(Events(1,:));
    Time = Events(1,nextEvent);
    Type = Events(2,nextEvent);
    Attr = Events(3,nextEvent);
    Customer = Events(4,nextEvent);

    waitingTimes = zeros(2,2);
    nExits=0;

    while (Time <= maxTime)

        if Type ==1

            if (x(1) >0)
                x(1)=x(1)-1;
                sTime = enttime(i,ent(i));
                ent(i) = ent(i)+1;
                u = unifrnd(0,1,1);
                if (u <=0.2)
                    Events = [Events [(Time+sTime);3;100;Customer]];
                else
                    Events = [Events [(Time+sTime);4;Attr;Customer]];
                end
            else
                entqueue = [entqueue [Customer;Time;Attr]];
            end

            nArrivals = nArrivals +1;
            next = Time + norarrtime(i,norarr(i));
            norarr(i) = norarr(i) + 1;
            Events = [Events [next;Type;Attr;nArrivals]];
            waitingTimes = [waitingTimes zeros(2,1)];

            if(next >= warmupTime && firstArrivalAfterWarmUp ==0)
                firstArrivalAfterWarmUp = nArrivals;
            end

        elseif Type==2

            if (x(2) > 0)
                x(2)=x(2)-1;
                sTime = rectime(i,rec(i));
                rec(i)=rec(i)+1;
                u = unifrnd(0,1,1);
                if (u <=0.2)
                    Events = [Events [(Time+sTime);3;1000;Customer]];
                elseif (u <= 0.6)
                    Events = [Events [(Time+sTime);5;Attr;Customer]];
                else
                    Events = [Events [(Time+sTime);6;Attr;Customer]];
                end
            else
                recqueue = [recqueue [Customer;Time;Attr]];
            end

            nArrivals = nArrivals +1;
            next = Time + ambarrtime(i,ambarr(i));
            ambarr(i) = ambarr(i)+1;
            Events = [Events [next;Type;Attr;nArrivals]];
            waitingTimes = [waitingTimes zeros(2,1)];

            if(next >= warmupTime && firstArrivalAfterWarmUp ==0)
                firstArrivalAfterWarmUp = nArrivals;
            end

        elseif Type==3
            waitingTimes(2,Customer) = 4;
            nExits = nExits +1;

            if Attr == 100
                if (length(entqueue(1,:)) > 1)
                    waitingTimes(1,entqueue(1,2)) = waitingTimes(1,entqueue(1,2))+(Time - entqueue(2,2));
                    sTime = enttime(i,ent(i));
                    ent(i)=ent(i)+1;
                    u = unifrnd(0,1,1);
                    if (u <=0.2)
                        Events = [Events [(Time+sTime);3;100;entqueue(1,2)]];
                    else
                        Events = [Events [(Time+sTime);4;entqueue(3,2);entqueue(1,2)]];
                    end
                    entqueue(:,2)=[];
                else
                    x(1) = x(1) + 1;
                end
            else
                if (length(recqueue(1,:)) > 1)
                    waitingTimes(1,recqueue(1,2)) = waitingTimes(1,recqueue(1,2))+(Time - recqueue(2,2));
                    sTime = rectime(i,rec(i));
                    rec(i)=rec(i)+1;
                    u = unifrnd(0,1,1);
                    if (u <=0.2)
                        Events = [Events [(Time+sTime);3;1000;Customer]];
                    elseif recqueue(3,2)==2 && (u <= 0.6)
                        Events = [Events [(Time+sTime);5;recqueue(3,2);recqueue(1,2)]];
                    else
                        Events = [Events [(Time+sTime);6;recqueue(3,2);recqueue(1,2)]];
                    end
                    recqueue(:,2)=[];
                else
                    x(2) = x(2) + 1;
                end
            end

        elseif Type==4
            if (length(entqueue(1,:)) > 1)
                waitingTimes(1,entqueue(1,2)) = waitingTimes(1,entqueue(1,2))+(Time - entqueue(2,2));
                sTime = enttime(i,ent(i));
                ent(i)=ent(i)+1;
                u = unifrnd(0,1,1);
                if (u <=0.2)
                    Events = [Events [(Time+sTime);3;100;entqueue(1,2)]];
                else
                    Events = [Events [(Time+sTime);4;entqueue(3,2);entqueue(1,2)]];
                end
                entqueue(:,2)=[];
            else
                x(1) = x(1) + 1;
            end

            if (x(2) > 0)
                x(2) = x(2)-1;
                sTime = rectime(i,rec(i));
                rec(i)=rec(i)+1;
                Events = [Events [(Time+sTime);6;Attr;Customer]];
            else
                recqueue = [recqueue [Customer;Time;Attr]];
            end

        elseif Type==6
            if (length(recqueue(1,:)) > 1)
                waitingTimes(1,recqueue(1,2)) = waitingTimes(1,recqueue(1,2))+(Time - recqueue(2,2));
                sTime = rectime(i,rec(i));
                rec(i)=rec(i)+1;
                u = unifrnd(0,1,1);
                if (u <=0.2)
                    Events = [Events [(Time+sTime);3;1000;Customer]];
                elseif recqueue(3,2)==2 && (u <= 0.6)
                    Events = [Events [(Time+sTime);5;recqueue(3,2);recqueue(1,2)]];
                else
                    Events = [Events [(Time+sTime);6;recqueue(3,2);recqueue(1,2)]];
                end
                recqueue(:,2)=[];
            else
                x(2) = x(2)+1;
            end

            if (x(3) > 0)
                x(3) = x(3)-1;
                sTime = examtime(i,exam(i));
                exam(i)=exam(i)+1;
                u = unifrnd(0,1,1);
                if (u <= 0.5)
                    Events = [Events [(Time+sTime);7;Attr;Customer]];
                elseif Attr == 2 ||  (u <= 0.86)
                    Events = [Events [(Time+sTime);8;Attr;Customer]];
                elseif (u <= 0.9)
                    Events = [Events [(Time+sTime);5;Attr;Customer]];
                else
                    Events = [Events [(Time+sTime);9;Attr;Customer]];
                end
            else
                examqueue = [examqueue [Customer;Time;Attr;0]];
            end

        elseif Type == 7
            if (length(examqueue(1,:)) > 1)
                waitingTimes(1,examqueue(1,2)) = waitingTimes(1,examqueue(1,2))+(Time - examqueue(2,2));

                if (examqueue(3,2)==1)
                    sTime = reexamtime(i,reexam(i));
                    reexam(i) = reexam(i)+1;
                else
                    sTime = examtime(i,exam(i));
                    exam(i) = exam(i)+1;
                end

                u = unifrnd(0,1,1);
                if (u <= 0.5)
                    Events = [Events [(Time+sTime);7;Attr;Customer]];
                elseif Attr == 2 ||  (u <= 0.86)
                    Events = [Events [(Time+sTime);8;Attr;Customer]];
                elseif (u <= 0.9)
                    Events = [Events [(Time+sTime);5;Attr;Customer]];
                else
                    Events = [Events [(Time+sTime);9;Attr;Customer]];
                end
                examqueue(:,2) = [];
            else
                x(3) = x(3)+1;
            end

            if (x(4) > 0)
                x(4) = x(4)-1;
                sTime = labtime(i,lab(i));
                lab(i) = lab(i)+1;
                Events = [Events [(Time+sTime);10;Attr;Customer]];
            else
                labqueue = [labqueue [Customer;Time;Attr]];
            end

        elseif Type == 10

            if(length(labqueue(1,:))>1)
                sTime = labtime(i,lab(i));
                lab(i) = lab(i)+1;
                Events = [Events [(Time+sTime);10;labqueue(3,2);labqueue(1,2)]];
                labqueue(:,2)=[];
            else
                x(4) = x(4)+1;
            end

            if(x(3)>0)
                x(3)=x(3)-1;
                sTime = reexamtime(i,reexam(i));
                reexam(i) = reexam(i)+1;
                u = unifrnd(0,1,1);
                if (u <= 0.5)
                    Events = [Events [(Time+sTime);7;Attr;Customer]];
                elseif Attr == 2 ||  (u <= 0.86)
                    Events = [Events [(Time+sTime);8;Attr;Customer]];
                elseif (u <= 0.9)
                    Events = [Events [(Time+sTime);5;Attr;Customer]];
                else
                    Events = [Events [(Time+sTime);9;Attr;Customer]];
                end
            else
                examqueue = [examqueue [Customer;Time;Attr;1]];
            end

        elseif Type == 8
            if(length(examqueue(1,:))>1)
                waitingTimes(1,examqueue(1,2)) = waitingTimes(1,examqueue(1,2)) + (Time - examqueue(2,2));

                if(examqueue(3,2)==1)
                    sTime = reexamtime(i,reexam(i));
                    reexam(i) = reexam(i)+1;
                else
                    sTime = examtime(i,exam(i));
                    exam(i) = exam(i)+1;
                end

                u = unifrnd(0,1,1);
                if (u <= 0.5)
                    Events = [Events [(Time+sTime);7;Attr;Customer]];
                elseif Attr == 2 ||  (u <= 0.86)
                    Events = [Events [(Time+sTime);8;Attr;Customer]];
                elseif (u <= 0.9)
                    Events = [Events [(Time+sTime);5;Attr;Customer]];
                else
                    Events = [Events [(Time+sTime);9;Attr;Customer]];
                end
                examqueue(:,2) = [];
            else
                x(3) = x(3)+1;
            end

            if (x(5) > 0)
                x(5) = x(5) -1;
                sTime = treattime(i,treat(i));
                treat(i) = treat(i) +1;
                Events = [Events [(Time+sTime);11;Attr;Customer]];
            else
                treatqueue = [treatqueue [Customer;Time;Attr]];
            end

        elseif Type == 5

            if Attr == 2
                if (length(recqueue(1,:)) > 1)
                    waitingTimes(1,recqueue(1,2)) = waitingTimes(1,recqueue(1,2))+(Time - recqueue(2,2));
                    sTime = rectime(i,rec(i));
                    rec(i)=rec(i)+1;
                    u = unifrnd(0,1,1);
                    if (u <=0.2)
                        Events = [Events [(Time+sTime);3;1000;Customer]];
                    elseif recqueue(3,2)==2 && (u <= 0.6)
                        Events = [Events [(Time+sTime);5;recqueue(3,2);recqueue(1,2)]];
                    else
                        Events = [Events [(Time+sTime);6;recqueue(3,2);recqueue(1,2)]];
                    end
                    recqueue(:,2)=[];
                else
                    x(2) = x(2) + 1;
                end
            else
                if(length(examqueue(1,:))>1)
                    waitingTimes(1,examqueue(1,2)) = waitingTimes(1,examqueue(1,2)) + (Time - examqueue(2,2));

                    if(examqueue(3,2)==1)
                        sTime = reexamtime(i,reexam(i));
                        reexam(i) = reexam(i)+1;
                    else
                        sTime = examtime(i,exam(i));
                        exam(i) = exam(i)+1;
                    end

                    u = unifrnd(0,1,1);
                    if (u <= 0.5)
                        Events = [Events [(Time+sTime);7;Attr;Customer]];
                    elseif Attr == 2 ||  (u <= 0.86)
                        Events = [Events [(Time+sTime);8;Attr;Customer]];
                    elseif (u <= 0.9)
                        Events = [Events [(Time+sTime);5;Attr;Customer]];
                    else
                        Events = [Events [(Time+sTime);9;Attr;Customer]];
                    end
                    examqueue(:,2) = [];
                else
                    x(3) = x(3)+1;
                end
            end

            if (x(6) > 0)
                x(6) = x(6) -1;
                sTime = emertime(i,emer(i));
                emer(i) = emer(i)+1;
                Events = [Events [(Time+sTime);12;Attr;Customer]];
            else
                emerqueue = [emerqueue [Customer;Time;Attr]];
            end

        elseif Type==9

            waitingTimes(2,Customer) = 3;
            nExits = nExits +1;

            if(length(examqueue(1,:))>1)
                waitingTimes(1,examqueue(1,2)) = waitingTimes(1,examqueue(1,2)) + (Time - examqueue(2,2));

                if(examqueue(3,2)==1)
                    sTime = reexamtime(i,reexam(i));
                    reexam(i) = reexam(i)+1;
                else
                    sTime = examtime(i,exam(i));
                    exam(i) = exam(i)+1;
                end

                u = unifrnd(0,1,1);
                if (u <= 0.5)
                    Events = [Events [(Time+sTime);7;Attr;Customer]];
                elseif Attr == 2 ||  (u <= 0.86)
                    Events = [Events [(Time+sTime);8;Attr;Customer]];
                elseif (u <= 0.9)
                    Events = [Events [(Time+sTime);5;Attr;Customer]];
                else
                    Events = [Events [(Time+sTime);9;Attr;Customer]];
                end
                examqueue(:,2) = [];
            else
                x(3) = x(3)+1;
            end

        elseif Type == 11
            nExits = nExits +1;
            waitingTimes(2,Customer) = 2;

            if (length(treatqueue(1,:))>1)
                waitingTimes(1,treatqueue(1,2)) = waitingTimes(1,treatqueue(1,2)) + (Time - treatqueue(2,2));
                sTime = treattime(i,treat(i));
                treat(i) = treat(i) +1;
                Events = [Events [(Time+sTime);11;treatqueue(3,2);treatqueue(1,2)]];
                treatqueue(:,2)=[];
            else
                x(5) = x(5)+1;
            end

        elseif Type == 12
            nExits = nExits + 1;
            waitingTimes(2,Customer) = 1;

            if (length(emerqueue(1,:))>1)
                waitingTimes(1,emerqueue(1,2)) = waitingTimes(1,emerqueue(1,2)) + (Time - emerqueue(2,2));
                sTime = emertime(i,emer(i));
                emer(i) = emer(i)+1;
                Events = [Events [(Time+sTime);12;emerqueue(3,2);emerqueue(1,2)]];
                emerqueue(:,2)=[];
            else
                x(6)=x(6)+1;
            end
        end

        Events(:,nextEvent) = [];
        [~,nextEvent] = min(Events(1,:));
        Time = Events(1,nextEvent);
        Type = Events(2,nextEvent);
        Attr = Events(3,nextEvent);
        Customer = Events(4,nextEvent);
    end

    totalWT = zeros(4,1);
    nCustomers = zeros(4,1);
    for j = firstArrivalAfterWarmUp:length(waitingTimes)
        typ = waitingTimes(2,j);
        if typ~=0
            totalWT(typ) = totalWT(typ) + waitingTimes(1,j);
            nCustomers(typ) =  nCustomers(typ)+1;
        end
    end
    avgWT = totalWT./nCustomers;
    throughput = nExits/(maxTime/60);
end
estmean = avgWT(1,:);
end