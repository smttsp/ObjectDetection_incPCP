function [out] = coloring_cars(S1, curFrame)

    [B,L] = bwboundaries(S1,'noholes');
    %     hold on
    cnt = 0;
    for k = 1:length(B)
        mx=max(B{k});
        mn=min(B{k});
        d=mx-mn;

        a = mn(1):mx(1);
        b = mn(2):mx(2);

        if d(1)*d(2) < 50 || sum(sum(S1(a,b))) < 30
            cnt=cnt+1;
    %             out(a, b, :)=cat(3, zeros(d(1)+1, d(2)+1), zeros(d(1)+1, d(2)+1), zeros(d(1)+1, d(2)+1));
            S1(a,b)=0;
        else
    %             [(mx(1)+mx(1))/2 (mn(2)+mn(2))/2]

    %             idx = kmeans(reshape(S1(a,b), 1, []),2);
    %             [lb, center]=adaptcluster_kmeans(out(a,b,:));
    %             boundary = B{k};
    %             plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
        end
    end
    %     out=imfill(out,'holes');
    %     tic
    %     out=my_kmeans(out,2);
    %     toc

    out= reg_dilate(S1,4);
    [B,L] = bwboundaries(out,'noholes');
    out = label2rgb(L, @jet, [0 0 0]);

    figure(6), imshow(curFrame) 

    rects = zeros(length(B),4);
    cnt = 1;
    for k = 1:length(B)
        mx=max(B{k});
        mn=min(B{k});
        d=mx-mn;

        a = mn(1):mx(1);
        b = mn(2):mx(2);

        if d(1)<15 || d(2)<15 || d(1)*d(2) < 400 || sum(sum(S1(a,b))) < 250
    %             out(a, b, :)=cat(3, zeros(d(1)+1, d(2)+1), zeros(d(1)+1, d(2)+1), zeros(d(1)+1, d(2)+1));
    %             S1(a,b)=0;
        else
            rect =[mn(2) mn(1) mx(2)-mn(2) mx(1)-mn(1)];
            rects(cnt, :)= rect;
            cnt=cnt+1;
        end
    end
    rects = rects(1:cnt-1,:);
    rects=combine_vert_box(rects);
    rects = remove_inner_rect(rects);
    %     figure(4), imshow(out)
    for i = 1:size(rects,1)
        if any(rects(i,:)~=0)
            color=rand(1,3);
            rectangle('Position',rects(i,:),'EdgeColor','r','LineWidth',2 );

        end
        
    end
end

function [new_rect] = combine_vert_box(rects)
    for i = 1:size(rects,1)
        for j = i+1:size(rects,1)
            if all(rects(i,:)==0)||all(rects(j,:)==0)
                continue;
            end
            [res,new_dim]=vert_cont(rects(i,:), rects(j,:));
            
            if res == 1
                rects(i,:) = new_dim;
                rects(j,:) = 0;
            elseif res == 2
                rects(i,:) = new_dim;
                rects(j,:) = new_dim;
            end
        end
    end
    new_rect = rects;

end
    
function [res,new_dim]= vert_cont(r1,r2)
    new_dim = 0;
    perc=0.65;
    max_dist = 5;
    r1_b = r1(2)+r1(4);
    r2_b = r2(2)+r2(4);
    
    r1_r = r1(1)+r1(3);
    r2_r = r2(1)+r2(3);

    res = 0;
    wid = max(r1(3),r2(3));
    if min(r1(3),r2(3))/wid < perc
        return;
    end
    
    if r2(2)-r1_b < max_dist && r1(2)<=r2(2)
        if (min(r1_r,r2_r) - max(r1(1),r2(1)) )/wid > perc
            new_dim = [ min(r1(1),r2(1)) r1(2)  max(r1_r,r2_r) - min(r1(1),r2(1))  r2_b-r1(2)];
            res=1;
        end
    elseif r1(2)-r2_b < max_dist && r2(2)<=r1(2)
        if (min(r1_r,r2_r) - max(r1(1),r2(1)) )/wid > perc
            new_dim = [min(r2(1),r1(1)) r2(2) max(r2_r,r1_r)-min(r2(1),r1(1)) r1_b-r2(2) ];
            res=2;
        end
    end
    
%     if res~=0
%         res
%     end
end



function new_rects = remove_inner_rect(rects)
    for i = 1:size(rects,1)
        for j = i+1:size(rects,1)
            if all(rects(i,:)==0)||all(rects(j,:)==0)
                continue;
            end
            res = intersection(rects(i,:), rects(j,:)); 
            if res == 1
                rects(i,:)=0;
            elseif res == 2
                rects(j,:)=0;
            end
        end
    end
    new_rects=rects;
end

function [res] = intersection(r1,r2)
    thresh = 0.6;
    l = max(r1(1),r2(1));
    r = min(r1(1)+r1(3),r2(1)+r2(3));
    b = max(r1(2),r2(2));
    t = min(r1(2)+r1(4),r2(2)+r2(4));
    r3 = [l r r-l t-b];

    res=0;
    if all(r3>0)
        if ((r3(3)*r3(4))/(r1(3)*r1(4)))>thresh
            res=1;
        elseif ((r3(3)*r3(4))/(r2(3)*r2(4)))>thresh
            res=2;
        end
    end
end

