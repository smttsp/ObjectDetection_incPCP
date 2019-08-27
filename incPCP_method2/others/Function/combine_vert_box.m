

function [new_rect] = combine_vert_box(rects)
rects = remove_inner_rect(rects);
% rects = remove_small_rects(rects);
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

function new_rects = remove_small_rects(rects)
    th = 10;
    for i = 1:size(rects,1)
        if rects(i,3) < th || rects(i,4) < th
            rects(i,:) = 0;
        end
    end

    new_rects = rects;
end

