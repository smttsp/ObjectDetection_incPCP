function [group_x_left,group_x_right,group_y_top,group_y_bottom] =combine_blobs(S,threshold)

    [B,L,N,A] =  bwboundaries(S,'noholes');
    s = regionprops(L,'centroid');
    centroids = cat(1, s.Centroid);        
    
    D = zeros(size(centroids,1),size(centroids,1));
    groups = {};    
    groupsx = {};
    groupsy = {};
    used = [];
    
    group_x_left = [];
    group_x_right =  [];
    group_y_top = [];
    group_y_bottom =  [];
    
    for i  = 1: size(centroids,1)
        if isempty(find(used == i))               
        group = {};   
        groupx = [];
        groupy = [];
        if i == 6
        end        
        y_bottom =  max(B{i}(:,1));
        y_top = min(B{i}(:,1));
        x_left = min(B{i}(:,2));
        x_right = max(B{i}(:,2));
        status = 0;
        for j  = 1: size(centroids,1)
            if isempty(find(used == j))                        
                status = 1;
                p1 = [centroids(i,1); centroids(i,2)];
                p2 = [centroids(j,1); centroids(j,2)];
               d = norm(p1 - p2); 
                D(i,j) = d;
                if d < threshold                                             
                    
                    group{end+1} = p2;
                    
                    groupx(end+1) = centroids(j,1);
                    groupy(end+1) = centroids(j,2);
                    
                    used(end+1) = j;                                                            
                    if max(B{j}(:,1)) >= y_bottom
                        y_bottom = max(B{j}(:,1));
                    end                    
                    if  min(B{j}(:,1)) <= y_top
                        y_top = min(B{j}(:,1));
                    end                    
                    if min(B{j}(:,2)) <= x_left
                         x_left = min(B{j}(:,2));
                    end                    
                    if max(B{j}(:,2)) >= x_right
                        x_right = max(B{j}(:,2));
                    end                                                            
                end
            end
            
        end
        
        if status == 1
            groups{end+1} = group;
            groupsx{end+1} = groupx;
            groupsy{end+1} = groupy;
            group_x_left(end+1) = x_left;
            group_x_right(end+1) = x_right;
            group_y_top(end+1) = y_top;
            group_y_bottom(end+1) = y_bottom;
        end
        end
    end
      
   
    main_centroidsx = [];
    main_centroidsy = [];

    main_centroidsx_left = [];
    main_centroidsx_right = [];

    main_centroidsy_top = [];
    main_centroidsx_bottom = [];   
    
    

   for  i  = 1:size(groups,2)
       
       if size(groups{i},2) < 1
       groups{i} = [];
       else

       xs = 0;
       ys = 0;
       xs_elements = [];
       ys_elements = [];
       
       num_of_centroids = size(groups{i},2);
       
        for j = 1:size(groups{i},2)
            xs = xs + groups{i}{j}(1);
            ys = ys + groups{i}{j}(2);
            
            xs_elements(end+1) = groups{i}{j}(1);
            ys_elements(end+1) = groups{i}{j}(2);
        end
        
        xs_centroid = round(xs / num_of_centroids); 
        ys_centroid = round(ys / num_of_centroids); 
        
        main_centroidsx_left(end+1) = round(min(xs_elements));
        main_centroidsx_right(end+1) = round(max(xs_elements));
   
       main_centroidsy_top(end+1) = round(min(ys_elements));
       main_centroidsx_bottom(end+1) = round(max(ys_elements));
        
        main_centroidsx(end+1) = xs_centroid;
        main_centroidsy(end+1) = ys_centroid;
       
       end
   end
                 
   
   
   try
for i = 1:size(groupsx,2)        
    for j = 1:size(groupsx{i},2)        
        tf = intersect(groupsx{i},groupsx{j});
        if i~=j & ~isempty(tf)
            groupsx{i} = union(groupsx{i},groupsx{j});
            groupsx{j} = [-1];
            
            if group_x_left(i) > group_x_left(j)
                group_x_left(i) = group_x_left(j);
            else
                group_x_left(i) = group_x_left(i);
            end
            
            group_x_left(j) = [];
            
            if group_x_right(i) < group_x_right(j)
                group_x_right(i) = group_x_right(j);
            else
                group_x_right(i) = group_x_right(i);
            end
            
            group_x_right(j)= [];
            
            if group_y_top(i) > group_y_top(j)
                group_y_top(i) = group_y_top(j);
            else
                group_y_top(i) = group_y_top(i);
            end
            
            group_y_top(j) = [];
            
            if group_y_bottom(i) < group_y_bottom(j)
                group_y_bottom(i) = group_y_bottom(j);
            else
                group_y_bottom(i) = group_y_bottom(i);
            end
            
            group_y_bottom(j) = [];
           
        end
    end
end
   catch
   end


for i = 1:size(groups,2)
       if  isempty(groups{i})
           try
            group_x_left(i) = [];
            group_x_right(i) = [];
            group_y_top(i) = [];
            group_y_bottom(i) = [];
           catch
           end
       end
   end

%  figure,imshow(L)
% hold on       
% plot(centroids(:,1),centroids(:,2),'r*')
% plot(main_centroidsx,main_centroidsy,'b*')
% for i =1:size(group_x_left,2)
% rectangle('Position',[group_x_left(i),group_y_top(i),...
%     group_x_right(i)-group_x_left(i),...
%     group_y_bottom(i) - group_y_top(i)], 'EdgeColor','r','LineWidth',2 );
% 
% end
% hold off
    



