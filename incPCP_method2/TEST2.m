p1 = [1; 2];
p2 = [6; 7];
p3 = [4; 5];

p4 = [10; 11];
p5 = [12; 22];
p6 = [6; 7];


A = {};
b = {};
c = {};

b{end+1} = p1;
b{end+1} = p2;
b{end+1} = p3;
c{end+1} = p4;
c{end+1} = p5;
A{end+1} = b;
A{end+1} = c;

pt1 = [6; 7];
pt2 = [8; 9];

test = {};
test{end+1} = pt1;
test{end+1} = pt2;




stat = 0;
for k = 1:size(test,2)
    for i = 1:numel(A)
        for j = 1:size(A{i},2)
            disp(A{i}{j})
            if A{i}{j} == test{k}
                disp('yess')
                stat =1;
                for l = 1:size(test,2)
                    A{i}{end+1} =test{l};
                end                
                break
            end
        end
        if stat == 1
            break
        end
    end
    if stat == 1
        break
    end
end






