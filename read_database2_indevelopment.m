
% A is numbers, B is text

[A,B] = xlsread(dbname);

names = B(1,:);
B = B(2:end,:);

for i = 1:length(names)
    
    
    