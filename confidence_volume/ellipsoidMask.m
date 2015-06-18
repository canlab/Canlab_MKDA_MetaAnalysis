function vox = ellipsoidMask(a,b,c)

xP = [];
yP = [];
zP = [];

for z = -c:c
    
    ymax = sqrt(b.^2 * (1 - z.^2 ./ c.^2));
    
    for y = -ymax:ymax
        
            xmax = sqrt(a ^ 2 * (1 - z^2 / c^2 - y^2 / b^2));
            
            for x = -xmax:xmax
                
                xP = [xP;x];
                yP = [yP;y];
                zP = [zP;z];
                
            end
            
    end
    
end

vox = [xP yP zP];

return