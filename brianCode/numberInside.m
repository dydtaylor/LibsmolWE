function x = numberInside(coords,particleLength,radius)
x = 0 ;

for i = 1:2:particleLength
    dist = sqrt(coords(i).^2 + coords(i+1).^2);
    
    if dist < radius
        
        x = x +1;
        
    end
    
end
