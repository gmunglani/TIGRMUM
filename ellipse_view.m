function ellipse_view(center,phin,axes)

phi = linspace(0,2*pi,180);
cosphi = cos(phi);
sinphi = sin(phi);
xbar = center(1);
ybar = center(2);

for k = 1:length(phi)
    if (axes(2) > axes(1))
        a = axes(1);
        b = axes(2);
    else
        a = axes(2);
        b = axes(1);
    end
            
    R = [ cos(phin)  sin(phin)
          -sin(phin)   cos(phin)];

    xy = [a*cosphi; b*sinphi];
    xy = R'*xy;

    x = xy(1,:) + xbar;
    y = xy(2,:) + ybar;
    
    plot(y,x,'y','LineWidth',2);
end

%Lines=[(1:size(x,2)); (2:size(x,2)+1)]'; Lines(end,2)=1;
%curve = LineCurvature2D([x y],Lines);