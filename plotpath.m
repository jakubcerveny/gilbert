# Octave helper function to plot a 2D or 3D colored curve

function h = plotpath(P)

    x = P(:,1)';
    y = P(:,2)';

    if (size(P,2) >= 3)
        z = P(:,3)';
    else
        z = zeros(size(x));
    endif

    col = 1:size(x,2);

    h = surface([x;x],[y;y],[z;z],[col;col],...
                'facecolor','none',...
                'edgecolor','interp',...
                'linewidth',2);
