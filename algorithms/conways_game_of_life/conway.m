function conway(grid_size, time_steps, config)

switch config
    case 1 % random initial distribution
        grid_mat = randi([0,1], grid_size);
    case 2 % mosaic initial distribution: virus
        grid_mat = [1 1 0; 1 1 0; 0 0 0];
        grid_mat = repmat(grid_mat, round(grid_size / 3));
        qs = round((grid_size-5)/2);
        qz = round((grid_size-6)/2);
        grid_mat(qz,qs) = 1;
    case 3 % initial distribution: pentomino
        grid_mat = zeros(grid_size);
        q = round(grid_size/2);
        grid_mat(q,q-1:q) = 1;
        grid_mat(q-1,q:q+1) = 1;
        grid_mat(q+1,q) = 1;
    case 4 % glider gun
        grid_mat = zeros(grid_size);
        block1 = ones(2);
        block2 = ones(3,2);
        cannon1 = [0,0,0,1; 0,0,1,1; 0,1,1,0; 1,1,1,0; 0,1,1,0; 0,0,1,1; 0,0,0,1];
        cannon2 = [0,0,0,0,1; 0,0,1,0,1; 0,1,0,1,0; 1,0,0,1,0; 0,1,0,1,0; 0,0,1,0,1; 0,0,0,0,1];
        glider = [0,1,0; 0,0,1; 1,1,1];
        grid_mat(7:8,2:3) = block1;
        grid_mat(5:6,36:37) = block1;
        grid_mat(7:9,18:19) = block2;
        grid_mat(5:11,11:14) = cannon1;
        grid_mat(3:9,21:25) = cannon2;
        grid_mat(12:14,25:27) = glider;
end

p = size(grid_mat,1);
i = [p 1:p-1];
j = [2:p 1];

x = sparse(grid_mat);
live = zeros(time_steps);

for t=1:time_steps
    z = nnz(x);
    live(t) = z;
    
    y = x + x(:,i) + x(:,j) + x(i,:) + x(j,:) + x(i,i) + x(j,j) + x(i,j) + x(j,i);
    x = (x & (y == 4)) | (y == 3);
    
    image(logical(x))
    title({'Conway Game of Life'; num2str(t)});
    map = [1 1 1;0 0 0];
    colormap(map);
    pause(10^-1);
end

density = live / (grid_size^2);
figure
loglog(density)
xlabel('Time Steps')
ylabel('Number of Living Cells / Number of Dead Cells')
title('Ratio of Living to Dead Cells')
grid on

end


