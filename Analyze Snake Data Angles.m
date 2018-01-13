clear;

% f = fopen('SnakeRobotTestOmni2.txt', 'r');
f = fopen('SnakeRobotTestPhantomBronchoUpperPathOmni2.txt', 'r');

% % SnakeRobotTestOmni2.txt
% % tip section start
% x0 = 62.247666804;
% y0 = 9.07466665904;
% z0 = -68.0326660156;
% 
% % tip point
% x0start = 62.7576666514;
% y0start = 9.06766668955;
% z0start = -58.0306663513;
% 
% % Starting rotation
% ozeroRotation = [0.170738305648 0.953356126944 0.248903998236; 0.97954557538 -0.191525345544 0.0616533161451; 0.106450044612 0.233286563555 -0.966563516855];

% SnakeRobotTestPhantomBronchoUpperPathOmni2.txt
% tip section start
x0 = 44.6949994405;
y0 = 3.50233333111;
z0 = -68.6476671855;

% tip point
x0start = 42.5656669617;
y0start = 5.15233330727;
z0start = -58.3190006256;

% Starting rotation
ozeroRotation = [0.758324853579 -0.603657660882 0.246044113735; -0.613446325064 -0.788513358434 -0.0438968043774; 0.220507711669 -0.117646974574 -0.968263787031];

% % SnakeRobotTestPhantomBronchoLowerPathOmni2.txt
% % tip section start
% x0 = -7.5960000515;
% y0 = 16.2443333308;
% z0 = -76.4306673686;
% 
% % tip point
% x0start = -2.45133337975;
% y0start = 13.2919999758;
% z0start = -70.4213343302;

% % Starting rotation
% ozeroRotation = [-0.830665459236 0.311133266489 -0.461716670791; 0.147198076546 0.922498273849 0.356819146872; 0.536951629321 0.228433859348 -0.812094219526];

% % SnakeRobotTestPhantomBronchoHighestPathOmni.txt
% % tip section start
% x0 = 62.247666804;
% y0 = 9.07466665904;
% z0 = -68.0326660156;
% 
% % tip point
% x0start = 62.7576666514;
% y0start = 9.06766668955;
% z0start = -58.0306663513;

% % Starting rotation
% ozeroRotation = [0.351590564847 0.0284305817758 0.935721979539; 0.927066872517 0.128353581329 -0.352238501112; -0.130117607613 0.991320641836 0.0187707973644];

hold on
colorMap = [[1 0 1]; [0 1 1]; [0 1 0]; [1 0 0]; [0 0 1]];

ox0 = x0;
oy0 = y0;
oz0 = z0;
%scatter3(0, 0, 0 , 40, 'black', 'filled')

newOrigin = [0,0,0];

zeroVector = [100;100;100];
ozeroVector = [x0start - x0, y0start - y0, z0start - z0];
ozeroVector = [ozeroVector(1)/norm(ozeroVector); ozeroVector(2)/norm(ozeroVector); ozeroVector(3)/norm(ozeroVector)];

% loop through all data points
while (true)
    % read next line
    line = fgetl(f);
    if line == -1
        break;
    end
    empty = fgetl(f);
 
    rowdata1 = fgetl(f);
    rowdata2 = fgetl(f);
    rowdata3 = fgetl(f);

    row1 = strsplit(rowdata1);
    row2 = strsplit(rowdata2);
    row3 = strsplit(rowdata3);
    
    % Line info
    lineData = strsplit(line, ';');
    data = lineData(1);
    data2 = data{1:1};
    TBP = str2double(data2(5:end));
    data = lineData(2);
    data2 = data{1:1};
    MBP = str2double(data2(6:end));
    data = strtok(lineData(3),'(');
    data2 = data{1:1};
    TA = str2double(data2(6:end));
    data = strtok(lineData(4),'(');
    data2 = data{1:1};
    MA = str2double(data2(9:end));
    % End line info
    
    rotMat = [str2double(row1(1)) str2double(row1(2)) str2double(row1(3)); str2double(row2(1)) str2double(row2(2)) str2double(row2(3)); str2double(row3(1)) str2double(row3(2)) str2double(row3(3))];
    

    if TBP == 0.0 && TA == 0.0
        rotChange = rotMat * ozeroRotation';
        zeroVector = rotChange * ozeroVector;
        
        x0start = str2double(row1(4));
        y0start = str2double(row2(4));
        z0start = str2double(row3(4));

        x0 = x0start - zeroVector(1)*10;
        y0 = y0start - zeroVector(2)*10;
        z0 = z0start - zeroVector(3)*10;
        
        zeroRotation = rotMat;
        
        newOrigin = [x0 - ox0, y0 - oy0, z0 - oz0];
        scatter3(newOrigin(1), newOrigin(2), newOrigin(3) , 40, 'black', 'filled')
        
        planeZeroVector = zeroVector;
        if MA == 0.0 && 1
            plotPlane(planeZeroVector, newOrigin);
        end
        if MA == 31.0 && MBP == 0.0 && 1
            plotPlane(planeZeroVector, newOrigin);
        end                
        if MA == 31.0 && MBP == 120.0 && 1
            plotPlane(planeZeroVector, newOrigin);
        end
        if MA == 31.0 && MBP == 240.0 && 1
            plotPlane(planeZeroVector, newOrigin);
        end
    end    
    
    rotChange = rotMat * zeroRotation';
    diVec = rotChange * zeroVector;
    
    x = str2double(row1(4)) - ox0; 
    y = str2double(row2(4)) - oy0;
    z = str2double(row3(4)) - oz0;
    
    if MA == 0.0 && 1
        scatter3(x,y,z, 30, 'red', 'filled')
        
         v = [x y z];
         v2 = [x+diVec(1), y+diVec(2), z+diVec(3)];
         vc = [v2;v];
         plot3(vc(:,1),vc(:,2),vc(:,3),'red');
         v2 = [x+zeroVector(1), y+zeroVector(2), z+zeroVector(3)];
         vc = [v2;v];
         plot3(vc(:,1),vc(:,2),vc(:,3),'black');

        linep1 = [newOrigin(1), newOrigin(2), newOrigin(3)];
        linep2 = [x y z];
        lineps = [linep2; linep1];
        plot3(lineps(:,1),lineps(:,2),lineps(:,3),'red');
        
        % Calculate theta (real angle error)
        zv = [0; 0; 1];
        zvr = rotChange * zv;
        theta = radtodeg(atan2(norm(cross(zvr,zv)),dot(zvr,zv)));
        text(v2(1), v2(2), v2(3), num2str(theta));
        
        % Calculate phi (projected angle error)
        if TA == 0.0
            if TBP == 0.0
                TA0zerovector = zvr;
                %TA0projectedzerovector = TA0zerovector - dot(TA0zerovector, planeZeroVector) * planeZeroVector;
                TA0projectedzerovector = TA0zerovector - dot(TA0zerovector, zv) * zv;
                TA0projectedzerovector = Threshold(TA0projectedzerovector);
            end
            TA0vector = zvr;
            TA0projectedvector = TA0vector - dot(TA0vector, zv) * zv;
            TA0projectedvector = Threshold(TA0projectedvector);
            phi = radtodeg(atan2(norm(cross(TA0projectedvector,TA0projectedzerovector)),dot(TA0projectedvector,TA0projectedzerovector)));
        end        
            
        if TA == 34.0
            if TBP == 0.0
                TA34zerovector = zvr;
                TA34projectedzerovector = TA34zerovector - dot(TA34zerovector, zv) * zv;
            end
            TA34vector = zvr;
            TA34projectedvector = TA34vector - dot(TA34vector, zv) * zv;
            phi = radtodeg(atan2(norm(cross(TA34projectedvector,TA34projectedzerovector)),dot(TA34projectedvector,TA34projectedzerovector)));
        end
        
        if TA == 68.0
            if TBP == 0.0
                TA68zerovector = zvr;
                TA68projectedzerovector = TA68zerovector - dot(TA68zerovector, zv) * zv;
            end
            TA68vector = zvr;
            TA68projectedvector = TA68vector - dot(TA68vector, zv) * zv;
            phi = radtodeg(atan2(norm(cross(TA68projectedvector,TA68projectedzerovector)),dot(TA68projectedvector,TA68projectedzerovector)));
        end
        
        display(['TBP: ', num2str(TBP), '; MBP: ', num2str(MBP), '; Tip: ', num2str(TA), '; Middle: ', num2str(MA), '; Theta: ', num2str(theta), '; Phi: ', num2str(phi)]);
    end
    
    if MA == 31.0 && MBP == 0.0 && 1
        scatter3(x,y,z, 30, 'blue', 'filled')
        v = [x y z];
        v2 = [x+diVec(1), y+diVec(2), z+diVec(3)];
        vc = [v2;v];
        plot3(vc(:,1),vc(:,2),vc(:,3),'blue');
        v2 = [x+zeroVector(1), y+zeroVector(2), z+zeroVector(3)];
        vc = [v2;v];
        plot3(vc(:,1),vc(:,2),vc(:,3),'black');        
        
        linep1 = [newOrigin(1), newOrigin(2), newOrigin(3)];
        linep2 = [x y z];
        lineps = [linep2; linep1];
        plot3(lineps(:,1),lineps(:,2),lineps(:,3),'blue');        
        
        % Calculate theta (real angle error)
        zv = [0; 0; 1];
        zvr = rotChange * zv;
        theta = radtodeg(atan2(norm(cross(zvr,zv)),dot(zvr,zv)));
        text(v2(1), v2(2), v2(3), num2str(theta));        
        
        % Calculate phi (projected angle error)
        if TA == 0.0
            if TBP == 0.0
                TA0zerovector = zvr;
                %TA0projectedzerovector = TA0zerovector - dot(TA0zerovector, planeZeroVector) * planeZeroVector;
                TA0projectedzerovector = TA0zerovector - dot(TA0zerovector, zv) * zv;
                TA0projectedzerovector = Threshold(TA0projectedzerovector);
            end
            TA0vector = zvr;
            TA0projectedvector = TA0vector - dot(TA0vector, zv) * zv;
            TA0projectedvector = Threshold(TA0projectedvector);
            phi = radtodeg(atan2(norm(cross(TA0projectedvector,TA0projectedzerovector)),dot(TA0projectedvector,TA0projectedzerovector)));
        end        
            
        if TA == 34.0
            if TBP == 0.0
                TA34zerovector = zvr;
                TA34projectedzerovector = TA34zerovector - dot(TA34zerovector, zv) * zv;
            end
            TA34vector = zvr;
            TA34projectedvector = TA34vector - dot(TA34vector, zv) * zv;
            phi = radtodeg(atan2(norm(cross(TA34projectedvector,TA34projectedzerovector)),dot(TA34projectedvector,TA34projectedzerovector)));
        end
        
        if TA == 68.0
            if TBP == 0.0
                TA68zerovector = zvr;
                TA68projectedzerovector = TA68zerovector - dot(TA68zerovector, zv) * zv;
            end
            TA68vector = zvr;
            TA68projectedvector = TA68vector - dot(TA68vector, zv) * zv;
            phi = radtodeg(atan2(norm(cross(TA68projectedvector,TA68projectedzerovector)),dot(TA68projectedvector,TA68projectedzerovector)));
        end
        
        display(['TBP: ', num2str(TBP), '; MBP: ', num2str(MBP), '; Tip: ', num2str(TA), '; Middle: ', num2str(MA), '; Theta: ', num2str(theta), '; Phi: ', num2str(phi)]);        
    end
    
    if MA == 31.0 && MBP == 120.0 && 0
        scatter3(x,y,z, 30, 'green', 'filled')
        v = [x y z];
        v2 = [x+diVec(1), y+diVec(2), z+diVec(3)];
        vc = [v2;v];
        plot3(vc(:,1),vc(:,2),vc(:,3),'green');
        v2 = [x+zeroVector(1), y+zeroVector(2), z+zeroVector(3)];
        vc = [v2;v];
        plot3(vc(:,1),vc(:,2),vc(:,3),'black');           
        
        linep1 = [newOrigin(1), newOrigin(2), newOrigin(3)];
        linep2 = [x y z];
        lineps = [linep2; linep1];
        plot3(lineps(:,1),lineps(:,2),lineps(:,3),'green');        
        
        % Calculate theta (real angle error)
        zv = [0; 0; 1];
        zvr = rotChange * zv;
        theta = radtodeg(atan2(norm(cross(zvr,zv)),dot(zvr,zv)));
        text(v2(1), v2(2), v2(3), num2str(theta));
        
        % Calculate phi (projected angle error)
        if TA == 0.0
            if TBP == 0.0
                TA0zerovector = zvr;
                %TA0projectedzerovector = TA0zerovector - dot(TA0zerovector, planeZeroVector) * planeZeroVector;
                TA0projectedzerovector = TA0zerovector - dot(TA0zerovector, zv) * zv;
                TA0projectedzerovector = Threshold(TA0projectedzerovector);
            end
            TA0vector = zvr;
            TA0projectedvector = TA0vector - dot(TA0vector, zv) * zv;
            TA0projectedvector = Threshold(TA0projectedvector);
            phi = radtodeg(atan2(norm(cross(TA0projectedvector,TA0projectedzerovector)),dot(TA0projectedvector,TA0projectedzerovector)));
        end        
            
        if TA == 34.0
            if TBP == 0.0
                TA34zerovector = zvr;
                TA34projectedzerovector = TA34zerovector - dot(TA34zerovector, zv) * zv;
            end
            TA34vector = zvr;
            TA34projectedvector = TA34vector - dot(TA34vector, zv) * zv;
            phi = radtodeg(atan2(norm(cross(TA34projectedvector,TA34projectedzerovector)),dot(TA34projectedvector,TA34projectedzerovector)));
        end
        
        if TA == 68.0
            if TBP == 0.0
                TA68zerovector = zvr;
                TA68projectedzerovector = TA68zerovector - dot(TA68zerovector, zv) * zv;
            end
            TA68vector = zvr;
            TA68projectedvector = TA68vector - dot(TA68vector, zv) * zv;
            phi = radtodeg(atan2(norm(cross(TA68projectedvector,TA68projectedzerovector)),dot(TA68projectedvector,TA68projectedzerovector)));
        end
        
        display(['TBP: ', num2str(TBP), '; MBP: ', num2str(MBP), '; Tip: ', num2str(TA), '; Middle: ', num2str(MA), '; Theta: ', num2str(theta), '; Phi: ', num2str(phi)]);        
    end
    
    if MA == 31.0 && MBP == 240.0 && 0
        scatter3(x,y,z, 30, 'magenta', 'filled')
        v = [x y z];
        v2 = [x+diVec(1), y+diVec(2), z+diVec(3)];
        vc = [v2;v];
        plot3(vc(:,1),vc(:,2),vc(:,3),'magenta');
        v2 = [x+zeroVector(1), y+zeroVector(2), z+zeroVector(3)];
        vc = [v2;v];
        plot3(vc(:,1),vc(:,2),vc(:,3),'black');           
        
        linep1 = [newOrigin(1), newOrigin(2), newOrigin(3)];
        linep2 = [x y z];
        lineps = [linep2; linep1];
        plot3(lineps(:,1),lineps(:,2),lineps(:,3),'magenta');        
        
        % Calculate theta (real angle error)
        zv = [0; 0; 1];
        zvr = rotChange * zv;
        theta = radtodeg(atan2(norm(cross(zvr,zv)),dot(zvr,zv)));
        text(v2(1), v2(2), v2(3), num2str(theta));
        
        % Calculate phi (projected angle error)
        if TA == 0.0
            if TBP == 0.0
                TA0zerovector = zvr;
                %TA0projectedzerovector = TA0zerovector - dot(TA0zerovector, planeZeroVector) * planeZeroVector;
                TA0projectedzerovector = TA0zerovector - dot(TA0zerovector, zv) * zv;
                TA0projectedzerovector = Threshold(TA0projectedzerovector);
            end
            TA0vector = zvr;
            TA0projectedvector = TA0vector - dot(TA0vector, zv) * zv;
            TA0projectedvector = Threshold(TA0projectedvector);
            phi = radtodeg(atan2(norm(cross(TA0projectedvector,TA0projectedzerovector)),dot(TA0projectedvector,TA0projectedzerovector)));
        end        
            
        if TA == 34.0
            if TBP == 0.0
                TA34zerovector = zvr;
                TA34projectedzerovector = TA34zerovector - dot(TA34zerovector, zv) * zv;
            end
            TA34vector = zvr;
            TA34projectedvector = TA34vector - dot(TA34vector, zv) * zv;
            phi = radtodeg(atan2(norm(cross(TA34projectedvector,TA34projectedzerovector)),dot(TA34projectedvector,TA34projectedzerovector)));
        end
        
        if TA == 68.0
            if TBP == 0.0
                TA68zerovector = zvr;
                TA68projectedzerovector = TA68zerovector - dot(TA68zerovector, zv) * zv;
            end
            TA68vector = zvr;
            TA68projectedvector = TA68vector - dot(TA68vector, zv) * zv;
            phi = radtodeg(atan2(norm(cross(TA68projectedvector,TA68projectedzerovector)),dot(TA68projectedvector,TA68projectedzerovector)));
        end
        
        display(['TBP: ', num2str(TBP), '; MBP: ', num2str(MBP), '; Tip: ', num2str(TA), '; Middle: ', num2str(MA), '; Theta: ', num2str(theta), '; Phi: ', num2str(phi)]);        
    end    
    
    empty = fgetl(f);
    empty = fgetl(f);
end


%{
for MP = [0 60 120 180 240 300]
    for MA = [0 40 80 -40 -80]
        MP_index = MP/60 + 1;

        if MA < 0
            MA_index = -MA/40 + 3;
        elseif MA >= 0
            MA_index = MA/40 + 1;
        end        

        avgPoint = [mean(points(MP_index, MA_index, :, 1)),mean(points(MP_index, MA_index, :, 2)), mean(points(MP_index, MA_index, :, 3))];
        sphereRadius = crr([std(points(MP_index, MA_index, :, 1)), std(points(MP_index, MA_index, :, 2)), std(points(MP_index, MA_index, :, 3))]);
        [sphX, sphY, sphZ] = sphere(10);
        m = mesh(sphX*sphereRadius + avgPoint(1), sphY*sphereRadius + avgPoint(2), sphZ*sphereRadius + avgPoint(3));
        set(m, 'facecolor', 'none', 'edgecolor', colorMap(MA_index, :))
        
        
        text(avgPoint(1), avgPoint(2), avgPoint(3), ['MA: ' num2str(MA)]); 
        
    end
    t = text(avgPoint(1), avgPoint(2), avgPoint(3) - 4, ['MP: ' num2str(MP)]); 
    t.FontSize = 2;
end



% zeroed robot for reference (straight)
xs = zeros(1,30);
ys = zeros(1,30);
zs = linspace(0, 30, 30);
plot3(xs, ys, zs, 'b')

% line in 0 degrees bending plane for reference
xr = linspace(0, 2, 30);
yr = zeros(1, 30);
plot(xr, yr, 'r')

xlim([-30 30])
ylim([-30 30])
zlim([0 35])
daspect([1 1 1])

hold off
rotate3d on
fclose('all');

%}

fclose('all');
