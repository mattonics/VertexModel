close all;

origHex = [[1, 0];[1/2,sqrt(3)/2];[-1/2, sqrt(3)/2];[-1 0];[-1/2 -sqrt(3)/2];[1/2 -sqrt(3)/2]];

Vertices = [];
VertexIDs = {};

rowsize = 30;
columnsize = 30;

%Make row of hexagons
for i = 1:rowsize
    rightTranslation = (i - 1)*3;
    newVertices = origHex + [rightTranslation 0];
    Vertices = [Vertices;newVertices];
    VertexIDs = [VertexIDs, {(1:6)+(i-1)*6}];
end

%Fill in gaps between gaps using a shifted row
for i = 1:length(Vertices)
    Vertices = [Vertices;Vertices(i,:)+[1.5, -sqrt(3)/2]];
    if (mod(i,6) == 0)
       VertexIDs = [VertexIDs {VertexIDs{end}+6}];
    end
end

%Replicate row pair for many rows
for i = 1:columnsize
   Vertices = [Vertices;[Vertices(1:360,:) + [0, -sqrt(3)]*i]];
   for j = 1:rowsize*2
      VertexIDs = [VertexIDs {VertexIDs{end}+6}];
   end
end

%Join vertices sharing a vertex into a single vertex
for i = 1:length(VertexIDs)
    currIDs = VertexIDs{i};
    for j = 1:length(currIDs)
        currVert = Vertices(currIDs(j),:);
        pos = find(abs(sum(Vertices - currVert,2)) < 10e-8);
        if (~isempty(pos))
            VertexIDs{i}(j) = pos(1);
        end
    end
end

%Select cells within a bounding circle
center = [30 -25];
delinds = zeros(length(VertexIDs),1);
r = 15;

for i = 1:length(VertexIDs)
    currIDs  = VertexIDs{i};
    for j = 1:length(currIDs)
        point = Vertices(currIDs(j),:);
        if ((center(1) - point(1))^2 + (center(2) - point(2))^2 > r^2)
            delinds(i) = 1;
            break;
        end
    end
end

v = VideoWriter('D:\Vertex Model\vm.avi');
open(v);

VertexIDs(logical(delinds)) = [];

lambda = 7;
a0 = 3*sqrt(3)/2*1.25;
c0 = 2*sqrt(pi*a0);
beta  = 4;
rng(1)

gamma = 1;
nu = 1;
deltat = 0.0055;
v0 = 2;
gammarandom = 0.25;
betaP = 2;
alphaP = 0.3;

mu_gla = 1;

hi = 1;

T1min = 0.1;

AllVels = {};

h = VisualizeVertexModel(VertexIDs,Vertices);

neighbors = {};
for i = 1:length(VertexIDs)
    currIDs = VertexIDs{i};
    currNeighbors = [];
    for j = 1:length(currIDs)
        if (j == length(currIDs))
            nextInd = 1;
        else
            nextInd = j + 1;
        end
        
        for k = 1:length(VertexIDs)
            if k == i
                continue;
            end
            
            tempIDs = VertexIDs{k};

            
            for l = 1:length(tempIDs)
                if(l == length(tempIDs))
                    tempNextInd = 1;
                else
                    tempNextInd = l + 1;
                end
                
                if (all([tempIDs(tempNextInd) tempIDs(l)] == [currIDs(j) currIDs(nextInd)]) | all([tempIDs(tempNextInd) tempIDs(l)] == [currIDs(nextInd) currIDs(j)]))
                    currNeighbors = [currNeighbors k];
                end                
            end           
        end
    end
    
    neighbors = [neighbors {currNeighbors}];
end

frames = [];
thetajs = rand(length(VertexIDs),1)*2*pi;

veljs = zeros(length(VertexIDs),2);

for step = 1:1000
    disp(string(step));
    
    for i = 1:length(thetajs)
        curr_neighbors = neighbors{i};
        gla = 0;
        for j = 1:length(curr_neighbors)
            diff_vec = veljs(curr_neighbors(j),:) - [cos(thetajs(i)) sin(thetajs(i))];
            angle = atan(diff_vec(2)/diff_vec(1));
            gla = gla + mu_gla*angle/length(curr_neighbors);
        end
        
        thetajs(i) = thetajs(i) + gla;
    end
    
    oldcentroids = zeros(length(VertexIDs),2);
    for i = 1:length(VertexIDs)
        oldcentroids(i,:) = mean(Vertices(VertexIDs{i},:));
    end
    
    oldPositions = Vertices;
    forces = zeros(size(Vertices,1),size(Vertices,2));
    areas = zeros(length(VertexIDs),1);
    for i = 1:length(VertexIDs)
        areas(i) = area(polyshape(Vertices(VertexIDs{i},1),Vertices(VertexIDs{i},2)));
    end
    
    perims = zeros(length(VertexIDs),1);
    for i = 1:length(VertexIDs)
        perims(i) = perimeter(polyshape(Vertices(VertexIDs{i},1),Vertices(VertexIDs{i},2)));
    end
    
    for i = 1:length(VertexIDs)
        
        currIDs = VertexIDs{i};
        for j = 1:length(currIDs)
            k = currIDs(j);
            currVert = Vertices(k,:);
            if (j == length(currIDs))
                nextInd = 1;
            else
                nextInd = j+1;
            end
            if (j == 1)
                prevInd = length(currIDs);
            else
                prevInd = j-1;
            end
            
            nextVert = Vertices(currIDs(nextInd),:);
            prevVert = Vertices(currIDs(prevInd),:);
            
            deltaiAl = 1/2*[nextVert(2)-prevVert(2),prevVert(1)-nextVert(1)];
           
            AreaTerm = 2*lambda*deltaiAl*(areas(i)-a0);   
            
            deltaidl1 = 1/norm(prevVert - currVert) * (currVert - prevVert);
            deltaidl2 = 1/norm(nextVert - currVert) * (currVert - nextVert);
            
            PerimeterTerm = 2*beta*(perims(i) - c0)*(deltaidl1+deltaidl2);

            ContractileTerm = gamma*(deltaiAl)*areas(i);
            
            RandomTerm = gammarandom * randn(1,2)/length(neighbors{i});
            
            PersistenceTerm = [0 0];
            
            for m = 1:length(AllVels)
                PersistenceTerm = PersistenceTerm + exp(-betaP * (step - m))*AllVels{m}(k);
            end
            
            if (~all(PersistenceTerm == 0))
                PersistenceTerm = PersistenceTerm/norm(PersistenceTerm)*alphaP;
            end
            
            currneighbors = neighbors{i};
            PolarityTerm = [0 0];
            for m = 1:length(currneighbors)
                PolarityTerm = PolarityTerm + v0 * [cos(thetajs(currneighbors(m))) sin(thetajs(currneighbors(m)))];
            end
            
            PolarityTerm = PolarityTerm/length(currneighbors);
            
            forces(k,:) = forces(k,:) - AreaTerm - PerimeterTerm - RandomTerm + PersistenceTerm + PolarityTerm;
        end
    end
        
    Vertices = Vertices + forces*deltat/nu;  
    
    AllVels = [AllVels {forces*deltat/nu}];
    
    %Check boundary conditions
    
    if (step > 1)
        
        testVertices = unique([VertexIDs{:}]);
        
        for i = 1:length(testVertices)
            currPoint = Vertices(testVertices(i),:);
            oldPoint = oldPositions(testVertices(i),:);
            if ((currPoint(1) - center(1))^2 + (currPoint(2) - center(2))^2 > r^2)
                syms xvar;

                p1x = currPoint(1);
                p1y = currPoint(2);

                p2x = oldPoint(1);
                p2y = oldPoint(2);

                slope = (p2y-p1y) / (p2x-p1x);

                b = p2y - p2x * slope;

                [xvalues, yvalues] = linecirc(slope, b, center(1), center(2), r);
                
                if (norm(currPoint - [xvalues(1) yvalues(1)]) < norm(currPoint - [xvalues(2) yvalues(2)]))
                    Vertices(testVertices(i),:) = [xvalues(1) yvalues(1)];
                else
                    Vertices(testVertices(i),:) = [xvalues(2) yvalues(2)];
                end
            end
        end
    end
    %Implement T1 Vertex Transitions
    prevT1s = [0 0];

    for i = 1:length(VertexIDs)
        currIDs = VertexIDs{i};
        for j = 1:length(currIDs)
            if (j == length(currIDs))
                next = 1;
            else
                next = j + 1;
            end
            dist = norm(Vertices(currIDs(j),:) - Vertices(currIDs(next),:));
            
            if (dist < T1min && sum([VertexIDs{:}] == currIDs(j)) >= 2 && sum([VertexIDs{:}] == currIDs(next)) >= 2)
                if (any(sum([currIDs(j) currIDs(next)] == prevT1s, 2) == 2) || any(sum([currIDs(next) currIDs(j)] == prevT1s, 2) == 2))
                    continue;
                end
                                
                cells2points = [];
                cells1point = [];
                for m = 1:length(VertexIDs)
                    tempIDs = VertexIDs{m};
                    if (sum(ismember([currIDs(j) currIDs(next)], tempIDs)) == 2)
                        cells2points = [cells2points m];
                    elseif (sum(ismember([currIDs(j) currIDs(next)], tempIDs)) == 1)
                        cells1point = [cells1point m];
                    end
                end
                
                p1x = Vertices(currIDs(j),1);
                p1y = Vertices(currIDs(j),2);
                
                p2x = Vertices(currIDs(next),1);
                p2y = Vertices(currIDs(next),2);
                
                mix = mean([p1x p2x]);
                miy = mean([p1y p2y]);
                
                slope = (p2y-p1y) / (p2x-p1x);
                slope = -1/slope;
                
                ydist = abs(p1y - p2y);
                xvalues = [mix - ydist/2 * hi, mix + ydist/2 * hi];
                
                yvalues = slope * (xvalues - mix) + miy;
                
                Vertices = [Vertices;[xvalues(1) yvalues(1)]];
                Vertices = [Vertices;[xvalues(2) yvalues(2)]];
                
                for m = 1:length(cells2points)
                    
                    tempIDs = VertexIDs{cells2points(m)};
                    
                    jindex = find(tempIDs == currIDs(j));
                    nextindex = find(tempIDs == currIDs(next));
                    
                    centroid = mean(Vertices(tempIDs,:));
                    
                    tempIDs1 = tempIDs;
                    tempIDs1(jindex) = length(Vertices) - 1;
                    tempIDs1(nextindex) = [];

                    tempIDs2 = tempIDs;
                    tempIDs2(jindex) = length(Vertices);
                    tempIDs2(nextindex) = [];
                    
                    oldpolygon = polyshape(Vertices(tempIDs,1),Vertices(tempIDs,2));
                    
                    if (inpolygon(Vertices(length(Vertices)-1,1),Vertices(length(Vertices)-1,2),Vertices(tempIDs,1),Vertices(tempIDs,2)))
                        tempIDs = tempIDs1;
                    else
                        tempIDs = tempIDs2;
                    end
%                     if (perimeter(polyshape(Vertices(tempIDs1,1),Vertices(tempIDs1,2))) < perimeter(polyshape(Vertices(tempIDs2,1),Vertices(tempIDs2,2))))
%                         tempIDs = tempIDs1;
%                     else
%                         tempIDs = tempIDs2;
%                     end
            
                    VertexIDs{cells2points(m)} = tempIDs;
                    
                    currNeighbors = neighbors{cells2points(m)};
                    if (length(cells2points) == 2)
                        if (m == 1)
                            currNeighbors(currNeighbors == cells2points(2)) = [];
                        else
                            currNeighbors(currNeighbors == cells2points(1)) = [];
                        end
                        neighbors{cells2points(m)} = currNeighbors;
                    end
                end
                
                for m = 1:length(cells1point)
                    
                    tempIDs = VertexIDs{cells1point(m)};
                    
                    pointindex = find((tempIDs == currIDs(j)) + (tempIDs == currIDs(next)));
                    
                    if (pointindex == 1)
                        prevInd = length(tempIDs);
                    else
                        prevInd = pointindex - 1;
                    end
                    
                    if (pointindex == length(tempIDs))
                        nextInd = 1;
                    else
                        nextInd = pointindex + 1;
                    end
                    
                    for n = 1:length(cells2points)
                        if (ismember(tempIDs(prevInd), VertexIDs{cells2points(n)}))
                            if (ismember(length(Vertices) - 1, VertexIDs{cells2points(n)}))
                                tempIDs = [tempIDs(1:pointindex - 1) length(Vertices) - 1 length(Vertices) tempIDs(pointindex + 1:end)];
                                break
                            elseif (ismember(length(Vertices), VertexIDs{cells2points(n)}))
                                tempIDs = [tempIDs(1:pointindex - 1) length(Vertices) length(Vertices) - 1 tempIDs(pointindex + 1:end)];
                                break
                            end
                        end
                    end
                    VertexIDs{cells1point(m)} = tempIDs;
                    
                    if (length(cells1point) == 2)
                        if (m == 1)
                            otherCell = cells1point(2);
                        else
                            otherCell = cells1point(1);
                        end
                        neighbors{cells1point(m)} = [neighbors{cells1point(m)}  otherCell];
                    end
                end
                
                prevT1s = [prevT1s;[length(Vertices)-1 length(Vertices)]];
                prevT1s = [prevT1s;[length(Vertices) length(Vertices)-1]];
                
                break;                        
            end
        end
    end
    
    for i = 1:length(VertexIDs)
        centroid = mean(Vertices(VertexIDs{i},:));
        veljs(i,:) = centroid - oldcentroids(i,:);
    end
    
    h = VisualizeVertexModel(VertexIDs,Vertices);
    frame = getframe(h);
    writeVideo(v, frame);
end

close(v);



