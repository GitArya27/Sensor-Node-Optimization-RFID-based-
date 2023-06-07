%this is BBO-DE algo.
% we make generation containing a solution of n population of positions of
% sensers

%initializing the variables that we are going to use (you can skip it )
sensingRange = 50;
communicationRange = 50;
maxIteration = 100;
npopulationSize = 100;
bsLocation = [250,250,250];
nTargetPoints = 50; %15 to 150
bbopopulationsize = 200;
fitness = zeros(bbopopulationsize,1);
targetpoints = zeros(nTargetPoints,3);% Targets in 3-D space with cordinates


%randomly distribute the target points in 3-d plane
for i = 1:nTargetPoints
    for k = 1:3
        targetpoints(i,k) = randi(500);
    end
end

%setAppropriatePositions contain set of all possible positions where
%sensers can be plased


generation = zeros(npopulationSize,3,bbopopulationsize); % Generation of diffrent Population

for i = 1:bbopopulationsize
    for j = 1: npopulationSize
        for k = 1:3
            generation(j,k,i)= randi(500);
        end
    end
end





% plot the initial generation of population on 3-d space
plot_nodes(generation(:,:,1),targetpoints,generation(:,:,1),sensingRange,"initial position");

%
%
%                                     main code(main code starts from here)

for i = 1:maxIteration
    
    % Update the population using the BBO-DE algorithm
    newGeneration = bbo_de(generation,nTargetPoints,sensingRange,communicationRange,bsLocation,targetpoints,bbopopulationsize,npopulationSize);
    generation = newGeneration;
    if i == 10 || i == 20 || i == 30 || i == 40 || i == 50 || i == 60 || i ==70 || i == 80 || i ==90
        str1 = "Position after " ;
        str2 = int2str(i);
        str3 = " Iteration";
        str12 = strcat(str1 , str2);
        str = strcat(str12,str3);
        plot_nodes(generation(:,:,1),targetpoints,generation(:,:,1),sensingRange,str );
    end
end


plot_nodes(generation(:,:,1),targetpoints,generation(:,:,1),sensingRange,"final position");



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculation for the final output (can skip it and come on it later)
population = generation(:,:,1);
positions = generation(:,:,1);
number_of_targets_covered =0;
for j = 1: nTargetPoints
    for i = 1:npopulationSize
        x1 = targetpoints(j,1);
        y1 = targetpoints(j,2);
        z1 = targetpoints(j,3);
        x2 = positions(i,1);
        y2 = positions(i,2);
        z2 = positions(i,3);
        distance = sqrt( ( (x1-x2)*(x1-x2)  ) + (  (y2-y1)*(y2-y1)  )+ (  (z1-z2)*(z1-z2)  )  );
        if distance <= sensingRange
            number_of_targets_covered = number_of_targets_covered +1;
            break;
        end
    end
end
number_of_connection =0;
for i = 1:size(positions,1)
    for j = 1:size(positions,1)
        x1 = positions(j,1);
        y1 = positions(j,2);
        z1 = positions(j,3);
        x2 = positions(i,1);
        y2 = positions(i,2);
        z2 = positions(i,3);
        distance = sqrt( ( (x1-x2)*(x1-x2)  ) + (  (y2-y1)*(y2-y1)  )+ (  (z1-z2)*(z1-z2)  )  );
        if distance <= 2* communicationRange
             number_of_connection = number_of_connection +1;
        end
    end
end
overlap = 0;
for i = 1:size(positions, 1)
    for j = i+1:size(positions, 1)
        x1 = positions(j,1);
        y1 = positions(j,2);
        z1 = positions(j,3);
        x2 = positions(i,1);
        y2 = positions(i,2);
        z2 = positions(i,3);
        distance = sqrt( ( (x1-x2)*(x1-x2)  ) + (  (y2-y1)*(y2-y1)  )+ (  (z1-z2)*(z1-z2)  )  );
        if distance <= 2*sensingRange
            x = distance;
            h = sensingRange - x/2;
            volofcurve = (h*h*(3*sensingRange - h)*3.14)/3;
            overlap = overlap + 2*(volofcurve);
        end
    end
end
popvol =( 3.14 *4* sensingRange * sensingRange * sensingRange * size(positions,1))/3;

number_of_node_covering_target =0;
    
    for i = 1: npopulationSize
        for j = 1:nTargetPoints
                x1 = targetpoints(j,1);
                y1 = targetpoints(j,2);
                z1 = targetpoints(j,3);
                x2 = positions(i,1);
                y2 = positions(i,2);
                z2 = positions(i,3);
            distance = sqrt( ( (x1-x2)*(x1-x2)  ) + (  (y2-y1)*(y2-y1)  )+ (  (z1-z2)*(z1-z2)  )  );
            if distance <= sensingRange
                number_of_node_covering_target = number_of_node_covering_target  +1;
            end
        end
    end
final_number_of_nodes =number_of_node_covering_target /npopulationSize
final_target_cover_ratio = number_of_targets_covered / nTargetPoints
final_connection_ratio = number_of_targets_covered/(size(positions,1)*(size(positions,1)-1))
final_overlap_ratio = overlap/popvol












%
%                                               BBO-DE ALGO
%

function newGeneration= bbo_de(generation,nTargetPoints,sensingRange,communicationRange,bsLocation,targetpoints,bbopopulationsize,npopulationSize)
    % Set the algorithm parameters
    %mutationprobability = 0.2
    crossoverprobability = 0.2;
    F = zeros(npopulationSize,3); % F belongs to 0.4 to 0.9
    for i = 1:npopulationSize
        for j = 1:3
            F(i,j)=0.4;
        end
    end
    
    % setting emigration rate
    emigrationrate = zeros(bbopopulationsize,1);
    immigrationrate = zeros(bbopopulationsize,1);
    for i = 1:bbopopulationsize
        emigrationrate(i,1) = (bbopopulationsize-i)/bbopopulationsize;
        immigrationrate(i,1)= 1- emigrationrate(i,1);
    end
    %calculating fitness of all generation
    fitness = zeros(bbopopulationsize,1);
    for i = 1:bbopopulationsize
        fitness(i) = evaluateFitness(generation(:,:,i),nTargetPoints, sensingRange, communicationRange, bsLocation,targetpoints,npopulationSize);
    end
    [~, sortIndex] = sort(fitness(1:bbopopulationsize), 'descend');
    emigrationrate = emigrationrate(sortIndex,1);
    immigrationrate = immigrationrate(sortIndex,1);
    for i = 1:bbopopulationsize
        for j = 1:npopulationSize
            jrand = rand();
            if jrand < immigrationrate(i,1)
                immigratingpos = -1;
                for k = 1:bbopopulationsize
                    if generation(j,:,i)==generation(j,:,k)
                        continue;
                    end
                    jjrand = rand();
                    if jjrand < emigrationrate(k,1)
                        immigratingpos = k;
                        break;
                    end
                end
                if immigratingpos ~=-1
                    generation(j,:,i)=generation(j,:,immigratingpos);
                end                
            end
        end
    end
    
    %Habitual Mutation
    for i = 1:bbopopulationsize
        r = randi(bbopopulationsize,3,1);
        mutantvector = generation(:,:,r(1,1)) + F.*(generation(:,:,r(2,1))- generation(:,:,r(3,1)));
        for k = 1:npopulationSize
            if mutantvector(k,1) > 500 || mutantvector(k,1) < 0 
                mutantvector(k,1) = randi(500);
            end
            if mutantvector(k,2) > 500 || mutantvector(k,2) < 0 
                mutantvector(k,2) = randi(500);
            end
            if mutantvector(k,3) > 500 || mutantvector(k,3) < 0 
                mutantvector(k,3) = randi(500);
            end
        end
        trialvector = generation(:,:,i);
        jrand = randi(npopulationSize);
        for j = 1:npopulationSize
            jjrand = rand();
            if jjrand < crossoverprobability || j==jrand 
                trialvector(j,:)=mutantvector(j,:);
            end
        end
        trialfitness = evaluateFitness(trialvector,nTargetPoints, sensingRange, communicationRange, bsLocation,targetpoints,size(trialvector,1));
        if trialfitness > fitness(i,1)
            generation(:,:,i)= trialvector;
            fitness(i)=trialfitness;
        end
    end
    [~, sortIndex] = sort(fitness(1:bbopopulationsize), 'descend');
    generation = generation(:,:,sortIndex);
    newGeneration = generation;
    
end





%fitness finction to get the accuricy of each solution

% Define the fitness function
function fitness = evaluateFitness(positions, nTargetPoints, sensingRange, communicationRange, bsLocation,targetpoints,possize)
    [f1, f2, f3, f4] = eof(positions, nTargetPoints, sensingRange, communicationRange, bsLocation,targetpoints,possize);
    w1 = 0.25;
    w2 = 0.25;
    w3 = 0.25;
    w4 = 0.25;
    fitness = w1*f1+w2*f2-w3*f3-w4*f4;
end
% Define the objective functions
function [f1, f2, f3, f4] = eof(positions, nTargetPoints, sensingRange, communicationRange, ~,targetpoints,npopulationSize)
        reward = 100;

    % Objective function 1: maximize the coverage of target points by sensors
    number_of_targets_covered =0;
    possize = npopulationSize;
    for j = 1: nTargetPoints
        for i = 1:npopulationSize
                x1 = targetpoints(j,1);
                y1 = targetpoints(j,2);
                z1 = targetpoints(j,3);
                x2 = positions(i,1);
                y2 = positions(i,2);
                z2 = positions(i,3);
            distance = sqrt( ( (x1-x2)*(x1-x2)  ) + (  (y2-y1)*(y2-y1)  )+ (  (z1-z2)*(z1-z2)  )  );
            if distance <= sensingRange
                number_of_targets_covered = number_of_targets_covered +1;
                break;
            end
        end
    end

    f1 =100*reward*number_of_targets_covered/nTargetPoints;
    % Objective function 2: maximize the connectivity of sensors
    % Objective function 3: minimize the overlap between sensors
    overlap = 0;
    number_of_connection =0;
    for i = 1:possize
        for j = i+1:possize
            x1 = positions(j,1);
            y1 = positions(j,2);
            z1 = positions(j,3);
            x2 = positions(i,1);
            y2 = positions(i,2);
            z2 = positions(i,3);
            distance = sqrt( ( (x1-x2)*(x1-x2)  ) + (  (y2-y1)*(y2-y1)  )+ (  (z1-z2)*(z1-z2)  )  );
            if distance <= 2 * sensingRange
                x = distance;
                h = sensingRange-x/2;
                volofcurve = h*h*(3*sensingRange - h)*3.14/3;
                overlap = overlap + 2*(volofcurve);
            end
            if distance <= 2*communicationRange
                 number_of_connection = number_of_connection +1;
            end
        end
    end
    number_of_connection = number_of_connection /2;
    f2 =reward*number_of_connection/(possize*(possize-1));

    popvol = 4*3.14 * sensingRange * sensingRange* sensingRange * possize/3;
    f3 =reward*overlap/popvol;

    % Objective function 4: minimize numer of sensing node covering target
    
    number_of_node_covering_target =0;
    
    for i = 1: npopulationSize
        for j = 1:nTargetPoints
                x1 = targetpoints(j,1);
                y1 = targetpoints(j,2);
                z1 = targetpoints(j,3);
                x2 = positions(i,1);
                y2 = positions(i,2);
                z2 = positions(i,3);
            distance = sqrt( ( (x1-x2)*(x1-x2)  ) + (  (y2-y1)*(y2-y1)  )+ (  (z1-z2)*(z1-z2)  )  );
            if distance <= sensingRange
                number_of_node_covering_target = number_of_node_covering_target  +1;
            end
        end
    end
    f4 =reward*number_of_node_covering_target /npopulationSize;
end



