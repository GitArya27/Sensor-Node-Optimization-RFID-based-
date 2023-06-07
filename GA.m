% Define the problem parameters
nAppropriatePositions = 200; % number of appropriate positions
sensingRange = 50; % sensing range of sensors
communicationRange = 100; % communication range of sensors
maxIteration = 100; % maximum number of iterations
nTargetPoints = 50; % number of target points
bsLocation = [250, 250]; % base station location
mutationProbability = 0.2; % mutation probability

% Define the GA parameters
npopulationSize = 200; % size of the population
population = zeros(npopulationSize,3);
GApopulationsize = 100;
generationSize = zeros(GApopulationsize,1);
fitness = zeros(GApopulationsize,1);
targetpoints = zeros(nTargetPoints,3);
eliteCount = 5; % number of elite solutions to keep
tournamentSize = 5; % size of the tournament for selection
crossoverProbability = 0.7; % probability of crossover
mutationRate = 0.1; % mutation rate

% Initialize the population of candidate solutions


for i = 1:nTargetPoints
    for k = 1:3
        targetpoints(i,k) = randi(500);
    end
end
generation = zeros(npopulationSize,3,GApopulationsize); % Generation of diffrent Population

for i = 1:GApopulationsize
    for j = 1: npopulationSize
        for k = 1:3
            generation(j,k,i)= randi(500);
        end
    end
end


%Initialising other variables


plot_nodes(generation(:,:,1),targetpoints,generation(:,:,1),sensingRange,"initial position");



for i = 1:maxIteration
    % Evaluate the fitness of each solution using the objective functions
    for j = 1:GApopulationsize
        fitness(j,1) = evaluateFitness(generation(:,:,j), nTargetPoints, sensingRange, communicationRange, bsLocation,targetpoints);
    end
    % Sort the population by fitness
    [~, idx] = sort(fitness, 'descend');
    generation = generation(:, :, idx);
    fitness = fitness(idx,1);
    
    % Select the elite solutions
    eliteGeneration = generation(:, :, 1:eliteCount);
    
    % Perform selection using tournament selection
    selectedGeneration = tournamentSelection(generation, fitness, tournamentSize, GApopulationsize - eliteCount,npopulationSize);
    
    % Perform crossover
    offspringGeneration = crossover(selectedGeneration, crossoverProbability,npopulationSize);
    
    % Perform mutation
    mutatedGeneration = mutation(offspringGeneration, mutationRate,fitness,GApopulationsize, nTargetPoints, sensingRange, communicationRange, bsLocation,targetpoints);
    
    % Combine the elite and mutated populations
    generation = cat(3, eliteGeneration, mutatedGeneration);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
final_connection_ratio = number_of_targets_covered/(size(positions,1)*(size(positions,1)-1));
final_overlap_ratio = overlap/popvol;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











% Define the tournament selection function
function selectedGeneration = tournamentSelection(Generation, fitness, tournamentSize, numSelections,npopulationSize)
    selectedGeneration = zeros(npopulationSize,3,numSelections);
    for i = 1:numSelections
        % Choose the individuals for the tournament
        tournamentIndividuals = randperm(size(Generation, 3), tournamentSize);
        
        % Select the best individual from the tournament
        [~, idx] = max(fitness(tournamentIndividuals));
        
        % Add the selected individual to the new population
        selectedGeneration(:, :, i) = Generation(:, :, tournamentIndividuals(idx));
    end
end

% Define the crossover function
function offspringGeneration = crossover(parentGeneration, crossoverProbability,npopulationSize)
    offspringGeneration = parentGeneration;
    for i = 1:2:size(parentGeneration, 3)
        % Perform crossover with a given probability
        if i == size(parentGeneration, 3)
            offspringGeneration(:,:,i) = parentGeneration(:, :, i);
            break;
        end
        if rand < crossoverProbability
            % Choose two parents for crossover
            parent1 = parentGeneration(:, :, i);
            parent2 = parentGeneration(:, :, i+1);
            breakpnt= randi(npopulationSize);
            parent11 = parent1(1:breakpnt,:);
            parent12 = parent2(breakpnt+1:npopulationSize,:);
            parent21 = parent2(1:breakpnt,:);
            parent22 = parent1(breakpnt+1:npopulationSize,:);
            offspringGeneration(:,:,i) = cat(1,parent11,parent12);
            offspringGeneration(:,:,i+1)= cat(1,parent21,parent22);
            % Perform crossover by selecting a random split point and
            % swapping the positions of the two parents at and after that
            % point
           
        end
    end
end

function mutatedGeneration = mutation(generation, mutationRate,fitness,GApopulationSize, nTargetPoints, sensingRange, communicationRange, bsLocation,targetpoints)
    mutatedGeneration = generation;
    npopulationSize = size(generation,1);
    for i = 1:GApopulationSize
        if rand()< mutationRate
            mutantvector = 500*rand(npopulationSize,3);
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
            fitness(i,1) = evaluateFitness(mutantvector, nTargetPoints, sensingRange, communicationRange, bsLocation,targetpoints);
        end
    end
end



%fitness finction to get the accuricy of each solution

% Define the fitness function
function fitness = evaluateFitness(positions, nTargetPoints, sensingRange, communicationRange, bsLocation,targetpoints)
    [f1, f2, f3, f4] = eof(positions, nTargetPoints, sensingRange, communicationRange, bsLocation,targetpoints,size(positions,1));
    w1 = 0.25;
    w2 = 0.25;
    w3 = 0.25;
    w4 = 0.25;
    fitness = w1*f1+w2*f2-w3*f3-w4*f4;
end
% Define the objective functions
function [f1, f2, f3, f4] = eof(positions, nTargetPoints, sensingRange, communicationRange, ~,targetpoints,possize)
        reward = 100;

    % Objective function 1: maximize the coverage of target points by sensors
    number_of_targets_covered =0;
    npopulationSize = possize;
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

    % Objective function 4: minimize numer of sensing node
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





