w1 = [0.2,  0.25,  0.3, 0.25,  0.2  ];
w2 = [0.35, 0.3,  0.3,  0.25,  0.3  ];
w3 = [0.25, 0.25, 0.2,  0.35,  0.2  ];
w4 = [0.2,  0.2,  0.2,  0.15,  0.3  ];
TCR_bbo_de = zeros(5,1); %Target coverage ratio
CR_bbo_de = zeros(5,1);    %Connectivity Ratio
OR_bbo_de = zeros(5,1);    %Overlap Ratio
non_bbo_de = zeros(5,1);
time_bbo_de = zeros(5,1);   %number of selected nodes

TCR_bbo = zeros(5,1); 
CR_bbo = zeros(5,1);    
OR_bbo = zeros(5,1);   
non_bbo = zeros(5,1);
time_bbo = zeros(5,1);   

TCR_ga = zeros(5,1); 
CR_ga = zeros(5,1);    
OR_ga = zeros(5,1);    
non_ga = zeros(5,1);
time_ga = zeros(5,1);   





%for 200 targets and 100 sensors
ntarget = 50;
napppos = 200;
for i = 1:5
    tic
    [TCR,CR,OR,NON]= bbo_de_stats(ntarget,napppos,w1(i),w2(i),w3(i),w4(i),0)
    time = toc
    non_bbo_de(i) = NON;
    TCR_bbo_de(i) = TCR;
    CR_bbo_de(i) = CR;    
    OR_bbo_de(i) = OR;    
    time_bbo_de(i) = time;   
    
    tic
    [TCR,CR,OR,NON]= bbo_stats(ntarget,napppos,w1(i),w2(i),w3(i),w4(i),0)
    time = toc
    non_bbo(i) = NON;
    TCR_bbo(i) = TCR;
    CR_bbo(i) = CR;    
    OR_bbo(i) = OR;    
    time_bbo(i) = time;   
    
    tic
    [TCR,CR,OR,NON]= ga_stats(ntarget,napppos,w1(i),w2(i),w3(i),w4(i),0)
    time = toc
    non_ga(i) = NON;
    TCR_ga(i) = TCR;
    CR_ga(i) = CR;    
    OR_ga(i) = OR;    
    time_ga(i) = time;   
    
end

BBO_DE_DATA = cat(2,TCR_bbo_de,CR_bbo_de);
BBO_DE_DATA = cat(2,BBO_DE_DATA ,OR_bbo_de);
BBO_DE_DATA = cat(2,BBO_DE_DATA ,non_bbo_de);
BBO_DE_DATA = cat(2,BBO_DE_DATA ,time_bbo_de);



BBO_DATA = cat(2,TCR_bbo,CR_bbo);
BBO_DATA = cat(2,BBO_DATA ,OR_bbo);
BBO_DATA = cat(2,BBO_DATA ,non_bbo);
BBO_DATA = cat(2,BBO_DATA ,time_bbo);



GA_DATA = cat(2,TCR_ga,CR_ga);
GA_DATA = cat(2,GA_DATA ,OR_ga);
GA_DATA = cat(2,GA_DATA ,non_ga);
GA_DATA = cat(2,GA_DATA ,time_ga);





