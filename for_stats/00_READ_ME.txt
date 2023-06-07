**How to Use the functions to gets stats**

1> you can use/ create any file inside the stats folder to record stats
2> the main_stats file is recommended to check stats for yourself

How to use it ?

1. you can use the algorithem as functions (like any other function) in main_stats folder .

2. the functions names are bbo_de_stats , bbo_stats , ga_stats 

3. as in many functions these function(also called methods sometimes) takes many input to set value for yourself

4. the parameters that these function take as input are 
        a>numoftarget- it is the total number of targets you wish to set according to you
        b>napppos - it is the number of sensors/positions for sensors 
        c>w1 - it is the weight value for 1st objective fitness
        d>w2 - it is the weight value for 1st objective fitness
        e>w3 - it is the weight value for 1st objective fitness
        f>w4 - it is the weight value for 1st objective fitness
        g>plotnode - it takes binary input that is either 1 or 0 (if u wish to plot the nodes then 1 else 0 )

5. so these function/methods can be called/used as 
        a> bbo_de_stats( numoftargets , napppos , w1 , w2 , w3 , w4 , plotnode )
        b> bbo_stats( numoftargets , napppos , w1 , w2 , w3 , w4 , plotnode )
        c> ga_stats( numoftargets , napppos , w1 , w2 , w3 , w4 , plotnode )

6. upon calling these function they will return a list with 4 outputs those are 
        a> TCR = Target Coverage Ratio
        b> CR = Connection Ratio
        d> OR = Overlap Ratio
        e> NON = number of sensors discovering any targets/ total number of sensors

7. these can be stored in a list as given
        ---->  [TCR,CR,OR,NON]=bbo_de_stats( numoftargets , napppos , w1 , w2 , w3 , w4 , plotnode )
        ---->  [TCR,CR,OR,NON]=bbo_stats( numoftargets , napppos , w1 , w2 , w3 , w4 , plotnode )
        ---->  [TCR,CR,OR,NON]=ga_stats( numoftargets , napppos , w1 , w2 , w3 , w4 , plotnode )




8. finally to measure the time there is a inbuild function in Matlab that is Tic-Toc method
        ---> you can start the stopwatch with tic as a command/code line
        ---> you can stop the running stopwatch with toc which returns the total time taken to execute the code 
        ---> for example :
                    tic
                        // piece of code 
                    time = toc
            this will help to measure time to execute the piece of code and store in a variable called time

Note : you can run a loop with diffrent weights and store the data in a matrix as have already done or just can change the value and record it 








