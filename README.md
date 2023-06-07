# Sensor-Node-Optimization-RFID-based-

Problem Statement :
Abstract:
Proposed a scheme that depends on Biogeographical based optimization (BBO) and Differential Evolution (DE) algorithms to solve a multi-objective optimization problem using a classical weighted sum approach fitness function that is derived from the combination of the following conflicting objectives:
    1. Maximizing target coverage 
    2. Maximizing network connectivity
    3. Minimizing the number of sensing nodes
    4. Minimizing overlapping between sensors 
    
The system model and problem formulation:

In the WSN model, we represent the target space that will be covered by a number of sensors as a 3-D grid. Sensors can be placed at these grid points and all the grid points (i.e. points of interest - PoIs) in the system need to be covered. Each point is represented by (x, y, z) dimensions. These sensor nodes will be placed into predetermined appropriate points (positions). There are some predefined appropriate positions where the sensor nodes will be placed to monitor the target. 
It is assumed that WSN has the following properties: 
    • Base station, target points, and deployed sensor nodes are static.
    • More than one target point can be sensed by a sensor node.
    • Two nodes are connected if they are within the communication range of each other.
    • When a target is present within the sensing range of a sensor, then we say it is covered.
    • All the sensor nodes are said to be homogeneous as all of them possess an equal amount of sensing and communication range.
    • Different potential positions for deploying sensors are identified randomly to monitor set of target points.
    • The sensors that are deployed in identified locations forward data to the base station directly or via other nodes.
