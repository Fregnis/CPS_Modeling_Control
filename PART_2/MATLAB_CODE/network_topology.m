function [N,Ad,Aug,D1,L,G,Gdiag] = network_topology(number)
%NETWORK_TOPOLOGY Select network topology

switch (number)
    case 1
        % Network topology n.1
        % THOSE ARE ALL SELF LOOPS, all values on the diag of the Ad matrix
        % The leader is connected to the first node but it is not connected to any
        % other node

        % Converge to 0
        % DOESNT converge to constant
        % DOESNT converge to ramp
        % DOESNT converge to sin
        N = 6;

        Ad = [0 0 0 0 0 0
              0 3 0 0 0 0
              0 0 1 0 0 0
              0 0 0 4 0 0
              0 0 0 0 2 0
              0 0 0 0 0 6]; %matrice delle adiacenze

        Aug = diag([1 1 3 1 4 2 6]);

        D1 = diag(ones(1,6));

        L = D1 - Ad;

        G = diag([1,0,0,0,0,0]);
        Gdiag = [1 0 0 0 0 0];
    case 2
        % Network topology n.2
        % All nodes connected in a line, leader is only connected to first node

        % DISTRIBUTED
        % Convergence to constant
        % Convergence to ramp
        % Convergence to sin

        % LOCAL
        % Convergence to constant
        % Convergence to ramp
        % Convergence to sin
        N = 6;

        Ad = [ 0 0 0 0 0 0
                2 0 0 0 0 0
                0 6 0 0 0 0
                0 0 1 0 0 0
                0 0 0 1 0 0
                0 0 0 0 3 0];

        Aug = diag([ 1,2,6,1,1,3 ]);

        D1 = diag(ones(1,6));

        L = D1 - Ad;

        G = diag([1,0,0,0,0,0]);
        Gdiag = [1 0 0 0 0 0];
    case 3
        % Network topology n.3
        % Leader nodes connected with two groups of follower nodes
        % In each one of the two groups the nodes are connected in a line
        % Representation, '->' represnts a directed connection 
        %          -> 1 -> 2 -> 3
        % Leader / 
        %        \
        %          -> 4 -> 5 -> 6

        % DISTRIBUTED
        % Convergence to constant
        % Convergence to ramp
        % Convergence to sin

        % LOCAL
        % Convergence to constant
        % Convergence to ramp
        % Convergence to sin
        N = 6;

        Ad = [ 0 0 0 0 0 0
                2 0 0 0 0 0
                0 1 0 0 0 0
                0 0 0 0 0 0
                0 0 0 1 0 0
                0 0 0 0 3 0];

        Aug = [[1 0 0 1 0 0]', Ad];

        D1 = diag(ones(1,6));

        L = D1 - Ad;

        G = diag([1,0,0,1,0,0]);
        Gdiag = [1 0 0 1 0 0];
    case 4
        % Network topology n.4
        % Leader nodes connected with two groups of follower nodes
        % In each one of the two groups the nodes are connected in a line
        % Each first node of the two lines is connected with the last node on the
        % other line
        % Representation, '->' represnts a directed connection 
        % Leader -> 1 -> 2 -> 3   
        %           |
        %           4 -> 5 -> 6
        % and 6 connected to 1
        % 
        % Convergence is slower probably due to the higher amount of information
        % exchanged BUT it is more robust, it is able to lower the effect of the
        % error
        
        N = 6;

        Ad = [ 0 0 0 0 0 0
                2 0 0 0 0 0
                0 1 0 0 0 0
                2 0 3 0 0 0
                0 0 0 1 0 0
                0 0 0 0 3 0];

        Aug = [[1 0 0 0 0 0]', Ad];

        D1 = diag(ones(1,6));

        L = D1 - Ad;

        G = diag([1,0,0,0,0,0]);
        Gdiag = [1 0 0 0 0 0];
end
        
end

