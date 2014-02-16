% MCM 2014
% Routine to simulate traffic flow in one lane

% N = number of cars that we want on the road
% T = total time to simulate

% RLpos = RL positions over time
% RLvel = RL velocities over time
% LLpos = LL positions over time
% LLvel = LL velocities over time

%   mph   to   m/s
%	40         17.9
%	45         20.1
%	50         22.4
%	55         24.6
%	60         26.8
%	65         29.1
%	70         31.3
%	75         33.5
%	80         35.8
%	85         38.0
%	90         40.2
%	95         42.5
%	100        44.7

% 1 mph = 0.447 m/s

% IMPORTANT: When there are less than N cars in either lane,
%   we push them all to the front of each vector.
%   The rest of the vector is zero-filled.

% ERROR CODES:
%   If err != 0, err represents the time at which a crash or other event
%   happened.

function [RLpos,RLvel,LLpos,LLvel,err,vPrefL,vPrefR,kR,LR,swTime,phK,phD] = traffic_beast_ODEs(N, T)
    mu_dry = 0.7;
    g_accel = 9.8;
    err = 0;
    del_t = 0.05; % 0.01 seconds per time step
    
    swTime = [0; 0];
    phK = 0;
    phD = 0;
    
    N = 80;
    NR = 40;
    NL = N - NR;
%--------------------------------------------------------------------------    
    % initialize right-lane vectors
    kR = zeros(N,1);
    xR = 10^20 * ones(N,1);
    vR = kR;
    vPrefR = kR;
    buffEnterCountR = kR;
    dR = kR;
    LR = kR; % little l in the model
    buffDistR = kR;
    swRtoL = kR; % keeps track of which cars want to switch lanes
%--------------------------------------------------------------------------    
    % initialize right lane output matrices
    RLpos = zeros(N+1, (T * del_t)+1);
    RLvel = zeros(N+1, (T * del_t)+1);
%--------------------------------------------------------------------------
    %initialize left lane vectors
    kL = zeros(N,1);
    xL = 10^20 * ones(N,1);
    vL = kL;
    vPrefL = kL;
    buffEnterCountL = kL;
    dL = kL;
    LL = kL; % little l in the model
    buffDistL = kL;
    swLtoR = kL; % keeps track of which cars want to switch lanes

%--------------------------------------------------------------------------
    % initialize left lane output matrices
    LLpos = zeros(N+1, (T * del_t)+1);
    LLvel = zeros(N+1, (T * del_t)+1);
%--------------------------------------------------------------------------
    % generate cars
    [xR,vR,vPrefR,xL,vL,vPrefL] = gen_cars(NR,NL);
    
    % insert initial state into output matrices
    RLpos(:,1) = [0;xR];
    RLvel(:,1) = [0;vR];
    LLpos(:,1) = [0;xL];
    LLvel(:,1) = [0;vL];
    
%--------------------------------------------------------------------------
    % MAIN LOOP for simulation
    counter = 1;
    for t = del_t:del_t:T
        counter = counter + 1;
        % calculate k buffers
        [LR,buffDistR,kR,LL,buffDistL,kL] = calc_K(NR,NL,mu_dry,g_accel,vR,vL);
        
        % calculate current distances between cars
        [dR,dL] = calc_D(xR,xL,NR,NL,dR,dL);
        
        % iterate through and perform any lane changes R to L
        for j1 = 1:NR-1
            if(LR(j1) < dR(j1) && dR(j1) < kR(j1))
                % check for crashes too in switch function
                [N,NR,NL,xR,vR,LR,buffDistR,kR,dR,xL,vL,LL,buffDistL,kL,dL,mu_dry,g_accel,swTime,buffEnterCountR,vPrefR,vPrefL] = switch_RtoL(N,NR,NL,xR,xL,vR,vL,swRtoL,dR,dL,kR,kL,LR,LL,mu_dry,g_accel,buffEnterCountR,buffDistR,buffDistL,t,vPrefR,vPrefL,swTime);
            end
        end
        
        % this function needs to update NR and NL:
        [dR,dL] = calc_D(xR,xL,NR,NL,dR,dL);
        [LR,buffDistR,kR,LL,buffDistL,kL] = calc_K(NR,NL,mu_dry,g_accel,vR,vL);
        
        % iterate through and perform any lane changes L to R
        for j2 = 1:NL
            [phK,phD,N,NR,NL,xR,vR,LR,buffDistR,kR,dR,xL,vL,LL,buffDistL,kL,dL,mu_dry,g_accel,swTime,buffEnterCountL,vPrefR,vPrefL] = switch_LtoR(phK,phD,N,NR,NL,xR,xL,vR,vL,swLtoR,dR,dL,kR,kL,LR,LL,mu_dry,g_accel,buffEnterCountL,buffDistR,buffDistL,t,vPrefR,vPrefL,swTime);
        end
        
        % this function needs to update NR and NL:
        [dR,dL] = calc_D(xR,xL,NR,NL,dR,dL);
        [LR,buffDistR,kR,LL,buffDistL,kL] = calc_K(NR,NL,mu_dry,g_accel,vR,vL);
        
        % determine which cars are in the k buffers and adjust velocities
        for j3 = 1:NR-1
            if(dR(j3) > kR(j3))
                tol = 0.01 * vPrefR(j3);
                % outside of the buffer region
                % adjust speed
                if(vR(j3) > vPrefR(j3) + tol || vR(j3) < vPrefR(j3) - tol)
                    vR(j3) = vR(j3) + 1/3*(vPrefR(j3) - vR(j3));
                end
            
            elseif(dR(j3) <= kR(j3))
                % in the buffer region where ODE is valid
                % REMEMBER TO RESET TO ZERO WHEN LANE CHANGE HAPPENS
                buffEnterCountR(j3) = buffEnterCountR(j3) + 1;
                vR(j3) = vR(j3+1)+(vR(j3)-vR(j3+1))*exp(-2*mu_dry*g_accel*del_t*buffEnterCountR(j3)/(3*(vR(j3)+vR(j3+1))));
            %elseif(dR(j3) < LR(j3)/50)
            %    err = t;
            %    dR(j3)
            %    LR(j3)
                % exit right loop
            %    break;
            %elseif(dR(j3) < LR(j3))
            %    vR(j3) = vR(j3+1);
            else
                err = -t;
                % other unknown error occurred
            end
        end
        
        if(err ~= 0)
            % exit main loop
            break;
        end
        
        % determine which cars are in the k buffers and adjust velocities
        for j4 = 1:NL-1
            if(dL(j4) > kL(j4))
                tol = 0.01 * vPrefL(j4);
                % outside of the buffer region
                % adjust speed
                if(vL(j4) > vPrefL(j4) + tol || vL(j4) < vPrefL(j4) - tol)
                    vL(j4) = vL(j4) + 1/3*(vPrefL(j4) - vL(j4));
                end
            
            elseif(dL(j4) <= kL(j4))
                % in the buffer region where ODE is valid
                 % REMEMBER TO RESET TO ZERO WHEN LANE CHANGE HAPPENS
                buffEnterCountL(j4) = buffEnterCountL(j4) + 1;
                vL(j4) = vL(j4+1)+(vL(j4)-vL(j4+1))*exp(-2*mu_dry*g_accel*del_t*buffEnterCountL(j4)/(3*(vL(j4)+vL(j4+1))));
                
            %elseif(dL(j4) < LL(j4)/50)
            %    err = t;
                % exit right loop
            %    break;
            %elseif(dL(j4) < LL(j4))
            %    vL(j4) = vL(j4+1);
            else
                err = -t;
                % other unknown error occurred
            end
        end
        
        if(err ~= 0)
            % exit main loop
            break;
        end
        
        for j5 = 1:N
            if(j5 <= NR)
                xR(j5) = xR(j5) + vR(j5)*del_t;
            end
            if(j5 <= NL)
                xL(j5) = xL(j5) + vL(j5)*del_t;
            end
        end
        
        RLpos(:,counter) = [t;xR];
        RLvel(:,counter) = [t;vR];
        LLpos(:,counter) = [t;xL];
        LLvel(:,counter) = [t;vL];
        
    end
end

function [xR,vR,vPrefR,xL,vL,vPrefL] = gen_cars(NR,NL)
    % preset for now
    numcars = NR + NL;
    
  %  xR = [10; 50; 150; 200; 244; 300; 350; 400; 450; 500; 550; 627; 670; 750; 800; 900; 1200; 1600; 1700; 1800; 1900; 2000; 2100; 2300; 2400; 2500; 2650; 2770; 2800; 2950;3000;3100;3200;3300;3400;3500;3600;3700;3800;3900;4000;4100;4200;4300] ;
    %vR = [26.1; 33; 31; 20; 40; 27; 45; 20; 34; 31; 32; 45; 47; 30; 32; 43; 54; 45; 56; 12;32;34;25;43;45;54;43;23;32;12;34;45;43;57;37;12;43;34;52;18];
  
  % vR = normrnd(30, 5,[N 1]);
   
   xR(1) = 10;
   for helpGeneration = 2:numcars
       if(helpGeneration <= NR)
            xR = [xR; 30*helpGeneration];
       else
           xR = [xR; 10^20];
       end
   end
    vR = normrnd(30, 5, [NR, 1]);
    for help3 = 1:NR
       if(vR(help3) < 10)
           vR(help3) = 10;
       end
    end
   for help4 = NR+1:numcars
       vR = [vR; 0];
   end
   vPrefR = vR;
    
   xL(1) = 23;
   for Lhelper = 2:numcars
       if(Lhelper <= NL)
          xL = [xL; 25*Lhelper];
       else
           xL = [xL; 10^20];
       end
   end
   vL = normrnd(35, 5, [NL, 1]);
   for help = 1:NL
       if(vL(help) < 10)
           vL(help) = 10;
       end
   end
   for help5 = NL+1:numcars
       vL = [vL; 0];
   end
   vPrefL = vL;
end


function [LR,buffDistR,kR,LL,buffDistL,kL] = calc_K(NR,NL,mu_dry,g_accel,vR,vL)
    % add some probability factor to change lR, lL
    pfactor = 3;
    
    LR = zeros(NR+NL, 1);
    buffDistR = LR;
    kR = LR;
    LL = LR;
    buffDistL = LR;
    kL = LL;
    
    for ylll = 1:NR-1
        LR(ylll) = 2*vR(ylll+1);
        buffDistR(ylll) = ((vR(ylll))^2 - (vR(ylll+1))^2)/(2*mu_dry*g_accel);
        kR(ylll) = LR(ylll) + pfactor*buffDistR(ylll); 
    end
    if(NR ~= 0)
        LR(NR) = 0;
        buffDistR(NR) = 0;
        kR(NR) = 0;
    end 

    for icalcK2 = 1:NL-1
        LL(icalcK2) = 2*vL(icalcK2+1);
        buffDistL(icalcK2) = ((vL(icalcK2))^2 - (vL(icalcK2+1))^2)/(2*mu_dry*g_accel);
        kL(icalcK2) = LL(icalcK2) + pfactor*buffDistL(icalcK2); 
    end
    if(NL ~= 0)
        LL(NL) = 0;
        buffDistL(NL) = 0;
        kL(NL) = 0;
    end
end

function [dR,dL] = calc_D(xR,xL,NR,NL,dR,dL)
    for icalcD = 1:NR-1
        dR(icalcD) = xR(icalcD+1) - xR(icalcD);
    end
    if(NR ~= 0)
        dR(NR) = 0;
    end
    
    for icalcD2 = 1:NL-1
        dL(icalcD2) = xL(icalcD2+1) - xL(icalcD2);
    end
    if(NL ~= 0)
        dL(NL) ~= 0;
    end
end


function [N,NR,NL,xR,vR,LR,buffDistR,kR,dR,xL,vL,LL,buffDistL,kL,dL,mu_dry,g_accel,swTime,buffEnterCountR,vPrefR,vPrefL] = switch_RtoL(N,NR,NL,xR,xL,vR,vL,swRtoL,dR,dL,kR,kL,LR,LL,mu_dry,g_accel,buffEnterCountR,buffDistR,buffDistL,t,vPrefR,vPrefL,swTime)
    % check to see if we can switch safely
    
    for iSEE1 = 1:NR
        swRtoL(iSEE1) = 0;
    end
    
    for iSEE = 1:NR-1
        % figure out which cars might want to switch i.e. are in k buffer
        if(dR(iSEE) < kR(iSEE))
            %indicator vector
            swRtoL(iSEE) = 1;
        
            % send probe and search
            if(NL ~= 0 && xR(iSEE) < 10^20)
                index = 1;
                while(index <= NL)
                    if(xR(iSEE) < xL(index))
                        break;
                    end
                    index = index + 1;
                end
            
                % now index is pointing to the first car in LL ahead of
                % RL car or index is equal to NL+1
                if(index == 1)
                    % no cars in the left lane are behind xR(iSEE)
                    tempL = 2*vL(index);
                    if(xR(iSEE) + 2*tempL < xL(index))
                        % can switch without a crash
                        % first, insert the car from the right system into the
                        % left system at the end of the vectors
                        swTime = [swTime, [t;0]];
                        NL = NL + 1;
                        NR = NR - 1;
                        xL(NL) = xR(iSEE);
                        xR(iSEE) = 10^20;
                        vL(NL) = vR(iSEE);
                        vR(iSEE) = 0;
                        vPrefL(NL) = vPrefR(iSEE);
                        vPrefR(iSEE) = 0;
                        buffEnterCountR(iSEE) = 0;
                        swRtoL(iSEE) = 0;
                        
                        % reorder all vectors appropriately
                        xR = shiftDown(xR,iSEE,N);
                        vR = shiftDown(vR,iSEE,N);
                        vPrefR = shiftDown(vPrefR,iSEE,N);
                        buffEnterCountR = shiftDown(buffEnterCountR,iSEE,N);
                        [xL,vL,vPrefL] = reSort(xL,vL,vPrefL,NL);
                        
                        % recalc d and k for the entire system
                        [dR,dL] = calc_D(xR,xL,NR,NL,dR,dL);
                        [LR,buffDistR,kR,LL,buffDistL,kL] = calc_K(NR,NL,mu_dry,g_accel,vR,vL);
                        
                    end
                else
                    % first, check that you won't crash into the car in LL
                    % at index
                    tempL = 2*vL(index);
                    if(xR(iSEE) + 2*tempL < xL(index))
                        % now, check ratios of phantom buffers to see if we
                        % will switch (check with car behind at index-1)
                        pfactor = 3;
                        phantomL = 2*vR(iSEE);
                        phantomBuff = ((vL(index-1))^2 - (vR(iSEE))^2)/(2*mu_dry*g_accel);
                        phantomK = phantomL + pfactor*phantomBuff;
                        phantomD = xR(iSEE) - xL(index-1);
                        
                        if(phantomD > 2/3*phantomK && phantomD > phantomL)
                            % switch lanes
                            swTime = [swTime, [t;0]];
                            NL = NL + 1;
                            NR = NR - 1;
                            xL(index) = xR(iSEE);
                            xR(iSEE) = 10^20;
                            vL(index) = vR(iSEE);
                            vR(iSEE) = 0;
                            vPrefL(NL) = vPrefR(iSEE);
                            vPrefR(iSEE) = 0;
                            buffEnterCountR(iSEE) = 0;
                            swRtoL(iSEE) = 0;
                            
                            % reorder all vectors appropriately
                            xR = shiftDown(xR,iSEE,N);
                            vR = shiftDown(vR,iSEE,N);
                            vPrefR = shiftDown(vPrefR,iSEE,N);
                            buffEnterCountR = shiftDown(buffEnterCountR,iSEE,N);
                            [xL,vL,vPrefL] = reSort(xL,vL,vPrefL,NL);
                            
                            % recalc d and k for the entire system
                            [dR,dL] = calc_D(xR,xL,NR,NL,dR,dL);
                            [LR,buffDistR,kR,LL,buffDistL,kL] = calc_K(NR,NL,mu_dry,g_accel,vR,vL);
                        end
                    end
                end
            else
                if(xR(iSEE) >= 10^19)
                    break;
                end
                swTime = [swTime, [t;0]];
                % can switch because there are no cars in left lane
                NL = NL + 1;
                NR = NR - 1;
                xL(NL) = xR(iSEE);
                xR(iSEE) = 10^20;
                vL(NL) = vR(iSEE);
                vR(iSEE) = 0;
                vPrefL(NL) = vPrefR(iSEE);
                vPrefR(iSEE) = 0;
                buffEnterCountR(iSEE) = 0;
                swRtoL(iSEE) = 0;
                
                % reorder all vectors appropriately
                xR = shiftDown(xR,iSEE,N);
                vR = shiftDown(vR,iSEE,N);
                vPrefR = shiftDown(vPrefR,iSEE,N);
                buffEnterCountR = shiftDown(buffEnterCountR,iSEE,N);
                [xL,vL,vPrefL] = reSort(xL,vL,vPrefL,NL);
                
                % recalc d and k for the entire system
                [dR,dL] = calc_D(xR,xL,NR,NL,dR,dL);
                [LR,buffDistR,kR,LL,buffDistL,kL] = calc_K(NR,NL,mu_dry,g_accel,vR,vL);
            end
        end
    end
end

function [xVec,vVec,vPrefVec] = reSort(xVec,vVec,vPrefVec,nonZeroInd)
    if(nonZeroInd > 1)
        for rSrS = nonZeroInd:-1:2
            if(xVec(rSrS) < xVec(rSrS-1))
                tmp = xVec(rSrS-1);
                xVec(rSrS-1) = xVec(rSrS);
                xVec(rSrS) = tmp;
                
                tmp = vVec(rSrS-1);
                vVec(rSrS-1) = vVec(rSrS);
                vVec(rSrS) = tmp;
                
                tmp = vPrefVec(rSrS-1);
                vPrefVec(rSrS-1) = vPrefVec(rSrS);
                vPrefVec(rSrS) = tmp;
            end
        end
    end
end

function [vec] = shiftDown(vec, index, size)
    if(index == size)
        % no need to reorder
        % for robustness, set the index val to zero anyway
        vec(index) = 0;
    else
        % shift all values one up
        for Ordi = index:size-1
            temp_sw = vec(Ordi+1);
            vec(Ordi+1) = vec(Ordi);
            vec(Ordi) = temp_sw;
        end
    end
end

function [phK,phD,N,NR,NL,xR,vR,LR,buffDistR,kR,dR,xL,vL,LL,buffDistL,kL,dL,mu_dry,g_accel,swTime,buffEnterCountL,vPrefR,vPrefL] = switch_LtoR(phK,phD,N,NR,NL,xR,xL,vR,vL,swLtoR,dR,dL,kR,kL,LR,LL,mu_dry,g_accel,buffEnterCountL,buffDistR,buffDistL,t,vPrefR,vPrefL,swTime)
    % check to see if we can switch safely
    
    for iSEE3 = 1:NL
        swLtoR(iSEE3) = 0;
    end
    
    for FOO = 1:NL
        % figure out which cars might want to switch: For left to right
        % switching, we just check buffers and phantom buffers and make
        % sure that no buffers would interfere if we switched.
        
        if(NR ~= 0 && xL(FOO) < 10^20)
            ind = 1;
            while(ind <= NR)
                if(xL(FOO) < xR(ind))
                    break;
                end
                ind = ind + 1;
            end
            
            % if ind is at 1, car FOO is behind all cars in the right lane           
            % this case really shouldn't happen much
            % check phantom buffers to see if we will switch
            if(ind == 1)
                pfactor = 3;
                phantL = 2*vR(ind); 
                phantBuff = ((vL(FOO))^2 - (vR(ind))^2)/(2*mu_dry*g_accel);
                % mult by 5 for an extra buffer for changing lanes
                phantK = phantL + pfactor*phantBuff;% + 5*phantL;
                phantD = xR(ind) - xL(FOO);
                
                if(phantD > phantK)
                    % switch
                    swTime = [swTime, [0;t]];
                    NR = NR + 1;
                    NL = NL - 1;
                    if(xL(FOO) == 0)
                        fprintf(1, 'WTF');
                    end
                    xR(NR) = xL(FOO);
                    xL(FOO) = 10^20;
                    vR(NR) = vL(FOO);
                    vL(FOO) = 0;
                    vPrefR(NR) = vPrefL(FOO);
                    vPrefL(FOO) = 0;
                    buffEnterCountL(FOO) = 0;
                    swLtoR(FOO) = 0;
                
                    % reorder all vectors appropriately
                    xL = shiftDown(xL,FOO,N);
                    vL = shiftDown(vL,FOO,N);
                    vPrefL = shiftDown(vPrefL,FOO,N);
                    buffEnterCountL = shiftDown(buffEnterCountL,FOO,N);
                    [xR,vR,vPrefR] = reSort(xR,vR,vPrefR,NR);
                
                    % recalc d and k for the entire system
                    [dR,dL] = calc_D(xR,xL,NR,NL,dR,dL);
                    [LR,buffDistR,kR,LL,buffDistL,kL] = calc_K(NR,NL,mu_dry,g_accel,vR,vL);
                end
            else
                % check phantom buffers generated by xL(FOO) and then check
                % phantom buffers generated by xR(ind)
                pfactor = 3;
                phanL = 2*vL(FOO);
                phanBuff = ((vR(ind-1))^2 - (vL(FOO))^2)/(2*mu_dry*g_accel);
                phanK = phanL + pfactor*phanBuff;
                phanD = xL(FOO) - xR(ind-1);
                
                if(phanD > phanK)
                    if(xL(FOO) < xR(ind))
                        phL = 2*vR(ind);
                        phBuff = ((vL(FOO))^2 - (vR(ind))^2)/(2*mu_dry*g_accel);
                        phK = phL + pfactor*phBuff;% + 5*phL;
                        phD = xR(ind) - xL(FOO);
                    
                        if( phD > phK )
                            % switch
                            swTime = [swTime, [0;t]];
                            NR = NR + 1;
                            NL = NL - 1;
                            xR(NR) = xL(FOO);
                            xL(FOO) = 10^20;
                            vR(NR) = vL(FOO);
                            vL(FOO) = 0;
                            vPrefR(NR) = vPrefL(FOO);
                            vPrefL(FOO) = 0;
                            buffEnterCountL(FOO) = 0;
                            swLtoR(FOO) = 0;
                
                            % reorder all vectors appropriately
                            xL = shiftDown(xL,FOO,N);
                            vL = shiftDown(vL,FOO,N);
                            vPrefL = shiftDown(vL,FOO,N);
                            buffEnterCountL = shiftDown(buffEnterCountL,FOO,N);
                            [xR,vR,vPrefR] = reSort(xR,vR,vPrefR,NR);
                
                            % recalc d and k for the entire system
                            [dR,dL] = calc_D(xR,xL,NR,NL,dR,dL);
                            [LR,buffDistR,kR,LL,buffDistL,kL] = calc_K(NR,NL,mu_dry,g_accel,vR,vL);
                        end
                    else
                        phL = 2*vL(FOO);
                        phBuff = ((vR(ind))^2 - (vL(FOO))^2)/(2*mu_dry*g_accel);
                        phK = phL + pfactor*phBuff;% + 5*phL;
                        phD = xL(FOO) - xR(ind);
                    
                        if( phD > phK )
                            % switch
                            swTime = [swTime, [0;t]];
                            NR = NR + 1;
                            NL = NL - 1;
                            xR(NR) = xL(FOO);
                            xL(FOO) = 10^20;
                            vR(NR) = vL(FOO);
                            vL(FOO) = 0;
                            vPrefR(NR) = vPrefL(FOO);
                            vPrefL(FOO) = 0;
                            buffEnterCountL(FOO) = 0;
                            swLtoR(FOO) = 0;
                
                            % reorder all vectors appropriately
                            xL = shiftDown(xL,FOO,N);
                            vL = shiftDown(vL,FOO,N);
                            vPrefL = shiftDown(vL,FOO,N);
                            buffEnterCountL = shiftDown(buffEnterCountL,FOO,N);
                            [xR,vR,vPrefR] = reSort(xR,vR,vPrefR,NR);
                
                            % recalc d and k for the entire system
                            [dR,dL] = calc_D(xR,xL,NR,NL,dR,dL);
                            [LR,buffDistR,kR,LL,buffDistL,kL] = calc_K(NR,NL,mu_dry,g_accel,vR,vL);
                        end
                    end
                end
            end
        end
    end
end









