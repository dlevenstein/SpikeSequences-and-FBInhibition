%% Integration of Hodgkin-Huxley equations w/ Euler method
    clear; clf; 
   
    %% Setting parameters
    % Maximal conductances (in units of mS/cm^2) ; 1=k, 2=Na, 3=R
    g(1)=36; g(2)=120; g(3)=0.3; 
    % reversal potential / "Battery voltage" (in mV); 1=n, 2=m, 3=h
    E(1)=-12; E(2)=115; E(3)=10.613; 

    %% note: when we refer to (1), (2), and (3), whether we are
    % referring to k/Na/R or n/m/h is determined by what unit we are
    % using (mS/cm^2 versus mV)
    %%
    %Initialization of some variables
    I_ext=0; V=-10; x=zeros(1,3); x(3)=1; t_rec=0;
    % I is current: we set external current to 0 (no incoming ions)
    % membrane potential (V) is -10
    
    %Time step for integration
       dt=0.01; % mS
       
    %% Integration with Euler method
    for t= -30:dt:50    %simulate from  t=-30ms to t=50ms in steps of dt
  
        if t==10; I_ext=10; end % turns external current on at t=10
        if t==40; I_ext=0; end % turns external current off at t=40
        
    % alpha functions used by H-H
    Alpha(1)=(10-V)/(100*(exp((10-V)/10)-1));
    %% can we say e.g. Alpha(1) is Alpha(n), because V is a unit of mV? 
    %%
    Alpha(2)=(25-V)/(10*(exp((25-V)/10)-1));
    % Alpha(2) is Alpha(m)?
    Alpha(3)=0.07*exp(-V/20); 
    % Alpha(3) is Alpha(h)?
    
    % beta functions used by H-H
    Beta(1)=0.125*exp(-V/80);
    % as with alpha, can we say Beta(1) is Beta(n)?
    Beta(2)=4*exp(-V/18);
    Beta(3)=1/(exp((30-V)/10)+1);
    
    % tau_x and x_0 (x=1,2,3) are defined with alpha and beta
        tau=1./(Alpha+Beta); 
        x_0=Alpha.*tau;
              %%
        % leaky integration with Euler method
        x=(1-dt./ tau) .* x+dt ./ tau .* x_0;    %DL: element-wise multiplication/division should be ".*" and "./"
        
        % calculate actual conductances g with given n, m, h
        gnmh(1)=g(1)*x(1)^4;
        gnmh(2)=g(2)*x(2)^3*x(3);
        gnmh(3)=g(3);
        
        % Ohm's Law
            I=gnmh .* (V-E);        %DL: element-wise multiplication/division should be ".*" and "./"
            % gnmh, V-E are both 1x3 vectors
       % update voltage of membrane
       V=V+dt*(I_ext-sum(I));
       % record some variables for plotting after equilibration 
            if t>=0;
                t_rec=t_rec+1;
                x_plot(t_rec)=t;
                y_plot(t_rec)=V;
                
            end
            
    end %time loop
    
    %% Plotting results
        plot(x_plot,y_plot); xlabel('Time'); ylabel('Voltage');
        
    