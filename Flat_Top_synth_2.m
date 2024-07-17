global N;
N = 10;

F_carrier = 2.4e9;
lambda =  physconst('LightSpeed')/F_carrier;
k = 2*pi/lambda;
Dist_fraction = 1;
Dist = Dist_fraction*lambda/2; % distance between sensors
R = Dist/sin(pi/N);

m_direc = 0.7;
fun = @(x,y) power((1+sin(x)),m_direc) .* power((1+cos(y)),m_direc) .* sin(x);
Directivity = (pi .* power(2, 2*m_direc+2)) ./ integral2(fun,0,pi,0,2*pi);
Element_azimuth = 2*pi*(0:N-1)'/N; %pi/2+(0:N-1)'*pi/N pi+zeros(N,1) 2*pi*(0:N-1)'/N
Element_azimuth_patter = Element_azimuth; %Element_azimuth; 

step_manifold = 1;
low_azimuth_manifold = 0;
upper_azimuth_manifold = 359;
low_elevation_manifold = 90;
upper_elevation_manifold = 90;
azimuth_manifold = low_azimuth_manifold:step_manifold:upper_azimuth_manifold;
elevation_manifold = low_elevation_manifold:step_manifold:upper_elevation_manifold;
manifold_len = length(azimuth_manifold) * length(elevation_manifold);

global a;
a = zeros(N,manifold_len);
i = 1;
for el = deg2rad(elevation_manifold)
    for az = deg2rad(azimuth_manifold)
        U_direct = (Directivity/2^(2*m_direc))*(( 1 + sin(el) )^m_direc) * (( 1 + cos(az - Element_azimuth_patter) ).^m_direc);
        % a(:,i) = U_direct .* exp(1j*(-(N-1)/2 : (N-1)/2)'.*(2 * pi * Dist' / lambda * sin(az)));
        a(:,i) = U_direct.*exp(1j*k*R*sin(el*ones(N,1)).*cos(az*ones(N,1)-Element_azimuth));
        i = i + 1;
    end
end


i = 1;
global az_mask; 
az_mask = 0;
for az = deg2rad(azimuth_manifold)
    az_mask(i) = amplitude_mask_butter(az);%amplitude_mask_flat amplitude_mask_butter
    i = i + 1;
end
% az_mask = az_mask/max(az_mask);
el_mask = elevation_manifold;

% W = ones(N,1);
% a_formed = W'*a;
% 
% e_pattern = a_formed - az_mask;
% e_pattern_cost = sqrt(sum(e_pattern.^2) / length(e_pattern));


% phase_1 = optimvar("phase_1","LowerBound",0,"UpperBound",2*pi);
% phase_2 = optimvar("phase_2","LowerBound",0,"UpperBound",2*pi);
% gain_1 = optimvar("gain_1","LowerBound",0,"UpperBound",2);
% gain_2 = optimvar("gain_2","LowerBound",0,"UpperBound",2);
% 
% prob = optimproblem("Objective",e_pattern_cost_f(gain_1, phase_1, gain_2, phase_2));
% 
% [sol,fval] = solve(prob,"Solver","ga","Options",options);

numb_of_variables = 2*N; %  + N + 1
lb = [zeros(1,N) -2*pi*ones(1,N)]; % 2*ones(1,N)  4.7
ub = [3*ones(1,N) 2*pi*ones(1,N)]; % 5*ones(1,N)  10 N amplitudes, N phases, distance in lambda/2 and m_direc
% x0 = x; %ga(@e_pattern_cost_f,numb_of_variables,[],[],[],[],lb,ub); % [rand(1,N) 2*pi*rand(1,N) 1+rand(1,N)]
% 
% x = fmincon(@e_pattern_cost_f,x0,[],[],[],[],lb,ub);
x = ga(@e_pattern_cost_f,numb_of_variables,[],[],[],[],lb,ub,[],optimoptions('ga','PlotFcn', @gaplotbestf)); %,optimoptions('ga','InitialPopulationMatrix',x)
% x = particleswarm(@e_pattern_cost_f,numb_of_variables,lb,ub,optimoptions('particleswarm','InitialSwarmMatrix',x));

weight = x(1:N) .* exp(1j*x(N+1:2*N));
a_formed = abs( weight * a );

function e = e_pattern_cost_f(x_weight)
    
    global a;
    global az_mask;
    N = size(a,1);
    manifold_len = size(a,2);
    weight = x_weight(1:N) .* exp(1j*x_weight(N+1:2*N));
    
    a_formed = abs( weight * a ); 
    e_pattern = a_formed/max(a_formed) - az_mask;
    e = sqrt( e_pattern * e_pattern' / manifold_len );
    
end

function a = get_linear_array_manifold(N, Dist_fraction, m_direc)
    arguments
        N (1,1) double = 10;
        Dist_fraction (1,:) double =  1*ones(1,N); % [2.0785    2.7747    3.1111    2.9061    2.7678    2.8738    2.9607    2.8655    2.9271    3.2293]
        m_direc (1,1) double = 4.7;
    end

    if isempty(Dist_fraction)
        Dist_fraction =  1*ones(1,N);
    end

    % create array manifold
    
    F_carrier = 2.4e9;
    lambda =  physconst('LightSpeed')/F_carrier;
%     k = 2*pi/lambda;
    Dist = Dist_fraction*lambda/2; % distance between sensors
    
%     m_direc = 4.7;
    fun = @(x,y) power((1+sin(x)),m_direc) .* power((1+cos(y)),m_direc) .* sin(x);
    Directivity = (pi .* power(2, 2*m_direc+2)) ./ integral2(fun,0,pi,0,2*pi);
    Element_azimuth = pi+zeros(N,1); % 2*pi*(0:N-1)'/N
    Element_azimuth_patter = Element_azimuth; %Element_azimuth; 
    
    step_manifold = 1;
    low_azimuth_manifold = 0;
    upper_azimuth_manifold = 359;
    low_elevation_manifold = 90;
    upper_elevation_manifold = 90;
    azimuth_manifold = low_azimuth_manifold:step_manifold:upper_azimuth_manifold;
    elevation_manifold = low_elevation_manifold:step_manifold:upper_elevation_manifold;
    manifold_len = length(azimuth_manifold) * length(elevation_manifold);
    
    a = zeros(N,manifold_len);
    i = 1;
    for el = deg2rad(elevation_manifold)
        for az = deg2rad(azimuth_manifold)
            U_direct = (Directivity/2^(2*m_direc))*(( 1 + sin(el) )^m_direc) * (( 1 + cos(az - Element_azimuth_patter) ).^m_direc);
            a(:,i) = U_direct .* exp(1j*(-(N-1)/2 : (N-1)/2)'.*(2 * pi * Dist' / lambda * sin(az)));
            i = i + 1;
        end
    end

end

function ampl = amplitude_mask_flat(theta)
    
    if theta >= deg2rad(170) && theta <= deg2rad(190)
        ampl = db2mag(20);
    elseif (theta >= deg2rad(90) && theta < deg2rad(170)) || (theta > deg2rad(190) && theta <= deg2rad(270))
        ampl = db2mag(15);
    else
        ampl = db2mag(-20);
    end
    
end

function ampl = amplitude_mask_butter(theta)
    theta_p = pi; % думаю, что это частота среза
    theta_center = pi; % думаю, что это центр второго элемента
    n = 2;
    epsilon = 0.5;
    ampl = 1 ./ sqrt(1 + epsilon^2*((theta-theta_center)./theta_p).^(2*n));
    %1
%     if theta > pi
%         ampl = 1 ./ sqrt(1 + epsilon^2*((2*pi-theta)./theta_p).^(2*n));
%     else
%         ampl = 1 ./ sqrt(1 + epsilon^2*(theta./theta_p).^(2*n));
%     end
    % ampl = 2*ampl;
end
