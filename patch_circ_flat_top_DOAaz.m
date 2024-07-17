L = 100;  % Length of signal
azimuth =  75; 
elevation = 90*ones(size(azimuth)); 
% elevation =  25:2:90; 
% azimuth = 0*ones(size(elevation)); 
true_doa = [azimuth' elevation'];
S = 1; % The number of signals

% create array manifold
N = 20;
Rings = 1;
Nsub = 2;
F_carrier = 5e9;
lambda =  physconst('LightSpeed')/F_carrier;
k = 2*pi/lambda;
Dist = 1*lambda/2; % distance between sensors
% R = (sqrt(2)+1)*lambda/4; %(sqrt(2)+1)*lambda/4 - 8-угол (sqrt(3)/2)*Dist-шестиугол
R = 1*Dist/sin(pi/N); % радиус АР
Element_azimuth = (0:N-1)'*2*pi/N; % угол расположения антенных элементов по кольцу
[x,y,z] = pol2cart(Element_azimuth,R,0);
Element_azimuth_patter = Element_azimuth; %[1.5010    1.4835    0.7679    1.6406    3.9270    6.1959]'; 
Ring_Shift = [0 0];
Height = lambda;
m_direc = 2.7;
fun = @(x,y) power((1+sin(x)),m_direc) .* power((1+cos(y)),m_direc) .* sin(x);
Directivity = (pi .* power(2, 2*m_direc+2)) ./ integral2(fun,0,pi,0,2*pi);

patchObject = patchMicrostrip.empty(N,0);
for n = 1:N
    patchObject(n) = design(patchMicrostrip('Substrate',dielectric('Air'), 'Conductor', metal('PEC')), F_carrier);
    patchObject(n).Tilt = [90 rad2deg(Element_azimuth(n))];
    patchObject(n).TiltAxis = [0 1 0; 0 0 1];
end
circArrayObject = conformalArray(ElementPosition=[x y zeros(size(x))],Element=patchObject);
mesh(circArrayObject,MaxEdgeLength=0.09);
S_params = sparameters(circArrayObject,F_carrier);
C_1 = S_params.Parameters;
C_mutual = C_1 + diag(1-abs(diag(C_1)));

load("circ_flat_top_conf_6ae_air.mat")
Wsub = zeros(Nsub,N);
Wsub(1,1:N/Nsub) = linArrayObject.AmplitudeTaper.*exp(1j*deg2rad(linArrayObject.PhaseShift));
Wsub(2,N/Nsub+1:end) = linArrayObject.AmplitudeTaper.*exp(1j*deg2rad(linArrayObject.PhaseShift));

step_manifold = 1/2;
low_azimuth_manifold = 50;
upper_azimuth_manifold = 100;
low_elevation_manifold = 90;
upper_elevation_manifold = 90;

azimuth_manifold = low_azimuth_manifold:step_manifold:upper_azimuth_manifold;
elevation_manifold = low_elevation_manifold:step_manifold:upper_elevation_manifold;
manifold_len = length(azimuth_manifold) * length(elevation_manifold);

a = zeros(N,manifold_len);
i = 1;
for el = elevation_manifold
    for az = azimuth_manifold
        for n = 1:N
            % U_direct = (Directivity/2^(2*m_direc))*(( 1 + sin(deg2rad(el)) )^m_direc) * (( 1 + cos(deg2rad(az) - Element_azimuth(n)) ).^m_direc);
            U_direct = db2mag(pattern(patchObject(n),F_carrier,az,el-90, "Type", "realizedgain"));
            % U_direct = db2mag(pattern(circArrayObject,F_carrier,az,el-90, 'ElementNumber',n, "Type", "realizedgain"));
            a(n,i) = U_direct*exp(1j*k*R*sin(deg2rad(el)).*cos(deg2rad(az)-Element_azimuth(n))); %U_direct.*
        end
        i = i + 1;
    end
end
asub = Wsub*a;

RMSEs = zeros(length(true_doa), 2);
for coordinate = 1:length(true_doa)
    % signals DOA
    A = zeros(N,S);
    for i = 1:S    
        for n = 1:N
            % U_direct = (Directivity/2^(2*m_direc))*(( 1 + sin(deg2rad(elevation(i,coordinate))) )^m_direc) * (( 1 + cos(deg2rad(azimuth(i,coordinate)) - Element_azimuth(n)) ).^m_direc);
            U_direct = db2mag(pattern(circArrayObject,F_carrier,azimuth(i,coordinate),elevation(i,coordinate)-90, 'ElementNumber',n, "Type", "realizedgain"));
            % U_direct = db2mag(pattern(patchObject(n),F_carrier,azimuth(i,coordinate),elevation(i,coordinate)-90, "Type", "realizedgain"));
            A(n,i) =  U_direct*exp(1j*k*R*sin(deg2rad(elevation(i,coordinate))) *cos(deg2rad(azimuth(i,coordinate))-Element_azimuth(n))); %db2mag(pattern(circArrayObject,F_carrier,azimuth(i,coordinate),elevation(i,coordinate)-90, 'ElementNumber',n, "Type", "realizedgain")) *... 
        end
    end

    trials = 100;
    skipped = 0;
    nose_power = db2pow(-5);
    sig_power = 1;
    doa_estimated = zeros(S,2);
    for trial = 1:trials
        % add noise
        n = sqrt(nose_power/2).*(randn([N,L]) + 1j*randn([N,L]));
        x = sqrt(sig_power/2).*(randn([S,L]) + 1j*randn([S,L]));
        % create output signals
        y = Wsub*C_mutual*A*x + Wsub*n; %C_mutual*
        Rrr = (1/L)*(y*y');
        % MUSIC
        [V,D] = eig(Rrr);
        [C,Inx]=sort(diag(D));
        En = V(:,Inx(1:(Nsub-S)));
        P_MUSIC = zeros(1,manifold_len);
        for i = 1:manifold_len
            P_MUSIC(i) = 1/abs((asub(:,i))'*En*(En)'*asub(:,i));
        end
        P_MUSIC = (P_MUSIC/max(P_MUSIC));
        plot(azimuth_manifold,P_MUSIC)
        [~,I_mu] = findpeaks(P_MUSIC,'SORTSTR','descend');
        if isempty(I_mu)
            [~,I_mu] = max(P_MUSIC);
        end
        doa_est_sigs = sort(elevation_manifold(I_mu(1:S)),'ascend');
        doa_estimated = doa_estimated + (true_doa(coordinate,2) - sort(doa_est_sigs, 1)).^2;
    end
    doa_estimated = sqrt(doa_estimated/trials);
    RMSEs(coordinate,:) = doa_estimated;
end
mean(RMSEs(:,1)) 
max(RMSEs(:,1)) 
min(RMSEs(:,1))
