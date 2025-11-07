fun = @system01; 
tamanio_paso = 0.1;
tiempo_simulacion = 50;
iteraciones = tiempo_simulacion / tamanio_paso;
cao = 0.082345311795785; % mol / L
To = 441.81; % K
qc = 100; % Coolant flow rate: L/min
sigma = 5;

y1 = cao;
y2 = To;
z1 = cao;
y1nn = cao;
y = zeros(1, iteraciones + 1);
z = zeros(1, iteraciones + 1);
u = zeros(1, iteraciones + 1);
ynn = zeros(1, iteraciones + 1);
datos = zeros(iteraciones,3);
u(1) = qc;
y(1) = y1;
z(1) = z1;
ynn(1) = y1nn;
tiempo = 0:tamanio_paso:tiempo_simulacion;
contador = 0;

for iter=1:iteraciones
    contador = contador + 1;
    if iter == 1
        qcr = qc + sigma * randn;
    elseif (contador > 10) && (abs(mean(y(iter-9:iter-3)) - mean(y(iter-2:iter-1))) < 0.0001)
        qcr = qc + sigma * randn;
        contador = 0;
    end
    [y1, y2] = runge_kutta(fun, y1, y2, qcr, tamanio_paso);
    z1 = -0.04253 + 0.0005342 * qcr + 0.8667 * z1;
    y1nn = myNeuralNetworkFunction([qcr; y1nn]);
    y(iter+1) = y1;
    u(iter+1) = qcr;
    z(iter+1) = z1;
    ynn(iter+1) = y1nn;
    datos(iter,:) = [qcr, y(iter), y(iter+1)];
end

figure
plot(tiempo, u)
figure
plot(tiempo, y, tiempo, ynn, tiempo, z)
hold on

t_m = tiempo(1:2);
y_m = y(1:2);
curvaturas = [];
for i=3:length(y)-1
    criterio_a = rand<0.33333333;%acepta_y(y(i-2), y(i-1), y(i), 0.01, 0.99, 0.5);
    criterio_b = (y(i)-y(i-1)) * (y(i) - y(i+1)) > 0;
    criterio_c =0;% angulo_3p([tiempo(i-1), y(i-1)], [tiempo(i), y(i)], [tiempo(i+1), y(i+1)]) < 179.9;
%     curvaturas = [curvaturas, angulo_3p([tiempo(i-1), y(i-1)], [tiempo(i), y(i)], [tiempo(i+1), y(i+1)])];
    if criterio_a || criterio_b || criterio_c
        t_m = [t_m, tiempo(i)];
        y_m = [y_m, y(i)];        
    end
end
plot(t_m, y_m, '.k')
% plot(tiempo(2:end-1), curvaturas/2000, '.r')
