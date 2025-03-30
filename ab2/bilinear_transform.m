% Exemplo 7.3 - Transformação Bilinear de um Filtro Butterworth
% Script MATLAB para validar diagramas de Bode e testar desempenho

% Parâmetros do filtro
N = 6;          % Ordem do filtro
Omega_c = 0.766; % Frequência de corte normalizada
Ts = 1;         % Período de amostragem (por conveniência)

% Gerar filtro Butterworth de tempo contínuo (analógico)
[num_c, den_c] = butter(N, Omega_c, 's');
H_c = tf(num_c, den_c);

% Aplicar transformação bilinear para obter filtro de tempo discreto
H_d = c2d(H_c, Ts, 'tustin');
[num_d, den_d] = tfdata(H_d, 'v');

% Frequências críticas para verificação
w_p = 0.2*pi;   % Frequência de passagem (banda passante)
w_s = 0.3*pi;   % Frequência de corte (banda de rejeição)

% Pré-distorção das frequências para correção do warping
Omega_p = 2*tan(w_p/2);
Omega_s = 2*tan(w_s/2);

% Configuração para plotagem
omega_c = logspace(-2, 1, 1000);  % Frequências para filtro analógico

% Calcular respostas em frequência do filtro analógico
[mag_c, phase_c] = bode(H_c, omega_c);
mag_c = squeeze(mag_c);
phase_c = squeeze(phase_c);

% Calcular resposta em frequência do filtro digital manualmente
omega_d = linspace(0, pi, 1000);
mag_d = zeros(size(omega_d));
phase_d = zeros(size(omega_d));

% Avalia a resposta em frequência diretamente para cada frequência
for i = 1:length(omega_d)
    z = exp(1j * omega_d(i));
    Hz = polyval(num_d, z) / polyval(den_d, z);
    mag_d(i) = abs(Hz);
    phase_d(i) = angle(Hz) * 180/pi;
end

% Valores para verificação nas frequências críticas
[mag_wp, ~] = bode(H_c, Omega_p);
mag_wp_dB = 20*log10(squeeze(mag_wp));
[mag_ws, ~] = bode(H_c, Omega_s);
mag_ws_dB = 20*log10(squeeze(mag_ws));

% Calcular resposta em frequências específicas
z_wp = exp(1j * w_p);
H_wp = polyval(num_d, z_wp) / polyval(den_d, z_wp);
mag_wp_d_dB = 20*log10(abs(H_wp));

z_ws = exp(1j * w_s);
H_ws = polyval(num_d, z_ws) / polyval(den_d, z_ws);
mag_ws_d_dB = 20*log10(abs(H_ws));

% Plotar diagramas de Bode comparativos
figure('Position', [100, 100, 1000, 800]);

% Magnitude do filtro analógico
subplot(2,2,1);
semilogx(omega_c, 20*log10(mag_c));
grid on;
title('Magnitude do Filtro Analógico H(s)');
xlabel('Frequência (rad/s)');
ylabel('Magnitude (dB)');
hold on;
plot(Omega_p, mag_wp_dB, 'ro', 'MarkerSize', 8);
plot(Omega_s, mag_ws_dB, 'rx', 'MarkerSize', 8);
legend('Resposta em magnitude', ['Em \Omega_p: ' num2str(mag_wp_dB,'%.2f') ' dB'], ...
       ['Em \Omega_s: ' num2str(mag_ws_dB,'%.2f') ' dB'], 'Location', 'southwest');

% Fase do filtro analógico
subplot(2,2,3);
semilogx(omega_c, phase_c);
grid on;
title('Fase do Filtro Analógico H(s)');
xlabel('Frequência (rad/s)');
ylabel('Fase (graus)');

% Magnitude do filtro digital
subplot(2,2,2);
plot(omega_d, 20*log10(mag_d));
grid on;
title('Magnitude do Filtro Digital H(z)');
xlabel('Frequência (rad)');
ylabel('Magnitude (dB)');
hold on;
plot(w_p, mag_wp_d_dB, 'ro', 'MarkerSize', 8);
plot(w_s, mag_ws_d_dB, 'rx', 'MarkerSize', 8);
legend('Resposta em magnitude', ['Em \omega_p: ' num2str(mag_wp_d_dB,'%.2f') ' dB'], ...
       ['Em \omega_s: ' num2str(mag_ws_d_dB,'%.2f') ' dB'], 'Location', 'southwest');

% Fase do filtro digital
subplot(2,2,4);
plot(omega_d, phase_d);
grid on;
title('Fase do Filtro Digital H(z)');
xlabel('Frequência (rad)');
ylabel('Fase (graus)');

% Análise adicional: Resposta ao degrau de ambos os filtros
figure('Position', [100, 100, 1000, 400]);

% Resposta ao degrau do filtro analógico
subplot(1,2,1);
step(H_c, 10);
title('Resposta ao Degrau do Filtro Analógico H(s)');
grid on;

% Resposta ao degrau do filtro digital através de simulação direta
subplot(1,2,2);
t = 0:20; % 21 amostras
u = ones(size(t)); % Sinal de entrada degrau
u(1) = 0; % Começa em t=1

% Simular a resposta ao degrau através de filtragem direta
y = filter(num_d, den_d, u);

% Plotar resposta
stem(t, y);
title('Resposta ao Degrau do Filtro Digital H(z)');
grid on;
xlabel('Amostras');
ylabel('Amplitude');

% Exibição dos resultados no console
disp('--- Resultados da Transformação Bilinear ---');
disp(['Ordem do filtro: N = ' num2str(N)]);
disp(['Frequência de corte analógica: Omega_c = ' num2str(Omega_c) ' rad/s']);
disp(' ');
disp('Valores na banda passante e de rejeição:');
disp(['Filtro analógico em Omega_p = ' num2str(Omega_p) ' rad/s: ' num2str(mag_wp_dB,'%.2f') ' dB']);
disp(['Filtro analógico em Omega_s = ' num2str(Omega_s) ' rad/s: ' num2str(mag_ws_dB,'%.2f') ' dB']);
disp(['Filtro digital em omega_p = ' num2str(w_p) ' rad: ' num2str(mag_wp_d_dB,'%.2f') ' dB']);
disp(['Filtro digital em omega_s = ' num2str(w_s) ' rad: ' num2str(mag_ws_d_dB,'%.2f') ' dB']);

% Análise de desempenho
disp(' ');
disp('--- Análise de Polos e Zeros ---');
disp('Filtro Analógico H(s):');
disp('Zeros:');
disp(zero(H_c));
disp('Polos:');
disp(pole(H_c));

disp('Filtro Digital H(z):');
disp('Zeros:');
disp(zero(H_d));
disp('Polos:');
disp(pole(H_d));

% Análise do mapeamento de frequência - ilustração do warping
figure;
w_range = linspace(0, pi, 100);
Omega_range = 2*tan(w_range/2);

plot(w_range/pi, Omega_range, 'b-', 'LineWidth', 2);
hold on;
plot([0 1], [0 pi], 'r--');
grid on;
title('Distorção de Frequência (Warping) na Transformação Bilinear');
xlabel('Frequência Digital Normalizada (\omega/\pi)');
ylabel('Frequência Analógica (rad/s)');
legend('Mapeamento de Frequência', 'Mapeamento Linear', 'Location', 'northwest');

% Verificação adicional: coeficientes da função de transferência
disp(' ');
disp('--- Função de Transferência do Filtro Digital H(z) ---');
disp('Numerador:');
disp(num_d);
disp('Denominador:');
disp(den_d);

% Fórmula da tranformação bilinear para os parâmetros do problema
disp(' ');
disp('Verificação da função de transferência conforme o Exemplo 7.3:');
disp(['|H(e^(j*0.2*π))| ≈ 0.89125 ou -0.36 dB (calculado: ' num2str(mag_wp_d_dB,'%.2f') ' dB)']);
disp(['|H(e^(j*0.3*π))| ≈ 0.17783 ou -15 dB (calculado: ' num2str(mag_ws_d_dB,'%.2f') ' dB)']);
