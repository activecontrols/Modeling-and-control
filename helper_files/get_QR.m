function [Q_matrix, R_matrix] = get_QR(max_x, max_u)
Q_matrix = zeros(length(max_x), length(max_x));
R_matrix = zeros(length(max_u), length(max_u));
for i = 1:length(max_x)
    Q_matrix(i, i) = 1/(max_x(i)^2);
end
for j = 1:length(max_u)
    R_matrix(j, j) = 1/(max_u(j)^2);
end
end