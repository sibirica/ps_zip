function [EKG, basket] = readBasket(filename)
A = dlmread(filename);
A = A';
EKG = A(1:12,:);
basket = A(13:76,:);
basket = reshape(QRSfilt(basket, EKG), 8, 8, []); % number, letter, t
end