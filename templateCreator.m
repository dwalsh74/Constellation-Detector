
clear all;
close all;
clc;

%% Template 

imColor = imread('Virgo_Googled.png');
imGray = rgb2gray(imColor);
binaryTemplate = imGray > 30;


figure;imshow(imColor);title('Color Image');
figure;imshow(binaryTemplate);title('Binary Image before opening/closing');

%Fill holes
struct = strel('disk',2);
binaryTemplate = imclose(binaryTemplate,struct);

%Remove small stars
struct = strel('square',2);
binaryTemplate = imopen(binaryTemplate,struct);

figure;imshow(binaryTemplate);title('Binary Image after opening/closing');

