clear
tic

%% File info

file_name = "third_ventricle_2.tif";

%% ROI info

first_point = [
    270 18
    381 242
    464 529
    529 866
    ];

second_point = [
    385 209
    461 492
    520 814
    537 1040
    ];

%% Parameters
roi_interval = 30; % pixel
expansion = 5; % +-pixel
figure_limit = 100; % Hz
mean_calculation_min = 10; % Hz
mean_calculation_max = 40; % Hz
length_of_video = 1; % sec

% FFT parameters
Fs = 1000;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = length_of_video * Fs;             % Length of signal
t = (0:L-1)*T;        % Time vector

subPlotNumber = 4;

%% ROI info confirm
if sum(size(first_point) == size(second_point)) ~= 2
    disp("ROI point 수 틀림");
    disp("중지됨");
    return
end

%% file read
data_name = split(file_name,".");
data_name = data_name(1);

info = imfinfo(file_name);
image_number_of_page = length(info);
image_bit = info(1).BitsPerSample;
image_max_signal = power(2,image_bit)-1;

if image_number_of_page ~= L
    disp("Frame 수 틀림");
    disp(image_number_of_page + " (이 영상의 frame) ~= " + L + " (설정된 frame)");
    disp("중지됨");
    return
end
image_temp = imread(file_name, 1);
image = zeros(size(image_temp,2), size(image_temp,1), image_number_of_page);

temp = 0;
for k = 1 : image_number_of_page
    if round((k*100)/image_number_of_page) > temp+1
        temp = temp+1;
        file_name + " loading " + round((k*100)/image_number_of_page) + " %"
    end
    image(:,:,k) = imread(file_name, k)';
end

%% make ROI list

[M,I] = sort(first_point(:,1));
first_point = first_point(I,:);
second_point = second_point(I,:);

roi_list = [];

for i = 1 : size(first_point,1)
    x = [first_point(i,1) second_point(i,1)];
    y = [first_point(i,2) second_point(i,2)];
    
    if x(1) > x(2)
        x = flip(x);
        y = flip(y);
    end
    
    slope = (y(2)-y(1)) / (x(2)-x(1));
    x_interval = roi_interval / sqrt(slope^2 + 1);
    y_interval = x_interval * slope;
    roi_x = x(1) : x_interval : x(2);
    roi_y = y(1) : y_interval : y(2);
    
    roi_list = [roi_list; round([roi_x' roi_y'])];
end

%% alnaysis & plot

fftOfAllRoi = [];

Fig = figure('Position', [0 0 1600 800]);

subplot(subPlotNumber,1,1);
hold on;

temp = 0;
for roiOrder = 1:size(roi_list,1)
    if round((roiOrder*100)/size(roi_list,1)) > temp+1
        temp = temp+1;
        file_name + " analysis " + round((roiOrder*100)/size(roi_list,1)) + " %"
    end
    roi = roi_list(roiOrder, :);
    fftOfAllExpansion = [];
    
    for x_cor = (roi(1)-expansion):(roi(1)+expansion)
        for y_cor = (roi(2)-expansion):(roi(2)+expansion)
            data = squeeze(image(x_cor, y_cor, :));
            
            Y = fft(data);
            P2 = abs(Y/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            f = Fs*(0:(L/2))/L;
            f = f(2:end);
            
            fftOfAllExpansion = [fftOfAllExpansion P1(2:end)];
        end
    end
    
    fftOfAllRoi = [fftOfAllRoi mean(fftOfAllExpansion, 2) / max(mean(fftOfAllExpansion, 2))];
    
    plot(f,mean(fftOfAllRoi,2),'color',[0 0 0 5/size(roi_list,1)])
    plot(f,smooth(mean(fftOfAllRoi,2)),'color',[1 0 0 5/size(roi_list,1)])
    xlabel('f (Hz)');
    ylabel('|P1(f)|');
    xlim([1 figure_limit]);
    
end

subplot(subPlotNumber,1,2);
hold on;
plot(f,mean(fftOfAllRoi,2), 'color',[0 0 0]);
plot(f,smooth(mean(fftOfAllRoi,2)), 'color',[1 0 0]);
xlabel('f (Hz)');
ylabel('|P1(f)|');
xlim([1 figure_limit]);

data_heatmap = fftOfAllRoi(1:figure_limit,:)'; % fps 바뀔시 수정 필요

subplot(subPlotNumber,1,3);
h = heatmap(data_heatmap);
h.Colormap = flipud(gray);
h.GridVisible = 'off';
cdx = h.XDisplayLabels;
cdy = h.YDisplayLabels;
h.XDisplayLabels = repmat(' ',size(cdx,1), size(cdx,2));
h.YDisplayLabels = repmat(' ',size(cdy,1), size(cdy,2));

subplot(subPlotNumber,1,4);
[M, I] = max(fftOfAllRoi);
[maxSorted, ISorted] = sort(I);
h = heatmap(data_heatmap(ISorted,:));
h.Colormap = flipud(gray);
h.GridVisible = 'off';
cdx = h.XDisplayLabels;
cdy = h.YDisplayLabels;
h.XDisplayLabels = repmat(' ',size(cdx,1), size(cdx,2));
h.YDisplayLabels = repmat(' ',size(cdy,1), size(cdy,2));

sgtitle({file_name
    size(roi_list,1) + " ROIs"
    "Mean frequency in"
    "whole range = " + mean(I) + " Hz"
    "between " + mean_calculation_min + " Hz ~ " + mean_calculation_max + " Hz = " + mean(I(I>mean_calculation_min & I<mean_calculation_max)) + " Hz"}, 'interpreter', 'none');

time_label = "_" + convertCharsToStrings(datestr(now,'yyyy-mm-dd_HH-MM-SS'));

warning('off')
mkdir('pdf_file')
warning('on')
disp("Saving pdf...")
print(Fig, pwd + "\pdf_file\" + data_name + time_label + "_result.pdf", '-dpdf', '-fillpage', '-painters');
disp("pdf saved.")

warning('off')
mkdir('png_file')
warning('on')
disp("Saving png...")
print(Fig, pwd + "\png_file\" + data_name + time_label + "_result.png", '-dpng');
disp("png saved.")
close all;

Fig = figure;
imshow(image(:,:,1)'/image_max_signal)
hold on
%line(x, y,'Color','red','LineWidth',1)
for i = 1 : size(roi_list,1)
    %rectangle('position', [roi_list(i,1)-(expansion+1) roi_list(i,2)-(expansion+1) (2*expansion+2) (2*expansion+2)]);
    line([roi_list(i,1)-(expansion+1) roi_list(i,1)-(expansion+1)], [roi_list(i,2)-(expansion+1) roi_list(i,2)+(expansion+1)],'Color','red','LineWidth',1)
    line([roi_list(i,1)-(expansion+1) roi_list(i,1)+(expansion+1)], [roi_list(i,2)-(expansion+1) roi_list(i,2)-(expansion+1)],'Color','red','LineWidth',1)
    line([roi_list(i,1)+(expansion+1) roi_list(i,1)-(expansion+1)], [roi_list(i,2)+(expansion+1) roi_list(i,2)+(expansion+1)],'Color','red','LineWidth',1)
    line([roi_list(i,1)+(expansion+1) roi_list(i,1)+(expansion+1)], [roi_list(i,2)+(expansion+1) roi_list(i,2)-(expansion+1)],'Color','red','LineWidth',1)
end

warning('off')
mkdir('roi_file')
warning('on')
disp("Saving roi...")
print(Fig, pwd + "\roi_file\" + data_name + time_label + "_roi.png", '-dpng');
disp("roi saved.")
close all;

warning('off')
mkdir('mat_file')
disp("Saving mat...")
save(pwd + "\mat_file\" + data_name + time_label + ".mat");
disp("mat saved.")
warning('on')

disp("Analysis of " + data_name + time_label + ".tif is finished.")
toc
disp("whole range = " + mean(I) + " Hz");
disp("between " + mean_calculation_min + " Hz ~ " + mean_calculation_max + " Hz = " + mean(I(I>mean_calculation_min & I<mean_calculation_max)) + " Hz");
