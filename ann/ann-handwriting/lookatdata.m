%% Read training set


%
% http://neuralnetworksanddeeplearning.com/chap1.html
%
% http://yann.lecun.com/eadb/mnist/
%
% Class function based around:
% https://machinelearningmastery.com/implement-backpropagation-algorithm-scratch-python/
%

% Training set images
filename = 'train-images-idx3-ubyte';
fid = fopen(filename,'r');
train = [];
train.magic = fread(fid,1,'uint32','ieee-be');
train.noimages = fread(fid,1,'uint32','ieee-be');
train.nrows = fread(fid,1,'uint32','ieee-be');
train.ncols = fread(fid,1,'uint32','ieee-be');
train.image = zeros(train.nrows,train.ncols,train.noimages);
for k = (1:train.noimages)
    train.image(:,:,k) = reshape(fread(fid,train.nrows*train.ncols,'ubit8','ieee-be'),train.nrows,train.ncols);
end
fclose(fid);
disp('Training data loaded.');

% Training set labels
filename = 'train-labels-idx1-ubyte';
fid = fopen(filename,'r');
train.magic_labels = fread(fid,1,'uint32','ieee-be');
train.noimages_labels = fread(fid,1,'uint32','ieee-be');
train.labels = fread(fid,train.noimages_labels,'ubit8','ieee-be');
fclose(fid);
disp('Training labels loaded.');

% Test set images
filename = 't10k-images-idx3-ubyte';
fid = fopen(filename,'r');
test = [];
test.magic = fread(fid,1,'uint32','ieee-be');
test.noimages = fread(fid,1,'uint32','ieee-be');
test.nrows = fread(fid,1,'uint32','ieee-be');
test.ncols = fread(fid,1,'uint32','ieee-be');
test.image = zeros(test.nrows,test.ncols,test.noimages);
for k = (1:test.noimages)
    test.image(:,:,k) = reshape(fread(fid,test.nrows*test.ncols,'ubit8','ieee-be'),test.nrows,test.ncols);
end
fclose(fid);
disp('Test data loaded.');

% Test set labels
filename = 't10k-labels-idx1-ubyte';
fid = fopen(filename,'r');
test.magic_labels = fread(fid,1,'uint32','ieee-be');
test.noimages_labels = fread(fid,1,'uint32','ieee-be');
test.labels = fread(fid,test.noimages_labels,'ubit8','ieee-be');
fclose(fid);
disp('Test labels loaded.');

% Cleanup
clear('k','fid','filename','ans');

%% Neural network implementation

addpath ~/Dropbox/Matlab/library/ANN/m/

% Nominate training set
dat = train;

% Specify smoothing kernel on INPUT image
kernel_smoothing = 1;           % 1 = none

% Specify network architecture
input_nodes = dat.nrows * dat.ncols;
hidden_nodes = [ 30 ];
output_nodes = 10;

% Convolve input images
kernel = ones(kernel_smoothing);
kernel = kernel/sum(kernel(:));
if ( kernel_smoothing > 1 )
    for k = (1:size(dat.image,3))
        dat.image(:,:,k) = conv2( dat.image(:,:,k), kernel, 'same' );
    end
end

% Restructure training set
nImages = 1000;%size(dat.image,3);
index = @(i,j,nSize) (j-1)*nSize(1)+i;
input_data = reshape(dat.image(:,:,1:nImages),input_nodes,nImages);
input_data = input_data/max(input_data(:));
output_labels = zeros(output_nodes,nImages); output_labels(index(dat.labels(1:nImages)+1,(1:nImages)',[output_nodes nImages])) = 1;

% Initialise network
network = ANN();
network.SetArchitecture( [ input_nodes hidden_nodes output_nodes ] );
%network.SetTransferFunction( 'ReLU' );
network.SetData( input_data, output_labels );
network.Initialise();

% Determine untrained accuracy
predicted = network.Predict();
[num,~] = find( predicted == repmat(max(predicted,[],1),size(predicted,1),1) );
if ( length(num) == size(predicted,2) )
    disp(['Untrained accuracy (training set) ' num2str(mean((dat.labels(1:nImages)+1)==num))]);
else
    disp('Untrained accuracy unknown - too many viable candidates for each response.');
end

% Batch train on training data
%network.Train();                           % Shows exceptional performance but very slow and the current SGD implementation suffers from a halting problem
%network.Train('minibatch',100);            % Train on all trials using SGD split into n=100 minibatches
%network.Train('onehit',1);                 % One-hit learning; learn after every trial. Slow and only runs through dataset once.
network.Train('random',1000);              % Randomly select n=1000 trials each training iteration; good but does not converge as well as full training

% Determine trained accuracy
predicted = network.Predict();
[labels,~] = find( output_labels == repmat(max(output_labels,[],1),size(output_labels,1),1) );
[num,~] = find( predicted == repmat(max(predicted,[],1),size(predicted,1),1) );
if ( length(num) == size(predicted,2) )
    disp(['Trained accuracy (training set) ' num2str(mean((dat.labels(1:nImages)+1)==num))]);
    correct = sum(num==labels);
    fprintf(' Number correct %g/%g (%g%%)\n\n',correct,length(num),100*correct/length(num));
end

% Check against test set
dat = test;
if ( kernel_smoothing > 1 )
    for k = (1:size(dat.image,3))
        dat.image(:,:,k) = conv2( dat.image(:,:,k), kernel, 'same' );
    end
end
nImages = size(dat.image,3);
input_data = reshape(dat.image(:,:,1:nImages),input_nodes,nImages);
input_data = input_data/max(input_data(:));
output_labels = zeros(output_nodes,nImages); output_labels(index(dat.labels(1:nImages)+1,(1:nImages)',[output_nodes nImages])) = 1;
network.SetData( input_data, output_labels );       % Set to Test data
predicted = network.Predict();
[labels,~] = find( output_labels == repmat(max(output_labels,[],1),size(output_labels,1),1) );
[num,~] = find( predicted == repmat(max(predicted,[],1),size(predicted,1),1) );
if ( length(num) == size(predicted,2) )
    disp(['Trained accuracy (training set) ' num2str(mean((dat.labels(1:nImages)+1)==num))]);
    correct = sum(num==labels);
    fprintf(' Number correct %g/%g (%g%%)\n\n',correct,length(num),100*correct/length(num));
    % Per digit accuracy and confidence
    accuracy = zeros(1,10);
    fpr = accuracy;     % False positive rate
    confDigits = accuracy;
    confAll = max(predicted,[],1)./sum(predicted,1);
    for k = (1:10)
        kk = labels==k;
        accuracy(k) = sum(num(kk)==labels(kk))/sum(kk);
        fpr(k) = sum(num(labels~=k)==k)/sum(labels~=k);
        confDigits(k) = mean( confAll(kk) );       % Average confidence for that digit
    end
    % Confidence, split by whether classification was correct or incorrect
    if ( false )
        addpath ~/Dropbox/Matlab/thirdparty/violinplot/
        figure;
        violinplot( confAll, num==labels );
        xlabel('Classification accuracy');
        ylabel('Confidence in classification');
    end
end

% Prune and test again
if ( 1 )
    threshold = 0.05; make_sparse = true;
    network.Prune( threshold, make_sparse );
    disp(['Pruned ' num2str(100*mean(full(network.weights{2}(:))==0)) '% of hidden weights']);
    predicted = network.Predict();
    [labels,~] = find( output_labels == repmat(max(output_labels,[],1),size(output_labels,1),1) );
    [num,~] = find( predicted == repmat(max(predicted,[],1),size(predicted,1),1) );
    if ( length(num) == size(predicted,2) )
        disp(['Trained accuracy (pruned trained set) ' num2str(mean((dat.labels(1:nImages)+1)==num))]);
        correct = sum(num==labels);
        fprintf(' Number correct %g/%g (%g%%)\n\n',correct,length(num),100*correct/length(num));
    end
end

%% Examine output weights

cmap = [ [ (0:0.1:1)' 1*ones(11,1) (0:0.1:1)' ];
repmat([1 1 1],10,1);
[ (1:-0.1:0)' (1:-0.1:0)' 1*ones(11,1) ] ];

% Ouput layer
figure;
W = network.weights{end};
imagesc( 1:size(network.weights{end},2), (1:size(network.weights{end},1))-1, W );
colormap(cmap);
set(gca,'clim',max(abs(get(gca,'clim')))*[-1 1]);
set(gca,'xtick',(1:size(network.weights{end},2)),'ytick',(0:size(network.weights{end},1)));
axis image;
ylabel('Output value');
xlabel('weight from hidden layers');
title('Output mapping from hidden layer');

%% Examine hidden node weights

% Hidden layers
figure; layer = 2;
cols = ceil(sqrt(size(network.weights{layer},1)));
rows=ceil(size(network.weights{layer},1)/cols);
for neuron = 1:size(network.weights{layer},1)
    subplot(cols,rows,neuron);
    imagesc( reshape( network.weights{layer}(neuron,:).', [ dat.ncols dat.nrows ] ) )
    colormap(cmap);
    set(gca,'clim',max(abs(get(gca,'clim')))*[-1 1]);
    colorbar;
    axis off; axis image;
    title(neuron);
end

% Correlate weight matrices for hidden nodes (find redundancies)
R = zeros(size(network.weights{2},1));
for i = (1:size(network.weights{2},1))
    for j = (1:size(network.weights{2},1))
        R(i,j) = diag(corrcoef( full(network.weights{2}(i,:)), full(network.weights{2}(j,:)) ),1);
    end
end
figure;
subplot(221);
    imagesc(R);
    title('Hidden layer correlations');
subplot(222);
    dendrogram(linkage(R));
    title('Hidden layer correlations');

%% Examine (output weighted) hidden nodes

kernel = ones(4); kernel = kernel/sum(kernel(:));

figure;
cols=4;
rows=3;
for output_value = (0:9)

    % Weighting of hidden layer (weighted by output)
    node = output_value+1;      % ( 0->9 ) -> index
    ws = network.weights{end}(node,:);
    M = zeros([dat.nrows dat.ncols]);
    for neuron = (1:size(network.weights{end-1},1))
        M = M + ws(neuron)*conv2(reshape( full(network.weights{2}(neuron,:)), [ dat.ncols dat.nrows ] ).',kernel,'same');
    end
    
    subplot(cols,rows,node);
    imagesc(conv2(M,kernel,'same'));
    colormap(cmap); axis image
    set(gca,'clim',max(abs(get(gca,'clim')))*[-1 1]);
    colorbar;
    title(sprintf('Digit %g: Hit-rate = %.1f%%, FPR = %.1f%%, Confidence = %.1f%%', output_value, 100*accuracy(node), 100*fpr(node), 100*confDigits(node) ));
end
