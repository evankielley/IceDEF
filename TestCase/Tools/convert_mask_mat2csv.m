clear;

root_dir = '/home/evankielley/IceDEF/TestCase/';
input_dir = strcat(root_dir,'Inputs/');
filename = 'mask.csv';

load(strcat(input_dir,'mask.mat'));

csvwrite(strcat(input_dir, filename), msk);

msk_copy = csvread(strcat(input_dir, filename));

if isequal(msk, msk_copy)
      disp('msk Data match');
else
      disp('msk Data mis-match');
end
