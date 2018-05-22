function [tensorinfo_vec,tensorinfo_w,tensorinfo_length] = tensorinfo(rule,deltagrid,ind)

% helper function for the Smolyak algorithm
%
% Bettina Schieche, 2011

tensorinfo_vec = {};
tensorinfo_w = {};
tensorinfo_length = [];
for i = 1:length(rule)-1
    I = rule{end}==i;
    tensorinfo_vec = [tensorinfo_vec,{deltagrid.(char(rule{i})){ind(I),1}}];
    tensorinfo_w = [tensorinfo_w,{deltagrid.(char(rule{i})){ind(I),2}}];
    tensorinfo_length = [tensorinfo_length,[deltagrid.(char(rule{i})){ind(I),3}]];
end
