function label_out = renumber_label_to_start_at_one(label_in)

label_out = zeros(size(label_in));

label_nums = nonzeros(unique(label_in(:)));

for index = 1:length(label_nums)
    label_out(label_in == label_nums(index)) = index;
end