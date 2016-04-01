function high_pass_image = apply_high_pass_filter(image,filter_size)

I_filt = fspecial('disk',filter_size);
blurred_image = imfilter(image,I_filt,'same',mean(image(:)));
high_pass_image = image - blurred_image;

end