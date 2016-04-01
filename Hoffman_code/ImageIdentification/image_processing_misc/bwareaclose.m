function filtered_image = bwareaclose(bw_image,max_size)

filtered_image = bwpropclose(bw_image,'Area',max_size);

end