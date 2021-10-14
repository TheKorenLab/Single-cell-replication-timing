function rgb_color = convert_hex(hex_color)

    assert(strcmp(hex_color(1), '#') == true, 'Invalid input: color string must begin with #.')
    assert(length(hex_color) == 7, 'Invalid input: wrong length.')
    
    hex_color = hex_color(2:end);
    rgb_color = sscanf(hex_color,'%2x%2x%2x',[1 3]) ./ 255;
    
    assert(all(rgb_color >= 0) & all(rgb_color <= 1), 'Invalid color output.')
    
end
