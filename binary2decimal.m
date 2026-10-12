function y = binary2decimal(code)
y = code*(2.^(length(code)-1:-1:0)');
end
