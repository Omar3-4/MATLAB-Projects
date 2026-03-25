clear,clc
reset(gpuDevice(1));
tic;
T_Start = tic;
Final_P = [];
K_Limit = uint32(2147483647); 
Chunk = uint32(100000000); 

for s = uint32(1) : Chunk : K_Limit
    e = min(s + Chunk - 1, K_Limit);
    k = gpuArray(uint32(s : e));
    
    k_plus_1 = k + uint32(1);
    mask = bitand(k_plus_1, k) == uint32(0);
    
    if any(mask)
        candidates_k = gather(k(mask));
        for i = 1:length(candidates_k)
            K_val = uint64(candidates_k(i));
            
            if isprime(double(K_val))
                % The Physics-Engine Move: Bit-Shift instead of Multiplication
                p = round(log2(double(K_val) + 1));
                % N = K * 2^(p-1) -> logic shift left
                N_val = bitshift(K_val, p - 1);
                
                Final_P = [Final_P, N_val];
                fprintf('\n[SUCCESS: %lu] | K = %u | Time: %.2f s', N_val, uint32(K_val), toc(T_Start));
                
                fid = fopen('P8_BitPerfect.txt', 'a');
                fprintf(fid, '%lu | %s\n', N_val, datestr(now));
                fclose(fid);
            end
        end
    end
    clear k k_plus_1 mask;
end

fprintf('\n\nScan Finished. Total Time: %.2f mins\n', toc(T_Start)/60);
disp('The 8th Perfect Number (100% Precise):');
fprintf('%lu\n', unique(Final_P));