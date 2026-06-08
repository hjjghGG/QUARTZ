import os

def process_large_file_optimized(input_file, output_file, max_lines=25000000):

    processed_lines = 0
    batch_size = 500000
    
    with open(input_file, 'r', encoding='utf-16') as infile, \
         open(output_file, 'w', encoding='utf-16') as outfile:
        
        prev_timestamp = None
        batch_results = []
        
        for line in infile:
            if processed_lines >= max_lines:
                break
                
            elements = line.strip().split(',')
            if len(elements) <= 3:
                continue
            
            merged = ''.join(elements[:-1])
            try:
                merged_long = int(merged)
            except ValueError:
                merged_long = hash(merged) & 0xFFFFFFFFFFFFFFFF
            
            timestamp_clean = elements[-1].replace('.', '')
            current_timestamp = int(timestamp_clean)
            
            timestamp_diff = current_timestamp - prev_timestamp if prev_timestamp is not None else 0
            prev_timestamp = current_timestamp
            
            batch_results.append(f"{merged_long} {timestamp_diff}\n")
            processed_lines += 1
            
            # 批量写入
            if len(batch_results) >= batch_size:
                outfile.writelines(batch_results)
                batch_results = []
                print(f"processed: {processed_lines} lines")
        
        if batch_results:
            outfile.writelines(batch_results)
    
    print(f"complete with {processed_lines} lines")

# 使用优化版本
process_large_file_optimized("out_2.txt", "output.txt")
