import csv
import hashlib
import sys

def hash_to_unsigned_longlong(key_str):
    hash_obj = hashlib.md5(key_str.encode('utf-8'))
    hash_bytes = hash_obj.digest()[:8]
    return int.from_bytes(hash_bytes, byteorder='big', signed=False)

def process_mawi_csv(input_file, output_file):
    
    unique_keys = set()
    total_lines = 0
    output_lines = 0
    
    try:
        with open(input_file, 'r', encoding='utf-8') as csvfile, \
             open(output_file, 'w', encoding='utf-8') as outfile:
            
            reader = csv.reader(csvfile)
            
            header = next(reader)

            processed = 0
            
            for row in reader:
                if len(row) < 2:
                    continue
                
                key_str, time_diff_str = row[0].strip(), row[1].strip()
                
                try:
                    time_diff = int(time_diff_str)
                    if time_diff == 0:
                        continue
                except ValueError:
                    continue
                
                hashed_key = hash_to_unsigned_longlong(key_str)
                
                outfile.write(f"{hashed_key},{time_diff}\n")
                
                unique_keys.add(key_str)
                total_lines += 1
                output_lines += 1
                
                processed += 1
                if processed % 1000000 == 0:
                    print(f"processed: {processed:,} lines")
            

    except Exception as e:
        print(f"error: {str(e)}")
        return


    stats_file = "processing_stats.txt"
    with open(stats_file, 'w', encoding='utf-8') as f:
        f.write("MAWI.csv\n")
        f.write("="*40 + "\n")
        f.write(f"unique keys:: {len(unique_keys):,}\n")
        f.write(f"total lines: {output_lines:,}\n")


def main():
    #input_file = "CAIDA.csv"
    #output_file = "CAIDA.txt"
    input_file = "MAWI.csv"
    output_file = "MAWI.txt"


    process_mawi_csv(input_file, output_file)
    

if __name__ == "__main__":
    main()
