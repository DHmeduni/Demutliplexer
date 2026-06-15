import argparse
import glob
import regex as re  # Use the regex library instead of re
import itertools
import os
import multiprocessing
from Bio.Seq import Seq

def load_barcodes(barcode_file, split=None):
    barcodes = {}
    with open(barcode_file, "r") as f:
        for line_number, line in enumerate(f, start=1):
            if line_number < 5:
                continue
            if line.strip():
                fields = line.strip().split("\t")
                if len(fields) >= 3 and split:
                    sample_name, barcode1, barcode2 = fields[:3]
                    barcodes[(barcode1, barcode2)] = sample_name
                elif len(fields) >= 2 and not split:
                    sample_name, barcode = fields[:2]
                    barcodes[(barcode, None)] = sample_name
                # Additional columns are ignored but allowed
    return barcodes

def find_barcode_pair(sequence, barcodes, max_mismatches, search_window, require_both):
    for (barcode1, barcode2), sample_name in barcodes.items():
        rc_barcode1 = str(Seq(barcode1).reverse_complement())
        rc_barcode2 = str(Seq(barcode2).reverse_complement()) if barcode2 else None

        # Patterns allowing mismatches
        pattern1 = fr'({re.escape(barcode1)}){{e<={max_mismatches}}}'
        rc_pattern1 = fr'({re.escape(rc_barcode1)}){{e<={max_mismatches}}}'
        #pattern2 = fr'({re.escape(barcode2)}){{e<={max_mismatches}}}' if barcode2 else None
        #rc_pattern2 = fr'({re.escape(rc_barcode2)}){{e<={max_mismatches}}}' if rc_barcode2 else None

        # Search at the beginning
        match1 = re.search(pattern1, sequence[:search_window + len(barcode1)], re.BESTMATCH)
        end_match1 = re.search(rc_pattern1, sequence[-(search_window + len(rc_barcode1)):], re.BESTMATCH)

        # If `require_both` is True, check if both start and end matches exist
        if require_both and not barcode2:
            if match1 and end_match1:
                    return sample_name
        if not barcode2 and not require_both:
            if match1 or end_match1:
                    return sample_name

        # If one of the barcodes is found at the beginning
        if match1 or end_match1:
            # If barcode2 is present and found at the beginning, check for the corresponding barcode at the end

            expected_end_barcode1 = None
            expected_end_barcode2 = None


            if match1:
                expected_end_barcode1 = rc_barcode2 if barcode2 else None
            else:
                expected_end_barcode2 = barcode2 if barcode2 else None

            if expected_end_barcode1:
                # Check for the second barcode at the end
                end_match1 = re.search(fr'({re.escape(expected_end_barcode1)}){{e<={max_mismatches}}}',
                                      sequence[-(search_window + len(expected_end_barcode1)):], re.BESTMATCH)
                if end_match1:
                    return sample_name
            if expected_end_barcode2:
                # Check for the second barcode at the end
                match2 = re.search(fr'({re.escape(expected_end_barcode2)}){{e<={max_mismatches}}}',
                                      sequence[:search_window + len(expected_end_barcode2)], re.BESTMATCH)
                if match2:
                    return sample_name


    return None

def process_fastq_chunk(chunk, barcodes, max_mismatches, search_window, require_both):
    sorted_reads = {sample_name: [] for sample_name in barcodes.values()}
    processed_count = 0

    for i in range(0, len(chunk), 4):
        header, sequence, plus, quality = [line.strip() + "\n" for line in chunk[i:i+4]]
        sample_name = find_barcode_pair(sequence.strip(), barcodes, max_mismatches, search_window, require_both)

        if sample_name:
            sorted_reads[sample_name].extend([header, sequence, plus, quality])
        processed_count += 1

    return sorted_reads, processed_count

def read_chunks(input_file, chunk_size, queue):
    with open(input_file, 'r') as f:
        while True:
            chunk = list(itertools.islice(f, chunk_size))
            if not chunk:
                break
            queue.put(chunk)
    queue.put(None)  # Sentinel value to signal end of file

def worker(queue, barcodes, max_mismatches, search_window, require_both, output_dir, read_counter, lock, total_reads, progress_step, worker_id):
    temp_files = {sample_name: open(os.path.join(output_dir, f'temp_{sample_name}_{worker_id}.fastq'), 'w') for sample_name in barcodes.values()}
    
    while True:
        chunk = queue.get()
        if chunk is None:
            break

        sorted_reads, chunk_reads = process_fastq_chunk(chunk, barcodes, max_mismatches, search_window, require_both)

        with lock:
            for sample_name, reads in sorted_reads.items():
                temp_files[sample_name].writelines(reads)
            read_counter.value += chunk_reads
            
            if read_counter.value >= progress_step.value:
                print(f"Processed {read_counter.value} reads...")
                progress_step.value += total_reads.value // 10

    for f in temp_files.values():
        f.close()

def estimate_total_reads(input_file):
    with open(input_file, 'r') as f:
        return sum(1 for _ in f) // 4  # Each read spans 4 lines

def merge_temp_files(output_dir, barcodes):
    for sample_name in barcodes.values():
        output_file_path = os.path.join(output_dir, f'{sample_name}.fastq')
        with open(output_file_path, 'w') as outfile:
            temp_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.startswith(f'temp_{sample_name}_')]
            for temp_file in sorted(temp_files):
                with open(temp_file, 'r') as infile:
                    outfile.write(infile.read())
                os.remove(temp_file)

def main():
    parser = argparse.ArgumentParser(description='Demultiplex FASTQ data based on barcodes with mismatch tolerance.')
    parser.add_argument('-i', '--input', required=True, help='Input FASTQ file')
    parser.add_argument('-b', '--barcodes', required=True, help='Barcode file')
    parser.add_argument('-o', '--output', required=True, help='Output folder')
    parser.add_argument('-m', '--mismatches', type=int, default=2, help='Max mismatches')
    parser.add_argument('-w', '--window', type=int, default=60, help='Search window size')
    parser.add_argument('-r', '--require_both', action='store_true', help='Require barcode at both ends')
    parser.add_argument('-p', '--processes', type=int, default=4, help='Number of processes')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    parser.add_argument('-s', '--split', action='store_true', help='Split barcodes ')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    barcodes = load_barcodes(args.barcodes, args.split)

    if os.path.isfile(args.input):
        pass  # If it's a file, do nothing, continue with args.input as is
    else:
        # If it's not a file, search for a .fastq file in the directory
        fastq_files = glob.glob(os.path.join(args.input, "*.fastq"))
    
        if not fastq_files:
            exit(1)  # Or handle the error accordingly if no .fastq file is found
        else:
            args.input = fastq_files[0]  # Use the first found .fastq file
    
    queue = multiprocessing.Queue(maxsize=args.processes * 2)
    lock = multiprocessing.Lock()
    read_counter = multiprocessing.Value('i', 0)
    total_reads = multiprocessing.Value('i', estimate_total_reads(args.input))
    progress_step = multiprocessing.Value('i', total_reads.value // 10)

    reader_process = multiprocessing.Process(target=read_chunks, args=(args.input, 4000, queue))
    reader_process.start()

    workers = []
    for worker_id in range(args.processes):
        worker_process = multiprocessing.Process(target=worker, args=(queue, barcodes, args.mismatches, args.window, args.require_both, args.output, read_counter, lock, total_reads, progress_step, worker_id))
        worker_process.start()
        workers.append(worker_process)

    reader_process.join()

    for _ in range(args.processes):
        queue.put(None)

    for worker_process in workers:
        worker_process.join()

    merge_temp_files(args.output, barcodes)

    if args.verbose:
        print(f"Processing complete: {read_counter.value} total reads processed.")

if __name__ == '__main__':
    main()
