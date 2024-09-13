#!/usr/bin/env python

import sys;
import argparse as ap;
import subprocess as sp;
import concurrent.futures
import threading
import multiprocessing

VERSION='0.0.1'

# functions
def CpG_pos(s):
	'''
	find all CpGs in a sequence
	'''
	last_found=-1;
	s=s.upper();
	while True:
		last_found=s.find('CG', last_found+1);
		if last_found == -1:
			break; # no more occurences
		yield last_found;


def get_seq(chrom, start, end):
	'''
	Get sequence for a given region
	'''
	seq=sp.run(['twoBitToFa', args.twobitFile, 'stdout',
		"-seq="+chrom,"-start="+start,"-end="+end], stdout=sp.PIPE,
		universal_newlines=True);
	if seq.stdout == '': return(None);
	seq=seq.stdout.split("\n")[1:];
	seq=''.join(seq);
	return(seq);

# Function to process a chunk of lines
def process_lines(chunk):
    sep, countOnly = chunk[-2], chunk[-1]
    chunk = chunk[:-2]
    results = []
    for l in chunk:
        l = l.strip().split(sep)
        (chrom, start, end, name) = l[:4]
        seq = get_seq(chrom, start, end)
        if seq is None:
            print(f"Cannot get sequence for [{l}]", file=sys.stderr)
            continue
        # get CpG positions now
        poses = list(CpG_pos(seq))
        if countOnly:
            # Count only
            results.append(sep.join(map(str, [chrom, start, end, name, len(poses)])))
            continue
        
        if len(poses) < 1:  # no CpGs
            results.append(sep.join(map(str, [chrom, start, end, name + '.CpG.NA'])))
        else:
            start = int(start)
            for i, p in enumerate(poses, start=1):
                results.append(sep.join(map(str, [chrom, start + p, start + p + 2, name + ".CpG." + str(i)])))
    return results

# producer function
def write_chunk_to_queue(region_file,sep, chunk_size, countOnly, consumer_count):
    with open(region_file, 'r') as r:
        counter = 0
        chunk = []
        for l in r:
            counter += 1
            chunk.append(l)

            if counter % chunk_size == 0:
                print(f"# Start processing line {counter}", file=sys.stderr)
                # Submit chunk for processing
                chunk.extend([sep, countOnly])  # Pass sep and countOnly with chunk
                #futures.append(executor.submit(limited_process_lines, chunk))
                q.put(chunk)
                chunk = []

        # last chunk
        if chunk:
            chunk.extend([sep, countOnly])
            q.put(chunk)
            #futures.append(executor.submit(limited_process_lines, chunk))

    # Put None as sentinel values for each consumer to stop processing
    for _ in range(consumer_count):
        q.put(None)

# consumer function
def process_chunk_from_queue(consumer_id):
    while True:
        chunk = q.get()
        if chunk is None: # Check for the sentinel to stop consuming
            break
        # process the chunk now
        results = process_lines(chunk)
        with write_lock:
            for line in results:
                print(line)
    print(f"Consumer {consumer_id} finished", file=sys.stderr)

# Main parallel processing function using concurrent.futures.ProcessPoolExecutor
def parallel_process(region_file, sep, countOnly, num_workers=4,
        chunk_size=1000):

    # Producer process
    producer_process = multiprocessing.Process(
            target=write_chunk_to_queue,
            args=(region_file,sep,chunk_size,countOnly,num_workers)
            )

    producer_process.start()

    # create consumers using a process pool
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        consumer_futures = [
                executor.submit(process_chunk_from_queue, i) 
                for i in range(num_workers)]    
        
        # wait for all consumers to finish
        for future in concurrent.futures.as_completed(consumer_futures):
            future.result()
            #print("Collected one consumer", file=sys.stderr)
            #consumer_futures.remove(future)  # Remove completed future to free memory

    # Wait for the producer to finish
    producer_process.join()


if __name__ == "__main__":

    desc='''
    This program finds all the CpG positions in given genomic regions.
    The output is in bed format.
    ''';
    
    authorInfo='''
    Author: Zhenguo Zhang
    Email: zhangz.sci@gmail.com
    ''';
    
    # set up arguments
    optParser=ap.ArgumentParser(
    		description=desc,
    		formatter_class=ap.RawTextHelpFormatter,
    		epilog=authorInfo);
    
    ## positional arguments
    optParser.add_argument("regionFile",
    		help="file containing genomic region in bed format, having at least 4 fields");
    
    optParser.add_argument("twobitFile",
    		help="the .2bit file from which sequences can be extracted");
    
    optParser.add_argument("-c", "--count-only",
    		help="if provided, only the number of CpGs in each region is returned",
    		action='store_true',
    		dest='countOnly')
    
    optParser.add_argument("--cpus", default=1, type=int,
    		help="the number of cpus to run in parallel [%(default)d]",
    		action='store',
    		dest='cpus')
    
    optParser.add_argument("--chunk-size", default=1000, type=int,
    		help="the number of regions to process in each batch. A higher value means more memory to use [%(default)d]",
    		action='store',
    		dest='chunkSize')
    
    optParser.add_argument("-v", "--version",
        action='version', version=VERSION,
        help="show version number and exit"
    )

    args=optParser.parse_args();
    sep="\t";

    #max_queue_size = args.cpus * 2  # Limit the number of active submissions
    # use a queue to control how many tasks can be in submission
    q = multiprocessing.Queue(maxsize=args.cpus*2)
    
    # also use a write lock to control writing to output
    write_lock = multiprocessing.Lock()


    # Usage
    parallel_process(args.regionFile, sep, args.countOnly,
            num_workers=args.cpus,
            chunk_size=args.chunkSize)

    print("Job is done", file=sys.stderr);

    sys.exit(0);

