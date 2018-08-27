#!/usr/bin/env python
import sys
import os
import re
import datetime

output_folder = 'assigned'

def write_line(barcode_table, barcode_table_2, barcode, lines, lines_2):
    if barcode not in barcode_table:
        fpath = os.path.join(output_folder, barcode + '_1_.fastq')
        barcode_table[barcode] = open(fpath, 'w',1048576)
    if barcode not in barcode_table_2:
        fpath_2 = os.path.join(output_folder, barcode + '_2_.fastq')
        barcode_table_2[barcode] = open(fpath_2, 'w',1048576) 

    outfile = barcode_table[barcode]
    for line in lines:
        outfile.write(line)

    outfile = barcode_table_2[barcode]
    for line in lines_2:
        outfile.write(line)

def barcode_buckets(paths):
    # hash table barcode -> file descriptor
    print(datetime.datetime.now())
    barcode_table = {}
    barcode_table_2 = {}
    try:
        os.mkdir(output_folder)
    except: 
        print('Could not create output folder')
        sys.exit(1)

    barcode_re = re.compile('barcode=([\S]+)')
    processed_file = open('./processed.txt','at+',1)
    processed_entries = processed_file.readlines()
    for path in paths:
        for dirpath, dirnames, filenames in os.walk(path):
            for filename in filenames:
                if filename in processed_entries:
                    continue

                print("Input file {}".format(filename))
                extension = os.path.splitext(filename)[1]
                if(extension != '.fastq_bc'):
                    continue
                print(extension)
                full_path = os.path.join(dirpath, filename)
                full_path_2 = full_path + '_pair'
                with open(full_path) as fastqfile:
                    with open(full_path_2) as fastqfile_2:
                        lines = fastqfile.readlines()                 #Read full file in order to work through with linenumbers to be able to write blocks of 4
                        lines_2 = fastqfile_2.readlines()
                        if len(lines) != len(lines_2):
                            print('Not same file length')
                            sys.exit(1)
                        for line_nr in range(0,(len(lines)-1)):
                            match = barcode_re.search(lines[line_nr]) #changed to search instead of match since not at beginning
                            if match == None:
                                continue
                            barcode = match.group(1)
                            write_line(barcode_table, barcode_table_2, barcode, lines[(line_nr-2):(line_nr+2)], lines_2[(line_nr-2):(line_nr+2)]) #Pass on a block of 4 lines (2 before and 1 after match)
                        processed_file.write(filename + '\n')
    for barcode in barcode_table:
        fd = barcode_table[barcode]
        fd.close()
    for barcode in barcode_table_2:
        fd_2 = barcode_table_2[barcode]
        fd_2.close()
    processed_file.close()
    print(datetime.datetime.now())   


if sys.argv[1:]:
    barcode_buckets(sys.argv[1:])
    sys.exit(0) 
print("Usage: python script.py inputpath [paths...]")
