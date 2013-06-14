#!/usr/bin/python


INPUTFILE = './data/normal_tissue.csv'
OUTFILE_FOLDER = './data/tissues/'

example_line = '"ENSG00000000003","appendix","glandular cells","Moderate","Staining","Supportive"'

# counts number of lines in a file
#  source: http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
#   where it is defined as "bufcount"
def linesInFile(filename):
    f = open(filename)                  
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization
    
    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)

    return lines

def cleanString(field):
    return field.replace(" ",".").replace('"','').replace("/",".")

def getDescriptor(line):
    fields = line.split(",")
    tissue = cleanString(fields[1])
    celltype = cleanString(fields[2])
    return tissue + "_" + celltype

def splitIntoTissues(filename):
    print "Getting number of lines ..."
    lineno = linesInFile(filename)
    print "Number of lines is = " + str(lineno)
    lines_processed = 0
    last_percent = 0
    
    with open(filename, 'r') as infile:
        # read header line of input file
        headerline = infile.readline()
        lines_processed += 1
        
        # initialize open-file-dict
        out_files = {}
        
        # iterate through all lines in the file
        for line in iter(infile.readline, b''):
            lines_processed += 1
            descriptor = getDescriptor(line)
            outfilename = OUTFILE_FOLDER + descriptor + ".csv"
            
            # check if the file is already opened
            if (not out_files.has_key(descriptor)):
                # if not, first write the header
                outfile = open(outfilename, 'w')
                outfile.write(headerline)
                # the the new line
                outfile.write(line)
                # and add the open file to the dict
                out_files[descriptor] = outfile
                print "Creating new output file:", outfilename
            else:
                # just write the line to the already opened file
                out_files[descriptor].write(line)
                
            # update progress counter
            new_percent = 100*lines_processed/lineno
            if (new_percent >= last_percent + 10):
                print "Processing file ... " + str(new_percent) + "% done"
                last_percent = last_percent + 10
                
    # close all open files
    for file in out_files.values():
        file.close()
        
def main():
    splitIntoTissues(INPUTFILE)
    
main()