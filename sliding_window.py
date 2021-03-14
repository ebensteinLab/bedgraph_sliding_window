import os

def GetChrLengthsLists(chromsizes_path):
    chr_names=[]
    chr_lengths=[]
    with open(chromsizes_path) as hg38_chrsize:
            for chromline in hg38_chrsize:
                chrom_name,chrom_size=chromline.split("\t")
                chrom_size=chrom_size.strip()
                chr_names.append(chrom_name)
                chr_lengths.append(int(chrom_size))
    return (chr_names,chr_lengths)

def GetLine(line):
    line_chr,pos_start,pos_end,bg_value=line.split("\t")
    pos_start=int(pos_start)
    pos_end=int(pos_end)
    bg_value=float(bg_value)
    return (line_chr,pos_start,pos_end,bg_value)

def getBinVals (line, positions, values, BGfile, current_chr, bin_end,bin_sum,bin_start):
    while line:
        line_chr,pos_start,pos_end,bg_value=GetLine(line)
        if line_chr==current_chr and pos_start>= bin_start:
            if pos_end<=bin_end:
                for i in range(pos_start,pos_end):
                    positions.append(i)
                    values.append(bg_value)
                    bin_sum+=bg_value
                line=BGfile.readline()
            else:
                for i in range(pos_start,bin_end+1):
                    positions.append(i)
                    values.append(bg_value)
                    bin_sum+=bg_value
                pos_start=bin_end+1
                line="\t".join([line_chr,str(pos_start),str(pos_end),str(bg_value)])+"\n"
                break
        else:
            break
    return (line, bin_sum)  


chromsizes_path=r'/scratch200/denalaadan/hg38_chromSize.txt'
chr_names,chr_lengths=GetChrLengthsLists(chromsizes_path)

input_bg_path=r'/scratch200/denalaadan/sample1_sort_cov_bl.bedgraph'
output_file_path=r'/scratch200/denalaadan/aoa/s1_s200aoa.bedgraph'
bin_size=200

with open(input_bg_path) as BGfile:
    prev_chrom = ""
    with open(output_file_path, 'w') as output_file:
        line=BGfile.readline()
        while line:
            line_chr,pos_start,pos_end,bg_value=GetLine(line)
            if line_chr not in chr_names or pos_end > chr_lengths[chr_names.index(line_chr)]: #if the chromosome is not in the chrom_size file, or the position is beyond the chrom legth, disregard.
                line = BGfile.readline()
                continue
            if line_chr != prev_chrom: 
                print("reading",line_chr)
                output_file.flush()
                prev_chrom = line_chr
                positions=[]
                values=[]
                bin_start=0
                bin_sum = 0
                write_start=0
                write_end=1
                write_avg=None
            for bin_end in range(bin_size-1, chr_lengths[chr_names.index(prev_chrom)]+1): 
            
                print(bin_start,bin_end)
                line,bin_sum=getBinVals(line, positions, values, BGfile, prev_chrom, bin_end,bin_sum,bin_start)
                if positions and positions[0]<bin_start: 
                    bin_sum -= values[0]
                    del positions[0]
                    del values[0]
                
                avg=bin_sum/bin_size
                midpoint=(bin_end+bin_start)//2
                bin_start += 1
                new_chr=prev_chrom
                new_start=midpoint 
                new_end=midpoint+1
                new_avg=avg                    
                if write_avg==None or write_avg==0: 
                    write_start=new_start
                    write_end=new_end
                    write_avg=new_avg
                elif new_avg==write_avg: 
                    write_end=new_end
                elif write_avg!=0: 
                    outputline= "\t".join([prev_chrom, str(write_start), str(write_end), format(write_avg,'.3f')])+"\n"
                    output_file.write(outputline)
                    
                    write_start=new_start
                    write_end=new_end
                    write_avg=new_avg
            
            if write_avg!=0:
                outputline= "\t".join([prev_chrom, str(write_start), str(write_end), format(write_avg,'.3f')])+"\n"
                output_file.write(outputline)
                

                
