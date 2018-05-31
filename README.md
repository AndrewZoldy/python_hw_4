# python_hw_4

Given script builds De Brujin graph from reads in your .fasta file.

##Arguments

###Required:

Arguments which are required for launching script:
  * **-i (--input)** - path to your .fasta formatted reads
    * *-i /path/to/your/fasta_file.fasta*
    * *--input /path/to/your/fasta_file.fasta*
  * **-o (--output)** - path to output .dot file contains graph
    * *-o /path/to/output_file.dot*
    * *--output /path/to/output_file.dot*

###Optional:

Options which you can switch:
  * **-s (--size)** - size of k-mers (by default = 15)
    * *-s 5*
    * *--size 5*
  * **-m (--mode)** - graph can be shown in full mode (with names of nodes and edges) and in cutted mode (with only values). Use **'full'** or **'cutted'** values of this argument.
    * *-m 'full'*
    * *--mode 'full'*
    * *-m 'cutted'*
    * *--mode 'cutted'*
    
    
##Results

As a result you will obtain image of maden graph.
