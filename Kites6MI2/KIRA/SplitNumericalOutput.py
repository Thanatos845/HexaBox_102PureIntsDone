
import sys, getopt, os
def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print ('SplitNumericalOutput.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('SplitNumericalOutput.py -i <inputfile> -o <outputfile>')
            sys.exit()
        if opt in ("-i", "--ifile"):
            inputfile = arg
        if opt in ("-o", "--ofile"):
            outputfile = arg
    return (inputfile,outputfile)

inputfile = ''
outputfile = ''
if __name__ == "__main__":
   inputfile,outputfile=main(sys.argv[1:])


print ('Input file is ', 1000*int(inputfile[-3])+1)
# print ('Output file is ', outputfile)



 

os.system(f"mkdir -p {outputfile}/")

linecount=1000*int(inputfile[-3])+1

file1 = open(f"{outputfile}/phasespace_{linecount:05d}.m", 'w')
fileR = open(inputfile, 'r') 


prevline=""

linecount=linecount-1

for i,line in enumerate(fileR):
    if i==0:
        primeline=line
        file1.writelines(primeline)

    if (not "{hxb" in line) and ("{{" in line):
        linecount+=1
        # print(linecount)
        if (linecount>1):
            file1.writelines(prevline[:-2]+"}")
            file1.close()
            file1 = open(f'{outputfile}/phasespace_{linecount:05d}.m', 'w')
            file1.writelines(primeline)
    elif i>0:
        if "}}}}," in prevline:
            prevline=prevline[:-2]
        file1.writelines(prevline)
    prevline=line

file1.writelines(prevline)
file1.close()