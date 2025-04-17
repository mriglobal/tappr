import os
import zipfile
import argparse


parser = argparse.ArgumentParser(description="Subtract identical archive members in zip archive B from zip archive A.")

subparsers = parser.add_subparsers(dest='subparser_name',help='Select from Union or Difference routines.')
parser_sub = subparsers.add_parser('subtract',help='Subtract identical archive members in zip archive B from zip archive A.')
parser_sub.add_argument('-a',required=True,type=str,help='Zip archive A to be filtered.')
parser_sub.add_argument('-b',required=True,type=str,help='Zip archive B whose members will be filtered from A.')
parser_sub.add_argument('-o',default=None,help='Output archive name. Default is to inherit names from A and B archive file names.')
parser_add = subparsers.add_parser('add',help='Aggregate archive members into a single zip.')
parser_add.add_argument('-z',required=True,nargs='+',help='List of zip archives to aggregate.')
parser_add.add_argument('-o',default=None,help='Output archive name. Default is to inherit name from first archive in list.')

myargs = parser.parse_args()

if myargs.subparser_name == 'subtract':
    if myargs.o:
        outname = myargs.o
    else:
        outname = os.path.basename(myargs.a).replace('.zip','')+"_less_"+os.path.basename(myargs.b)

    with zipfile.ZipFile(myargs.b,'r') as file1:
        file1_set = set(file1.namelist())
        with zipfile.ZipFile(myargs.a,'r') as file2:
            output = [s for s in file2.namelist() if s not in file1_set]
            with zipfile.ZipFile(outname,'w') as outfile:
                for o in output:
                    print("Writing {} to zip.".format(o))
                    outfile.writestr(o,file2.read(o),compress_type=zipfile.ZIP_DEFLATED,compresslevel=6)

elif myargs.subparser_name == 'add':
    if myargs.o:
        outname = myargs.o
    else:
        outname = os.path.basename(myargs.z[0]).replace('.zip','')+"_aggregated.zip"
    
    with zipfile.ZipFile(outname,'w') as outfile:
        for file in myargs.z:
            print(file)
            with zipfile.ZipFile(file,'r') as read_file:
                for n in read_file.namelist():
                    print("Writing {} to zip.".format(n))
                    outfile.writestr(n,read_file.read(n),compress_type=zipfile.ZIP_DEFLATED,compresslevel=6)
