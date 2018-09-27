# Import libraries
import hashlib
import os.path
import os
import subprocess
import string
from subprocess import call

# Declare variables - bacterial, archaeal & viral reads text files respectively. Contains FTP links & zipped fastq md5s.
bacReads = "./bacterial_Reads.txt"
arcReads = "./archaeal_Reads.txt"
virReads = "./viral_Reads.txt"
retCheck = "./run_Table_Check.txt"

# Make a list of these files for later directory creation
readFiles = [bacReads, arcReads, virReads]


# Define methods

# Check for all necessary files for each read
# def retCheck()

# Return checksum of file specified
def md5Checksum(filePath):
    with open(filePath, 'rb') as fh:
        m = hashlib.md5()
        while True:
            data = fh.read(8192)
            if not data:
                break
            m.update(data)
        return m.hexdigest()


# Make directiories if they don't exist
# Need to add conditionals for each filetype (cortex, bloom, bracken results) & md5 checks
def dirExists(inputDirList):
    for file in inputDirList:
        dir = os.path.dirname(file)
        if not os.path.exists(dir):
            os.makedirs(dir)
            return True
        # If filepath exists for read, check that a cortex file is present
        elif os.path.exists(dir):
            return False


# Kraken2 classification
def krakenClass(fastqExtLst):
    for ext in fastqExtLst:
        fileSubStr = ext[30:-8]
        if "_1" in ext:
            continue
        elif "_2" in ext:
            try:
                print(subprocess.check_output(
                    "/homes/richards/.local/bin/kraken2 --db /hps/nobackup/iqbal/phelim/projects/bigsi-bench-richards/krak2DB --threads 10 --report " + ext[
                                                                                                                                                        0:-8] + ".kreport --paired " + ext[
                                                                                                                                                                                       0:-8] + "_1.fastq " + ext[
                                                                                                                                                                                                             0:-8] + "_2.fastq > " + ext[
                                                                                                                                                                                                                                     0:-8] + ".kraken",
                    shell=True))
                with open(retCheck, "a") as retFile:
                    retFile.write("Y\t")
            except subprocess.CalledProcessError, e:
                print("Kraken did not run successfully, updated table for later processing.")
                with open(retCheck, "a") as retFile:
                    retFile.write("N\t")
            try:
                print(subprocess.check_output(
                    "/homes/richards/Bracken-master/bracken -d /hps/nobackup/iqbal/phelim/projects/bigsi-bench-richards/krak2DB -i " + ext[
                                                                                                                                       0:-8] + ".kreport -o " + ext[
                                                                                                                                                                0:-8] + ".bracken -r 100 -l G -t 10",
                    shell=True))
                with open(retCheck, "a") as retFile:
                    retFile.write("Y\n")
            except subprocess.CalledProcessError, e:
                print("Bracken did not run successfully, updated table for later processing.")
                with open(retCheck, "a") as retFile:
                    retFile.write("N\n")
        else:
            fileSubStr = ext[0:-6]
            try:
                print(subprocess.check_output(
                    "/homes/richards/.local/bin/kraken2 --db /hps/nobackup/iqbal/phelim/projects/bigsi-bench-richards/krak2DB --threads 10 --report " + ext[
                                                                                                                                                        0:-6] + ".kreport " + ext + " >" + fileSubStr + ".kraken",
                    shell=True))
                with open(retCheck, "a") as retFile:
                    retFile.write("Y\t")
            except subprocess.CalledProcessError, e:
                print("Kraken did not run successfully, updated table for later processing.")
                with open(retCheck, "a") as retFile:
                    retFile.write("N\t")
            try:
                print(subprocess.check_output(
                    "/homes/richards/Bracken-master/bracken -d /hps/nobackup/iqbal/phelim/projects/bigsi-bench-richards/krak2DB -i " + ext[
                                                                                                                                       0:-6] + ".kreport -o " + ext[
                                                                                                                                                                0:-6] + ".bracken -r 100 -l G -t 10",
                    shell=True))
                with open(retCheck, "a") as retFile:
                    retFile.write("Y\n")
            except subprocess.CalledProcessError, e:
                print("Bracken did not run successfully, updated table for later processing.")
                with open(retCheck, "a") as retFile:
                    retFile.write("N\n")


# Download read(s) when no. of entries for FTP column in text file is 7
def pairReadDown(readsListSplit, fastqDir, fastqFile, fastqFile1):
    print("Retrieving paired-end reads...\n\n")
    try:
        print(subprocess.check_output("wget -q " + readsListSplit[2] + " -O " + fastqDir + fastqFile, shell=True))
    except subprocess.CalledProcessError, e:
        print("First paired read not downloaded successfully, updating table")
    try:
        print(subprocess.check_output("wget " + readsListSplit[3] + " -O " + fastqDir + fastqFile1, shell=True))
    except subprocess.CalledProcessError, e:
        print("Second paired read not downloaded successfully, updating table")


# Download read(s) when no. of entries for FTP column in text file is 5
def pairReadDown1(readsListSplit, fastqDir, fastqFile, fastqFile1):
    print("Retrieving paired-end reads...\n\n")
    subprocess.call("wget -q " + readsListSplit[1] + " -O " + fastqDir + fastqFile, shell=True)
    subprocess.call("wget -q " + readsListSplit[2] + " -O " + fastqDir + fastqFile1, shell=True)


# Download single reads
def singReadDown(readsListSplit, fastqDir, fastqFile):
    print("Retrieving single read...\n\n")
    subprocess.call("wget -q " + readsListSplit[1] + " -O " + fastqDir + fastqFile, shell=True)


# Method for paired read cortex conversion
def cortexConv(fileList):
    for ext in fileList:
        if "_1" in ext:
            continue
        elif "_2" in ext:
            ctxFile = ext[0:-8] + ".ctx"
            fileSubStr = ext[30:-8]
            try:
                print(subprocess.check_output(
                    "/homes/richards/mccortex/bin/mccortex31 build -f -q -t 10 -m 15G -k 31 -s fastq -2 " + ext[
                                                                                                            0:-8] + "_1.fastq:" + ext[
                                                                                                                                  0:-8] + "_2.fastq " + ext[
                                                                                                                                                        0:-8] + ".ctx",
                    shell=True))
                with open(retCheck, "a") as retFile:
                    retFile.write("Y\t")
            except subprocess.CalledProcessError, e:
                print(ctxFile + "was not generated successfully! Adding to list for later processing")
                with open(retCheck, "a") as retFile:
                    retFile.write("\tN\t")
        else:
            ctxFile = ext[0:-6] + ".ctx"
            fileSubStr = ext[25:-6]
            try:
                print(subprocess.check_output(
                    "/homes/richards/mccortex/bin/mccortex31 build -f -q -t 10 -m 15G -k 31 -s fastq -1 " + ext + " " + ext[
                                                                                                                        0:-6] + ".ctx",
                    shell=True))
                with open(retCheck, "a") as retFile:
                    retFile.write("Y\t")
            except subprocess.CalledProcessError, e:
                print(ctxFile + "was not generated successfully! Adding to list for later processing")
                with open(retCheck, "a") as retFile:
                    retFile.write("\tN\t")


# Method for paired-end download
# Reads each line in text file into a string. Converts this string into a list which is split based on the delimiters
# Substring the elements of the index in order to get the directory in which to store later files & fast.gz's
# Continually attempts to download the fastq files - this will be changed later to 'n' attempts
# Returns boolean in order to add to the conversion counter (ctx, bloom, kraken/bracken)
# If a read exists already, won't add to conversion counter
def pairedRead(readLine):
    readLine = readLine.replace("\t", ";")
    readSplitList = readLine.split(";")
    with open(retCheck, "a") as retFile:
        retFile.write(readSplitList[0] + "\t")
    # This occurs when paired reads are present without single fastq ftp link (i.e: _1.fastq.gz & _2.fastq.gz, but no .fastq.gz
    if len(readSplitList) == 7:
        if (len(readSplitList[0]) == 10):
            # The directory we're going to create...
            createDir = "./" + readSplitList[1][23:51]
            # The checksums for each paired file
            checkSumGz = readSplitList[4]
            checkSumGz1 = readSplitList[5]
            dirList = [createDir]
            # If the directory exists or not
            if dirExists(dirList):
                print("New reads! Creating directory " + createDir + "\n\n")
                # Open table file to add check to filepath exists
                # The paired reads without absolute path
                readFile = readSplitList[1][51:72]
                readFile1 = readSplitList[2][51:72]
                # download the paired reads
                pairReadDown1(readSplitList, createDir, readFile, readFile1)
                # If the reads were not downloaded successfully
                while not fastqCheckPaired(createDir, readFile, readFile1, checkSumGz, checkSumGz1):
                    try:
                        print("Did not download paired reads successfully. Retrying...\n\n")
                        pairReadDown1(readSplitList, createDir, readFile, readFile1)
                    except:
                        print("Downloaded paired reads successfully after retrying\n\n")
                return True
            else:
                # If the reads exist already, skip this entry
                print(readSplitList[0] + " already exists. Skipping...\n\n")
                return False
        # Method depends on read string length i.e.: DRR006497 is length 9 whereas SRR7758677 is length 10)
        elif len(readSplitList[0]) == 9:
            createDir = "./" + readSplitList[1][23:46]
            checkSumGz = readSplitList[4]
            checkSumGz1 = readSplitList[5]
            dirList = [createDir]
            if dirExists(dirList):
                print("New reads! Creating directory " + createDir + "\n\n")
                readFile = readSplitList[1][46:68]
                readFile1 = readSplitList[2][46:68]
                pairReadDown1(readSplitList, createDir, readFile, readFile1)
                while not fastqCheckPaired(createDir, readFile, readFile1, checkSumGz, checkSumGz1):
                    try:
                        print("Did not download paired reads successfully. Retrying...\n\n")
                        pairReadDown1(readSplitList, createDir, readFile, readFile1)
                    except:
                        print("Downloaded paired reads successfully after retrying\n\n")
                return True
            else:
                print(readSplitList[0] + " already exists. Skipping...\n\n")
                return False
    # Return different elements of the list if list is different length (i.e.: sometimes more row entries than others in the text files for columns)
    elif len(readSplitList) == 9:
        if (len(readSplitList[0]) == 9):
            createDir = "./" + readSplitList[1][23:46]
            checkSumGz = readSplitList[6]
            checkSumGz1 = readSplitList[7]
            dirList = [createDir]
            if dirExists(dirList):
                print("New reads! Creating directory " + createDir + "\n\n")
                readFile = readSplitList[2][51:69]
                readFile1 = readSplitList[3][51:69]
                pairReadDown(readSplitList, createDir, readFile, readFile1)
                while not fastqCheckPaired(createDir, readFile, readFile1, checkSumGz, checkSumGz1):
                    try:
                        print("Did not download paired reads successfully. Retrying...\n\n")
                        pairReadDown(readSplitList, createDir, readFile, readFile1)
                    except:
                        print("Downloaded paired reads successfully after retrying\n\n")
                return True
            else:
                print(readSplitList[0] + " already exists. Skipping...\n\n")
                return False
        elif len(readSplitList[0]) == 10:
            createDir = "./" + readSplitList[1][23:51]
            checkSumGz = readSplitList[6]
            checkSumGz1 = readSplitList[7]
            dirList = [createDir]
            if dirExists(dirList):
                print("New reads! Creating directory " + createDir + "\n\n")
                readFile = readSplitList[2][51:69]
                readFile1 = readSplitList[3][51:69]
                pairReadDown(readSplitList, createDir, readFile, readFile1)
                while not fastqCheckPaired(createDir, readFile, readFile1, checkSumGz, checkSumGz1):
                    try:
                        print("Did not download paired reads successfully. Retrying...\n\n")
                        pairReadDown(readSplitList, createDir, readFile, readFile1)
                    except:
                        print("Downloaded paired reads successfully after retrying\n\n")
                return True
            else:
                print(readSplitList[0] + " already exists. Skipping...\n\n")
                return False
    else:
        print(readSplistList)


# Method for single read download
# Reads each line in text file into a string. Converts this string into a list which is split based on the delimiters
# Substring the elements of the index in order to get the directory in which to store later files & fast.gz's
# Continually attempts to download the fastq files - this will be changed later to 'n' attempts
# Returns boolean in order to add to the conversion counter (ctx, bloom, kraken/bracken)
# If a read exists already, won't add to conversion counter
def singRead(readLine):
    readLine = readLine.replace("\t", ";")
    readSplitList = readLine.split(";")
    with open(retCheck, "a") as retFile:
        retFile.write(readSplitList[0] + "\t")
    if len(readSplitList[0]) == 9:
        createDir = "./" + readSplitList[1][23:46]
        checkSumGz = readSplitList[2]
        dirList = [createDir]
        if dirExists(dirList):
            print("New read! Creating directory " + createDir + "\n\n")
            readFile = readSplitList[1][46:68]
            singReadDown(readSplitList, createDir, readFile)
            while not fastqCheckSing(createDir, readFile, checkSumGz):
                try:
                    print("Did not download single read successfully. Retrying...\n\n")
                    singReadDown(readSplitList, createDir, readFile)
                except:
                    print("Single read downloaded successfully\n\n")
            return True
        else:
            print(readSplitList[0] + " already exists. Skipping...\n\n")
            return False
    elif len(readSplitList[0]) == 10:
        createDir = "./" + readSplitList[1][23:51]
        checkSumGz = readSplitList[2]
        dirList = [createDir]
        if dirExists(dirList):
            print("New reads! Creating directory " + createDir + "\n\n")
            readFile = readSplitList[1][51:69]
            singReadDown(readSplitList, createDir, readFile)
            while not fastqCheckSing(createDir, readFile, checkSumGz):
                try:
                    print("Did not download single read successfully. Retrying...\n\n")
                    singReadDown(readSplitList, createDir, readFile, readFile1)
                except:
                    print("Single read downloaded successfully after retrying\n\n")
            return True
        else:
            print(readSplitList[0] + " already exists. Skipping...\n\n")
            return False


# Method to run through for each read
# Extracts all compressed fastqs, converts & classifies them
def readProcess(typeRead):
    typeRead = typeRead.replace("\t", ";")
    typeReadSplitList = typeRead.split(";")
    if len(typeReadSplitList[0]) == 9:
        createDir = "./" + typeReadSplitList[1][23:46]
    elif len(readSplitList[0]) == 10:
        createDir = "./" + readSplitList[1][23:51]
    if "_2" in typeRead:
        subprocess.call("gunzip -f " + createDir + createDir[-11:len(createDir) - 1] + "_1.fastq.gz", shell=True)
        subprocess.call("gunzip -f " + createDir + createDir[-11:len(createDir) - 1] + "_2.fastq.gz", shell=True)
    elif "_2" not in typeRead:
        subprocess.call("gunzip -f " + createDir + createDir[-10:len(createDir) - 1] + ".fastq.gz", shell=True)
    # Get list of fastq files for conversion & classification
    proc1 = subprocess.Popen("find " + createDir + " -type f -name '*.fastq'",
                             stdout=subprocess.PIPE, shell=True)
    fastqExt = proc1.stdout.read()
    fastqExtLst = fastqExt.split()
    print(fastqExtLst)
    # Iterate over list & kraken/bracken, mccortex31, bigsi bloom & cobs cortex all gunzipped fastqs
    print("Converting to cortex files...")
    cortexConv(fastqExtLst)
    krakenClass(fastqExtLst)


#    for ext in fastqExtLst:
#        subprocess.call("rm " + ext, shell=True)

# Method for reading in each text file (bacterial, viral & arachaeal)
# Perform if & else statement to accurately return the relevant read files for processing
def readFile(typeReads):
    krakConvCounter = 0
    with open(typeReads) as b:
        typeLines = [line.rstrip('\n') for line in b]
    for typeRead in typeLines:
        # Lines where no FTP URL is available
        if "gz" not in typeRead:
            continue
        # Lines which contain 2 read files (paired)
        elif "_2" in typeRead:
            pairedRead(typeRead)
            readProcess(typeRead)
        # Lines which contain only a single read
        elif "_2" not in typeRead:
            singRead(typeRead)
            readProcess(typeRead)


# Method to check whether fastq was retrieved successfully for paired reads
def fastqCheckPaired(readDir, fastqFile, fastqFile1, checkSumGz, checkSumGz1):
    fastqString = readDir + fastqFile
    fastqString1 = readDir + fastqFile1
    if (os.path.exists(readDir + fastqFile) and os.path.exists(readDir + fastqFile1) and md5Checksum(
            fastqString) == checkSumGz and md5Checksum(fastqString1) == checkSumGz1):
        print("Paired reads downloaded successfully\n\n")
        with open(retCheck, "a") as retFile:
            retFile.write("\tY\t")
        return True
    else:
        with open(retCheck, "a") as retFile:
            retFile.write("\tN\t")
        return False


# Method to check whether fastq was retrieve successfuly for single read(s)
def fastqCheckSing(readDir, fastqFile, checkSumGz):
    fastqString = readDir + fastqFile
    if (os.path.exists(readDir + fastqFile) and md5Checksum(fastqString) == checkSumGz):
        print("Single read downloaded successfully\n\n")
        with open(retCheck, "a") as retFile:
            retFile.write("\tY\t")
        return True
    else:
        with open(retCheck, "a") as retFile:
            retFile.write("\tN\t")
        return False


# Open & write columns in table file for checks
if os.path.isfile(retCheck):
    print("Table file already exists")
    with open(retCheck, "a") as retFile:
        retFile.write("\n")
else:
    with open(retCheck, "a") as retFile:
        retFile.write("READ\tFASTQ(S)\tCTXFILE\tKRAKFILE\tBRACKFILE\n")

# Check if txts containing bacterial, archaeal & viral reads exists
# If doesn't exist, download
# If exists, continue
for file in readFiles:
    if os.path.isfile(file):
        continue
    else:
        # Download reads
        print("Downloading bacterial reads text file...")
        subprocess.call(
            "curl -sA \"Chrome\" -L 'https://www.ebi.ac.uk/ena/data/warehouse/search?query=%22tax_tree(2)%20AND%20first_created%3E=2017-01-01%22&limit=358768&length=358768&offset=1&display=report&result=read_run&fields=fastq_ftp,fastq_md5&download=txt' -o bacterial_Reads.txt",
            shell=True)
        print("Downloading archaeal reads text file...")
        subprocess.call(
            "curl -sA \"Chrome\" -L 'https://www.ebi.ac.uk/ena/data/warehouse/search?query=%22tax_tree(2157)%20AND%20first_created%3E=2017-01-01%22&limit=0&length=0&offset=1&display=report&result=read_run&fields=fastq_ftp,fastq_md5&download=txt' -o archaeal_Reads.txt",
            shell=True)
        print("Downloading viral reads text file...")
        subprocess.call(
            "curl -sA \"Chrome\" -L 'https://www.ebi.ac.uk/ena/data/warehouse/search?query=%22tax_tree(10239)%20AND%20first_created%3E=2017-01-01%22&limit=0&length=0&offset=1&display=report&result=read_run&fields=fastq_ftp,fastq_md5&download=txt' -o viral_Reads.txt",
            shell=True)
        print("Done")

# Work through bacterial reads
# Open downloaded bacterial text file & assign each line to a list
# print("Downloading bacterial reads...")
# readFile(bacReads)

# Work through viral reads
# Open downloaded viral text file & assign each line to a list
print("Downloading viral reads...")
readFile(virReads)

# Work through archaeal reads
# Open downloaded archaeal text file & assign each line to a list
# print("Downloading archaeal reads...")
# readFile(arcReads)
