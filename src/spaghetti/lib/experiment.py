import urllib
from os import path
from config import fastq_target_dir
import hashlib
import gzip
import shutil


class Experiment:
    def __init__(self, run_accession, fastq_ftp, fastq_md5):
        print self.run_accession, 'INIT EXPERIMENT'
        self.run_accession = run_accession
        self.fastq_ftp = fastq_ftp
        self.fastq_md5 = fastq_md5

        # extract filename
        parts = fastq_ftp.split('/')
        self.filename = parts[-1]

        # do the experiment: download, convert, kraken, brack
        # at each step, record whether it succeeded or not
        self.fastq_dl_success = self.do_fastq_dl()

        if self.fastq_dl_success:
            self.cortex_conversion_success = self.do_cortex_conversion()
            self.kraken_file_success = self.do_kraken_file()

            if self.kraken_file_success:
                self.brack_file_success = self.do_brack_file()
            else:
                self.brack_file_success = False
        else:
            self.cortex_conversion_success = False
            self.kraken_file_success = False
            self.brack_file_success = False

        summary = [self.fastq_dl_success, self.cortex_conversion_success, self.kraken_file_success,
                   self.brack_file_success]
        print self.run_accession, 'COMPLETE EXPERIMENT', summary

    def do_fastq_dl(self):
        try:
            print self.run_accession, 'FETCH START'

            # 1 check whether file exists in the target directory
            if not path.isfile(path.join(fastq_target_dir, self.filename)):
                # 2 download and write the file to the target directory with the extracted filename
                urllib.urlretrieve(self.fastq_ftp, path.join(fastq_target_dir, self.filename))

            # 3 verify self.fastq_md5 matches md5 of the downloaded file
            if not self.md5_match():
                print self.run_accession, 'MD5 CHECK FAILED'
                return False

            # 4 extract the fastq file from gzip - could delete the gz file at this point
            with gzip.open(path.join(fastq_target_dir, self.filename), 'rb') as f_in:
                uncompressed_name = self.filename[:-2]
                with open(path.join(fastq_target_dir, uncompressed_name), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # 5 return True unless there was an error at any point above
            return True
        except Exception as e:
            print(e)
            return False
        finally:
            print self.run_accession, 'FETCH END'

    # TODO: Ed
    def do_cortex_conversion(self):
        try:
            print self.run_accession, 'CONVERT START'
            # 1 do the cortex conversion logic
            # 2 return True if success, False otherwise
            return True
        except Exception as e:
            print(e)
            return False
        finally:
            print self.run_accession, 'CONVERT END'

    # TODO: Ed
    def do_kraken_file(self):
        try:
            print self.run_accession, 'KRAKEN START'
            # 1 do the kraken file logic
            # 2 return True if success, False otherwise
            # since brack depends on kraken, I guess the brack file step needs some information about the
            # kraken step so you may want to store some information in the instance
            # i.e. self.some_kraken_data = 'whatever' so that the do_brack_file func can access it
            # look at self.filename for an example of using an instance variable
            return True
        except Exception as e:
            print(e)
            return False
        finally:
            print self.run_accession, 'KRAKEN END'

    # TODO: Ed
    def do_brack_file(self):
        try:
            print self.run_accession, 'BRACK START'
            # 1 do the brack file logic
            # 2 return True if success, False otherwise
            return True
        except Exception as e:
            print(e)
            return False
        finally:
            print self.run_accession, 'BRACK END'

    def md5_match(self):
        hash_md5 = hashlib.md5()
        with open(path.join(fastq_target_dir, self.filename), "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest() == self.fastq_md5
