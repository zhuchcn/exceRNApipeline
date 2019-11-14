from ftplib import FTP
import re
import os
import argparse
import math
import progressbar
from magic import from_file
from pathos.multiprocessing import ProcessingPool as Pool
from utils import log


class EnsemblFTP():

    def __init__(self):
        self.connect()

    def connect(self):
        self.baseurl = "ftp.ensemblgenomes.org"
        ftp = FTP(self.baseurl)
        ftp.login("anonymous", "")
        self.ftp = ftp
    
    def ls(self, path):
        return self.ftp.nlst(path)
    
    def has_valid_domain(self, domain):
        if domain in ['bacteria', 'fungi', 'metazoa', 'plants', 'protists']:
            return True
        return False
    
    def download(self, path, filename, verbose):
        with open(filename, "wb") as fh:
        #log(f"Using url: {self.baseurl + path}")
            self.ftp.retrbinary("RETR " + path, fh.write)
        return True

    def getGenome(self, url, outdir, verbose):
        # check if path is valid, create diractory if not
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        #species = url.split('/')[-1]
        # find the fasta.gz file
        files= self.ls(f"{url}/dna")
        file_names = [os.path.basename(x) for x in files]
        r_matches = [
            bool(re.search("dna.toplevel.fa.gz$", x)) \
                for x in file_names
        ]
        path = files[r_matches.index(True)]
        
        filename = os.path.join(outdir, os.path.basename(path))
        self.download(path, filename, verbose)
        return True


def getBacteriaGenomes(outdir, version, threads, verbose=False):
    ftp = EnsemblFTP()
    collections = ftp.ls(f"pub/bacteria/{version}/fasta")

    def function(collection):
        if verbose:
            log('start downloading: ' + collection)
        ftp = EnsemblFTP()
        for species_url in ftp.ls(collection):
            output_dir = os.path.join(outdir, os.path.basename(collection))
            ftp.getGenome(species_url, output_dir, verbose)
        return
    
    with Pool(threads) as p:
        p.map(function, collections)

    return True