#!/usr/bin/env python3
import urllib.request
import time
import os.path
import re


cur_folder = os.path.dirname(__file__)
REGISTRY_FILE = os.path.join(cur_folder, 'PSICQUIC_registry_imex.txt')
PPI_FOLDER = os.path.join(cur_folder, 'ppis')


def get_psicquic_services(filename):
    """
    Returns all PSICQUIC services from the registry file.
    TODO: extend this to also work with registry URLs.
    """
    services = {}
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            service, url = line.split("=")
            services[service] = url
    return services


def urlopen_retry(url, retry = 3, delay = 1):
    print("opening " + url)
    i = 0
    while (i < retry):
        try:
            cf = urllib.request.urlopen(url)
            return cf
        except urllib.error.HTTPError as e:
            # try again
            print("Failed to connect, retrying...")
            time.sleep(delay)
            i = i+1
            if (i == retry):
                return e


def download_all_ppis(folder, basename):
    """
    Downloads all PPIs provided by the PSICQUIC registry
    into the folder with the given basename as tsv file.
    E.g. for the PPI `MINT` the file name will be `folder/basename_MINT.tsv`.
    """
    # build PSICQUIC RESTFUL query (added to the base service url)
    # https://code.google.com/p/psicquic/wiki/PsicquicSpec_1_3_Rest
    # can be any of 'v1.0','v1.1','v1.2','v1.3' and `current`

    psicquic_version = 'current'
    psicquic_method = 'query'
    psicquic_query = 'taxidA:9606%20AND%20taxidB:9606'
    psicquic_parameters = 'format=tab25'   # use the PSI-MITAB 2.5 format

    psicquic_query_url = (psicquic_version + '/search/' + psicquic_method
                          + '/' + psicquic_query)

    # get dict of PSICQUIC services and service urls
    psiquic_services = get_psicquic_services(REGISTRY_FILE)
    for service_name, service_url in psiquic_services.items():
        # remove the last /psicquic from the service_url
        service_url = re.sub("/psicquic$", "", service_url)
        # first get row count
        count_url = service_url + '/' + psicquic_query_url + '?format=count'
        tab25_url = service_url + '/' + psicquic_query_url + '?format=tab25'

        # try opening URL
        #with urllib.request.urlopen(count_url) as cf:
        #    count = int(cf.read())
        try:
            cf = urlopen_retry(count_url)
            count = int(cf.read())
            cf.close()
        except:
            print("[ERROR] Failed to download PPI " + service_name + ", skipping...")
            continue

        print("Downloading from '" + service_name + "': " + str(count)
              + " rows ...")
        # get filename for output
        outfilename = os.path.join(folder, basename + "_" + service_name
                                   + ".tsv")
        # in case the file already exists
        if os.path.exists(outfilename):
            # check that all interactions are there
            with open(outfilename, "r") as f:
                num_lines = sum(1 for line in f)
            if num_lines == count:
                # skip this specific PPI
                print("PPI file already exists, skipping download.")
                continue
            else:
                # remove the file
                os.remove(outfilename)

        CHUNK_SIZE = 4*1024
        # open url and write out, chunk-by-chunk
        with urllib.request.urlopen(tab25_url) as tab25_f:
            with open(outfilename, "wb") as outfile:
                while True:
                    chunk = tab25_f.read(CHUNK_SIZE)
                    if not chunk:
                        break
                    outfile.write(chunk)


download_all_ppis(PPI_FOLDER, "PSICQUIC")
